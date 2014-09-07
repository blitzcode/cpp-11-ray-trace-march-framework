
#include "framebuffer.h"
#include "trace.h"
#include "bmp_writer.h"
#include "timer.h"

#include <random>
#include <cassert>

Framebuffer::Framebuffer()
    : m_num_cpus(std::max(1u, std::thread::hardware_concurrency()))
{
    Trace("Initializing framebuffer with %i threads", m_num_cpus);
}

Framebuffer::~Framebuffer()
{
    KillAllWorkerThreads();
}

void Framebuffer::CreateWorkerThreads()
{
    assert(m_threads.empty());

    m_threads_done = 0;
    m_render_start_time = TimerGetTick();

    for (uint i=0; i<m_num_cpus; i++)
        m_threads.push_back(std::thread(&Framebuffer::WorkerThread, this));
}

void Framebuffer::KillAllWorkerThreads()
{
    // Shut down all worker threads, we're single threaded after the call returns

    assert(m_threads_stop == false);

    m_threads_stop = true;
    for (auto& thread : m_threads)
        if (thread.joinable())
            thread.join();
    m_threads.clear();
    m_threads_stop = false;
}

Framebuffer::Tile * Framebuffer::GetNextTileFromQueue()
{
    // Remove and return the next tile from the work queue.
    // Returns null if there are no more tiles left

    std::lock_guard<std::mutex> guard(m_work_queue_mtx);

    if (m_work_queue.empty())
        return nullptr;

    Tile *tile = &m_tiles[m_work_queue.back()];
    m_work_queue.pop_back();

    return tile;
}

void Framebuffer::WorkerThread()
{
    // Keep rendering tiles till we're done or asked to stop

    // Trace("Worker thread started");

    while (m_threads_stop == false)
    {
        Tile *tile = GetNextTileFromQueue();
        if (tile == nullptr)
            break;

        // Lock the tile while we call RenderTile() to work on it
        std::lock_guard<std::mutex> guard(tile->GetMutex());

        RenderTile(* tile);

        // Tile has been changed, texture needs to be updated
        tile->SetDirty(true);
    }

    // Increment thread done counter and check if we were the last
    if (std::atomic_fetch_add(&m_threads_done, 1u) == m_num_cpus - 1)
    {
        // Only report rendering time if we weren't asked to cancel
        if (m_threads_stop == false)
        {
            const double render_end_time = TimerGetTick() - m_render_start_time;
            Trace("Finished rendering after %.2fs", render_end_time);
        }
    }

    // Trace("Worker thread finished");
}

void Framebuffer::Resize(uint width, uint height)
{
    if (width == m_width && height == m_height)
        return;

    // Abandon current rendering efforts and be single threaded till we're done
    KillAllWorkerThreads();

    m_width  = width;
    m_height = height;

    // Compute new tile positions
    const uint tile_wdh = width  / m_tiles_x;
    const uint tile_hgt = height / m_tiles_y;
    for (uint y=0; y<m_tiles_y; y++)
        for (uint x=0; x<m_tiles_x; x++)
        {
            const uint tile_idx = x + y * m_tiles_x;

            m_tiles[tile_idx].SetPosition(
                x * tile_wdh,
                y * tile_hgt,
                (x == m_tiles_x - 1) ? width  : (x + 1) * tile_wdh,
                (y == m_tiles_y - 1) ? height : (y + 1) * tile_hgt);
        }

    FillWorkQueue();

    // Now we can render again
    CreateWorkerThreads();
}

void Framebuffer::RestartRendering()
{
    KillAllWorkerThreads();

    // Make sure all tiles and their textures have been cleared
    for (auto& tile : m_tiles)
    {
        tile.Clear();
        tile.UpdateTexture();
    }

    FillWorkQueue();

    CreateWorkerThreads();
}

void Framebuffer::FillWorkQueue()
{
    // Fill work queue with tiles
    m_work_queue.clear();
    for (uint i=0; i<m_tiles_x * m_tiles_y; i++)
        m_work_queue.push_back(i);

    // Random order
    std::random_shuffle(m_work_queue.begin(), m_work_queue.end());
}

void Framebuffer::Draw(uint x, uint y, uint width, uint height)
{
    glEnable(GL_TEXTURE_2D);

    for (uint ty=0; ty<m_tiles_y; ty++)
        for (uint tx=0; tx<m_tiles_x; tx++)
        {
            const uint tile_idx = tx + ty * m_tiles_x;
            Tile& cur_tile = m_tiles[tile_idx];

            // Try updating the texture with the buffer's content if the mutex is not
            // currently locked by a worker thread. UpdateTexture() checks and resets
            // the dirty flag automatically
            std::mutex& mtx = cur_tile.GetMutex();
            if (mtx.try_lock())
            {
                cur_tile.UpdateTexture();
                mtx.unlock();
            }

            // Tile position after shifting / scaling
            uint tx0, ty0, tx1, ty1;
            cur_tile.GetPosition(tx0, ty0, tx1, ty1);
            const float x0 = float(x) + float(tx0) * (float(width ) / float(m_width ));
            const float y0 = float(y) + float(ty0) * (float(height) / float(m_height));
            const float x1 = float(x) + float(tx1) * (float(width ) / float(m_width ));
            const float y1 = float(y) + float(ty1) * (float(height) / float(m_height));

            // Draw tile
            glBindTexture(GL_TEXTURE_2D, cur_tile.GetTexture());
            glColor3f(1.0f, 1.0f, 1.0f);
            glBegin(GL_QUADS);
                glTexCoord2f(0.0f, 0.0f);
                glVertex2f(x0, y0);
                glTexCoord2f(1.0f, 0.0f);
                glVertex2f(x1, y0);
                glTexCoord2f(1.0f, 1.0f);
                glVertex2f(x1, y1);
                glTexCoord2f(0.0f, 1.0f);
                glVertex2f(x0, y1);
            glEnd();
        }

    glDisable(GL_TEXTURE_2D);
}

void Framebuffer::SaveToBMP(const char *filename)
{
    std::vector<uint32> bitmap(m_width * m_height, 0);

    for (auto& tile : m_tiles)
    {
        // Tiles which are currently being rendered are just left black
        std::mutex& mtx = tile.GetMutex();
        if (mtx.try_lock())
        {
            uint x0, y0, x1, y1;
            tile.GetPosition(x0, y0, x1, y1);
            uint32 *buf = tile.GetBuffer();

            // Copy tile into bitmap
            for (uint y=0; y<tile.GetHeight(); y++)
                for (uint x=0; x<tile.GetWidth(); x++)
                    bitmap[x0 + x + (y0 + y) * m_width] = buf[x + y * tile.GetWidth()];

            mtx.unlock();
        }
    }

    WriteBitmap(filename, m_width, m_height, &bitmap[0]);

    Trace("Saved screenshot to '%s'", filename);
}

Framebuffer::Tile::Tile()
{
    SetPosition(0, 0, 1, 1);

    glGenTextures(1, &m_tex);
    glBindTexture(GL_TEXTURE_2D, m_tex);

    // No MIP-maps
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // This is important so we don't get artifacts at the borders from texture filtering
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

void Framebuffer::Tile::SetPosition(uint x0, uint y0, uint x1, uint y1)
{
    m_x0 = x0;
    m_y0 = y0;
    m_x1 = x1;
    m_y1 = y1;

    m_bgra.resize(GetWidth() * GetHeight());
    Clear();
    UpdateTexture();
}

void Framebuffer::Tile::Clear()
{
    // Clear image storage, flag for texture upload
    std::memset(&m_bgra[0], 0, sizeof(uint32) * m_bgra.size());
    SetDirty(true);
}

void Framebuffer::Tile::UpdateTexture()
{
    // Update the texture if the buffer contents have changed
    if (GetDirty())
    {
        glBindTexture(GL_TEXTURE_2D, m_tex);
        glTexImage2D(
            GL_TEXTURE_2D,
            0,
            GL_RGBA8,
            GetWidth(),
            GetHeight(),
            0,
            GL_BGRA,
            GL_UNSIGNED_BYTE,
            &m_bgra[0]);
        SetDirty(false);
    }
}

