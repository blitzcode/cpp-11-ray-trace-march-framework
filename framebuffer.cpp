
#include "framebuffer.h"

Framebuffer::Framebuffer(uint width, uint height)
{
    Resize(width, height);
}

Framebuffer::~Framebuffer()
{
}

void Framebuffer::Resize(uint width, uint height)
{
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

    // Debug fill tiles
    for (Tile &t : m_tiles)
        RenderTile(t);
}

void Framebuffer::RenderTile(Tile& tile)
{
    uint x0, y0, x1, y1;
    tile.GetPosition(x0, y0, x1, y1);

    uint32 *buf = tile.GetBuffer();
    for (uint y=0; y<tile.GetHeight(); y++)
        for (uint x=0; x<tile.GetWidth(); x++)
        {
            buf[x + y * tile.GetWidth()] =
                (((y0 + y) % 2 == 0) ? ((x0 + x) % 2 == 0 ? 0xFF00FF00 : 0xFF000000) :
                                       ((x0 + x) % 2 == 1 ? 0xFF00FF00 : 0xFF000000))
                * (x0 + y0);
        }

    tile.SetDirty(true);
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

}

Framebuffer::Tile::Tile()
{ 
    SetPosition(0, 0, 1, 1);

    glGenTextures(1, &m_tex);
    glBindTexture(GL_TEXTURE_2D, m_tex);

    // No MIP-maps
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // This is important so we don't get artifacts at the borders from
    // texture filtering 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    UpdateTexture();
}

void Framebuffer::Tile::SetPosition(uint x0, uint y0, uint x1, uint y1)
{
    m_x0 = x0;
    m_y0 = y0;
    m_x1 = x1;
    m_y1 = y1;

    // Re-allocate and clear image storage
    m_bgra.resize(GetWidth() * GetHeight());
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

