
#ifndef FRAMEBUFFER_H
#define FRAMEBUFFER_H

#include <vector>
#include <array>
#include <cstdio>
#include <thread>
#include <mutex>
#include <atomic>

#include <OpenGL/gl.h>

#include "types.h"

class Framebuffer
{
public:
    Framebuffer();
    virtual ~Framebuffer();

    void Resize(uint width, uint height);
    void Draw(uint x, uint y, uint width, uint height);
    void SaveToBMP(const char *filename);
    void RestartRendering();

protected:
    uint m_width                = 1;
    uint m_height               = 1;
    static const uint m_tiles_x = 12;
    static const uint m_tiles_y = 9;
    const uint m_num_cpus;

    class Tile
    {
    public:
        Tile();
        ~Tile() { glDeleteTextures(1, &m_tex); }

        void GetPosition(uint& x0, uint& y0, uint& x1, uint& y1) const
            { x0 = m_x0; y0 = m_y0; x1 = m_x1; y1 = m_y1; }
        void SetPosition(uint x0, uint y0, uint x1, uint y1);

        uint GetWidth()  const { return m_x1 - m_x0; }
        uint GetHeight() const { return m_y1 - m_y0; }
        uint32 * GetBuffer()   { return &m_bgra[0];  }
        void Clear();

        bool GetDirty() const     { return m_dirty;  }
        void SetDirty(bool dirty) { m_dirty = dirty; }
        GLuint GetTexture()       { return m_tex;    }
        void UpdateTexture();

        std::mutex& GetMutex() { return m_mtx; }

    protected:
        std::mutex          m_mtx;
        GLuint              m_tex = -1;
        std::vector<uint32> m_bgra;
        bool                m_dirty; // Does the texture different from the buffer?
        uint                m_x0;
        uint                m_y0;
        uint                m_x1;
        uint                m_y1;
    };

    std::array<Tile, m_tiles_x * m_tiles_y> m_tiles;

    std::vector<uint> m_work_queue; // Indices of tiles left to render
    std::mutex m_work_queue_mtx;    // Mutex to control work queue access
    Tile * GetNextTileFromQueue();

    volatile bool m_threads_stop = false; // Stop signal for the worker threads
    std::vector<std::thread> m_threads;   // Array of worker threads
    std::atomic<uint> m_threads_done;     // Number of threads done
    double m_render_start_time = 0.0;     // Time at which the worker threads started rendering
    void WorkerThread();
    void KillAllWorkerThreads();
    void CreateWorkerThreads();
    void FillWorkQueue();

    virtual void RenderTile(Tile& tile) = 0;
};

#endif // FRAMEBUFFER_H

