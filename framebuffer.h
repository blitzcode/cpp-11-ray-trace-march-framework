
#ifndef FRAMEBUFFER_H
#define FRAMEBUFFER_H

#include <vector>
#include <array>
#include <cstdio>
#include <thread>
#include <mutex>

#include <OpenGL/gl.h>

#include "types.h"

class Framebuffer
{
public:
    Framebuffer(uint width, uint height);
    ~Framebuffer();

    void Resize(uint width, uint height);
    void Draw(uint x, uint y, uint width, uint height);
    void SaveToBMP(const char *filename);

protected:
    uint m_width;
    uint m_height;
    static const uint m_tiles_x = 8;
    static const uint m_tiles_y = 8;

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
        uint32 * GetBuffer() { return &m_bgra[0]; }

        bool GetDirty() const     { return m_dirty;  }
        void SetDirty(bool dirty) { m_dirty = dirty; }

        std::mutex& GetMutex() { return m_mtx; }

        GLuint GetTexture() { return m_tex; }
        void UpdateTexture();

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

    void RenderTile(Tile& tile);

    std::array<Tile, m_tiles_x * m_tiles_y> m_tiles;
};

#endif // FRAMEBUFFER_H

