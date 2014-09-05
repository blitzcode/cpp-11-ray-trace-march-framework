
#include "font_rendering.h"
#include <cstdlib>
#include <cassert>
#include <OpenGL/glu.h>
#include "misc_fixed_6x12.h"

void FontRendering::DrawStringFixed6x12(
    int x, int y,
    const char *str,
    uint32 color,
    bool vertical)
{
    // Keep a record of all text rendered during a frame
    m_text += str;
    m_text += "\n";

    // Draw text as textured quads
    uint xoffs = 0, yoffs = 0;
    const size_t length = std::strlen(str);
    for (uint i=0; i<length; i++)
    {
        if (str[i] == '\n')
        {
            xoffs = 0;
            yoffs -= 10;
            continue;
        }

        const int idx = str[i];

        // Position of character on the font grid
        const uint tx = (idx % font_grid_wdh);
        const uint ty =
            font_grid_hgt - ((idx - (idx % font_grid_wdh)) / font_grid_wdh + 1);

        // Texture coordinate origin of character
        const float ftx = float(tx * font_char_wdh) / float(font_tex_wdh);
        const float fty = float(ty * font_char_hgt) / float(font_tex_wdh);

        const uint idx_base = m_pos.size();

        // Position & texture coordinates
        //
        // TODO: Re-use last two vertices from the second character on
        m_tex_coords.push_back(Vec2f(
           ftx,
           fty));
        m_pos.push_back(Vec2f(x + xoffs, y + yoffs));
        m_tex_coords.push_back(Vec2f(
           ftx + float(font_char_wdh) / float(font_tex_wdh),
           fty));
        m_pos.push_back(Vec2f(x + xoffs + font_char_wdh, y + yoffs));
        m_tex_coords.push_back(Vec2f(
           ftx + float(font_char_wdh) / float(font_tex_wdh),
           fty + float(font_char_hgt) / float(font_tex_wdh)));
        m_pos.push_back(Vec2f(x + xoffs + font_char_wdh, y + font_char_hgt + yoffs));
        m_tex_coords.push_back(Vec2f(
           ftx,
           fty + float(font_char_hgt) / float(font_tex_wdh)));
        m_pos.push_back(Vec2f(x + xoffs, y + font_char_hgt + yoffs));

        // Colors
        for (uint rep=0; rep<4; rep++)
            m_colors.push_back(color);

        // Indices
        m_indices.push_back(idx_base + 0);
        m_indices.push_back(idx_base + 1);
        m_indices.push_back(idx_base + 3);
        m_indices.push_back(idx_base + 1);
        m_indices.push_back(idx_base + 2);
        m_indices.push_back(idx_base + 3);

        // For drawing vertical text
        if (vertical)
        {
            xoffs = 0;
            yoffs -= 9;
        }
        else
            xoffs += font_char_wdh;
    }
}

void FontRendering::Render(bool filter_font_texture)
{
    assert(m_pos.size() == m_tex_coords.size());
    assert(m_indices.size() % 6 == 0);
    assert(m_indices.size() / 6 * 4 == m_pos.size());

    if (m_indices.empty())
        return;

    // Preserve text
    m_last_frame_text = m_text;
    m_text.clear();

    // Initialize OpenGL font texture
    if (m_font_tex == (GLuint) -1)
    {
        // Create font texture
        glGenTextures(1, &m_font_tex);
        glBindTexture(GL_TEXTURE_2D, m_font_tex);

        // Convert bit packed font image to luminance texture
        uint tex_image[font_tex_wdh * font_tex_wdh];
        for (uint y=0; y<font_img_hgt; y++)
            for (uint x=0; x<font_img_wdh; x++)
            {
                uint *dst = &tex_image[x + y * font_tex_wdh];
                const uint src_idx = x + y * font_img_wdh;
                // Extract bit (reversed in byte), store white / back pixel
                (* dst) = (reinterpret_cast<const uchar *>(g_misc_fixed_6x12_data)[src_idx / 8]
                    & (1 << (7 - (src_idx % 8)))) ? 0xFFFFFFFF : 0x00FFFFFF;
            }

        // Upload texture image
        gluBuild2DMipmapLevels(
            GL_TEXTURE_2D,
            GL_RGBA8,
            font_tex_wdh,
            font_tex_wdh,
            GL_BGRA,
            GL_UNSIGNED_BYTE,
            0,
            0,
            8,
            tex_image);
    }

    // Bind font texture
    glBindTexture(GL_TEXTURE_2D, m_font_tex);
    glEnable(GL_TEXTURE_2D);

    // Texture filtering?
    if (filter_font_texture)
    {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    }
    else
    {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    }

    // Alpha blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(2, GL_FLOAT, 0, &m_pos[0]);
    glTexCoordPointer(2, GL_FLOAT, 0, &m_tex_coords[0]);
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, &m_colors[0]);
    glDrawElements(GL_TRIANGLES, m_indices.size(), GL_UNSIGNED_INT, &m_indices[0]);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);

    // Empty for the next frame
    m_pos.clear();
    m_tex_coords.clear();
    m_colors.clear();
    m_indices.clear();

    glDisable(GL_BLEND);
    glDisable(GL_TEXTURE_2D);
}

