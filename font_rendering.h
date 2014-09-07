
#ifndef FONT_RENDERING_H
#define FONT_RENDERING_H

#include <vector>
#include <string>

#include <OpenGL/gl.h>

#include "types.h"
#include "lin_alg.h"

class FontRendering
{
public:
    FontRendering() : m_font_tex((GLuint) -1) { }

    ~FontRendering()
    {
        // Delete OpenGL font texture
        if (m_font_tex != (GLuint) - 1)
            glDeleteTextures(1, &m_font_tex);
    }

    void DrawStringFixed6x12(
        int x, int y,
        const char *str,
        uint32 color = 0xFFFFFFFF,
        bool vertical = false);

    void Render(bool filter_font_texture = false);

protected:
    GLuint              m_font_tex;
    std::vector<Vec2f>  m_pos;
    std::vector<Vec2f>  m_tex_coords;
    std::vector<uint32> m_colors;
    std::vector<GLuint> m_indices;

    std::string m_text;
    std::string m_last_frame_text;

    const uint font_grid_wdh = 16 ;
    const uint font_grid_hgt = 16 ;
    const uint font_img_wdh  = 96 ;
    const uint font_img_hgt  = 192;
    const uint font_char_wdh = 6  ;
    const uint font_char_hgt = 12 ;
    const uint font_tex_wdh  = 256;

};

#endif // FONT_RENDERING_H

