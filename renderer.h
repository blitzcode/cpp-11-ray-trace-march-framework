
#ifndef RENDERER_H
#define RENDERER_H

#include "framebuffer.h"
#include "lin_alg.h"

class Renderer : public Framebuffer
{
public:
    Renderer();

protected:
    void RenderTile(Tile& tile) override;

    Matrix44f m_camera;
};

#endif // RENDERER_H

