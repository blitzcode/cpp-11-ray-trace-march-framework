
#ifndef RENDERER_H
#define RENDERER_H

#include "framebuffer.h"

class Renderer : public Framebuffer
{
public:
    Renderer(uint width, uint height) : Framebuffer(width, height) { }

protected:
    void RenderTile(Tile& tile) override;

};

#endif // RENDERER_H

