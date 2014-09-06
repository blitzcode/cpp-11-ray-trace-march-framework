
#ifndef RENDERER_H
#define RENDERER_H

#include "framebuffer.h"

class Renderer : public Framebuffer
{
public:

protected:
    void RenderTile(Tile& tile) override;

};

#endif // RENDERER_H

