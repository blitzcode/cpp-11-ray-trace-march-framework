
#ifndef RENDERER_H
#define RENDERER_H

#include "framebuffer.h"
#include "lin_alg.h"

class Scene;

class Renderer : public Framebuffer
{
public:
    Renderer(Scene *scene);

protected:
    bool RayMarch(Vec3f origin, Vec3f dir, float& t);
    void RenderTile(Tile& tile) override;

    Matrix44f m_camera;
    Scene *m_scene;
};

#endif // RENDERER_H

