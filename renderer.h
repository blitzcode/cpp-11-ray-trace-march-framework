
#ifndef RENDERER_H
#define RENDERER_H

#include <memory>

#include "framebuffer.h"
#include "lin_alg.h"
#include "scene.h"

class Renderer : public Framebuffer
{
public:
    Renderer(std::unique_ptr<Scene> scene);
    ~Renderer() { KillAllWorkerThreads(); }

protected:
    bool RayMarch(Vec3f origin, Vec3f dir, float& t);
    void RenderTile(Tile& tile) override;

    Matrix44f m_camera;
    std::unique_ptr<Scene> m_scene;
};

#endif // RENDERER_H

