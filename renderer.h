
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

    void SetSampleCount(uint cnt);
    uint GetSampleCount() const { return m_sample_count; }

protected:
    bool RayMarch(Vec3f origin, Vec3f dir, float& t);
    void RenderTile(Tile& tile) override;

    std::unique_ptr<Scene> m_scene;
    uint m_sample_count = 16;
};

#endif // RENDERER_H

