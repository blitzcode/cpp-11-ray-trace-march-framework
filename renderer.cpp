
#include "renderer.h"

#include <random>

#include "sampling.h"
#include "camera.h"
#include "scene.h"

Renderer::Renderer(Scene *scene)
    : m_scene(scene)
{
    m_camera.BuildLookAtMatrix(Vec3f(0.0f, 0.0f, -2.0f), Vec3f(0.0f));
}

bool Renderer::RayMarch(Vec3f origin, Vec3f dir, float& t)
{
    const uint max_steps = 128;
    const float min_dist = 0.001f;

    for (uint steps=0; steps<max_steps; steps++)
    {
        Vec3f pos = origin + t * dir;
        float dist = m_scene->Distance(pos);
        t += dist;

        if (dist < min_dist)
            return true;
    }

    return false;
}

void Renderer::RenderTile(Tile& tile)
{
    // Sample locations. For now we just sample with a fixed pattern at each pixel
    const uint num_smp = 16;
    std::array<Vec2f, num_smp> smp_loc;
    std::mt19937 eng;
    std::uniform_real_distribution<float> sample_offs(-0.5f, 0.5f);
    for (uint smp=0; smp<num_smp; smp++)
    {
        smp_loc[smp].x = SAMP::HammersleySequence<SAMP::ScrambleNone>(smp, 0, num_smp) - 0.5f;
        smp_loc[smp].y = SAMP::HammersleySequence<SAMP::ScrambleNone>(smp, 1, num_smp) - 0.5f;
        /*
        smp_loc[smp].x = sample_offs(eng);
        smp_loc[smp].y = sample_offs(eng);
        */
    }

    uint x0, y0, x1, y1;
    tile.GetPosition(x0, y0, x1, y1);

    uint32 *buf = tile.GetBuffer();

    for (uint y=0; y<tile.GetHeight(); y++)
    {
        if (m_threads_stop)
            return;

        for (uint x=0; x<tile.GetWidth(); x++)
        {
            Vec2ui pixel(x0 + x, y0 + y);

            // Accumulate samples
            Vec3f col(0.0f);
            for (uint smp=0; smp<num_smp; smp++)
            {
                Vec3f origin, dir;
                GenerateRay(m_camera,
                            pixel,
                            m_width,
                            m_height,
                            smp_loc[smp],
                            false,
                            60.0f,
                            origin,
                            dir);

                float t = 0.0f;
                //bool hit = RayMarch(origin, dir, t);
                bool hit = m_scene->Intersect(origin, dir, t);

                if (hit)
                    col += Vec3f(t / 3.0f);
                else
                    col += Vec3f(float(pixel.y) / float(m_height), 0.0f, 0.0f);
            }

            buf[x + y * tile.GetWidth()] = ToBGRA8(col / float(num_smp));
        }
    }
}

