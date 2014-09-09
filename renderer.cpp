
#include "renderer.h"

#include <random>
#include <algorithm>
#include <cassert>

#include "sampling.h"
#include "camera.h"

Renderer::Renderer(std::unique_ptr<Scene> scene)
    : m_scene(std::move(scene))
{

}

void Renderer::SetSampleCount(uint cnt)
{
    assert(!WorkerThreadsRunning());
    m_sample_count = std::max(1u, cnt);
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
    //
    // TODO: Use "Enumerating Quasi-Monte Carlo Point Sequences in Elementary Intervals"
    //
    std::vector<Vec2f> smp_loc(m_sample_count);
    std::mt19937 eng;
    std::uniform_real_distribution<float> sample_offs(-0.5f, 0.5f);
    for (uint smp=0; smp<m_sample_count; smp++)
    {
        smp_loc[smp].x = SAMP::HammersleySequence<SAMP::ScrambleNone>(smp, 0, m_sample_count) - 0.5f;
        smp_loc[smp].y = SAMP::HammersleySequence<SAMP::ScrambleNone>(smp, 1, m_sample_count) - 0.5f;
        /*
        smp_loc[smp].x = sample_offs(eng);
        smp_loc[smp].y = sample_offs(eng);
        */
    }

    // Tile parameters
    uint x0, y0, x1, y1;
    tile.GetPosition(x0, y0, x1, y1);
    uint32 *buf = tile.GetBuffer();

    // Camera parameters from the scene
    Matrix44f cam_mat;
    float fov;
    m_scene->GetCameraParameters(fov, cam_mat);

    for (uint y=0; y<tile.GetHeight(); y++)
    {
        if (m_threads_stop)
            return;

        for (uint x=0; x<tile.GetWidth(); x++)
        {
            Vec2ui pixel(x0 + x, y0 + y);

            // Accumulate samples
            //
            // TODO: Use adaptive sampling
            //
            Vec3f col(0.0f);
            for (uint smp=0; smp<m_sample_count; smp++)
            {
                Vec3f origin, dir;
                GenerateRay(cam_mat,
                            pixel,
                            m_width,
                            m_height,
                            smp_loc[smp],
                            false,
                            fov,
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

            buf[x + y * tile.GetWidth()] = ToBGRA8(col / float(m_sample_count));
        }
    }
}

