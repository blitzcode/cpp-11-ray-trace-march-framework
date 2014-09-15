
#include "renderer.h"

#include <random>
#include <algorithm>
#include <cassert>

#include "sampling.h"
#include "camera.h"
#include "triangle.h"

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

    t = 0.0f;
    for (uint steps=0; steps<max_steps; steps++)
    {
        Vec3f pos = origin + t * dir;
        float dist = DistanceBruteForce(pos);
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

    const Mesh *mesh = m_scene->GetGrid()->GetMesh();

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

                float t, u, v;
                uint32 tri_idx;
                //bool hit = RayMarch(origin, dir, t);
                //bool hit = IntersectBruteForce(origin, dir, t, u, v, tri_idx);
                bool hit = m_scene->GetGrid()->Intersect(origin, dir, t, u, v, tri_idx);

                if (hit)
                {
                    const Mesh::Triangle& tri = mesh->m_triangles[tri_idx];
                    const Vec3f n = Normalize(BarycentricInterpolate(
                        u,
                        v,
                        mesh->m_vertices[tri.v0].n,
                        mesh->m_vertices[tri.v1].n,
                        mesh->m_vertices[tri.v2].n));
                    //Vec3f n = tri.n;
                    col += Vec3f((n + 1.0f) * 0.5f);
                    //col += Vec3f(t / 3);
                }
                else
                    col += Vec3f(float(pixel.y) / float(m_height));
            }

            Vec3f final_col = col / float(m_sample_count);
#define GAMMA_CORRECTION
#ifdef GAMMA_CORRECTION
            const float gamma = 1.0f / 2.0f;
            final_col.x = std::pow(final_col.x, gamma);
            final_col.y = std::pow(final_col.y, gamma);
            final_col.z = std::pow(final_col.z, gamma);
#endif // GAMMA_CORRECTION

            buf[x + y * tile.GetWidth()] = ToBGRA8(final_col);
        }
    }
}

float Renderer::DistanceBruteForce(Vec3f pos)
{
    // Get distance from point by testing all triangles

    const Mesh *mesh = m_scene->GetGrid()->GetMesh();
    float dist = std::numeric_limits<float>::max();

    for (const auto& tri : mesh->m_triangles)
    {
        dist = std::min(dist, DistancePointTri(
            pos,
            mesh->m_vertices[tri.v0].p,
            mesh->m_vertices[tri.v1].p,
            mesh->m_vertices[tri.v2].p));
    }

    return dist;
}

bool Renderer::IntersectBruteForce(
    Vec3f origin,
    Vec3f dir,
    // Output
    float& t,
    float& u,
    float& v,
    uint32& tri_idx)
{
    // Get intersection for ray by testing all triangles

    const Mesh *mesh = m_scene->GetGrid()->GetMesh();
    t = std::numeric_limits<float>::max();

    for (uint32 cur_tri_idx=0; cur_tri_idx<uint32(mesh->m_triangles.size()); cur_tri_idx++)
    {
        // Intersect
        const Mesh::Triangle& tri = mesh->m_triangles[cur_tri_idx];
        float cur_t, cur_u, cur_v;
        const bool hit = IntersectRayTri(
            origin,
            dir,
            mesh->m_vertices[tri.v0].p,
            mesh->m_vertices[tri.v1].p,
            mesh->m_vertices[tri.v2].p,
            cur_t,
            cur_u,
            cur_v);

        // Closer?
        if (hit && cur_t < t)
        {
            t       = cur_t;
            u       = cur_u;
            v       = cur_v;
            tri_idx = cur_tri_idx;
        }
    }

    return t !=std::numeric_limits<float>::max();
}

