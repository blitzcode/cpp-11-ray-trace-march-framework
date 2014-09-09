
#include "scene.h"

#include <limits>

#include "triangle.h"
#include "grid.h"

Scene::Scene(std::unique_ptr<Mesh> mesh, float fov, Matrix44f cam_mat)
    : m_mesh(std::move(mesh)), m_fov(fov), m_cam_mat(cam_mat)
{
    Grid grid(* m_mesh.get(), 64);
}

float Scene::Distance(Vec3f pos)
{
    float dist = std::numeric_limits<float>::max();

    for (const auto& tri : m_mesh->m_triangles)
    {
        dist = std::min(dist, DistancePointTri(
                                  pos,
                                  m_mesh->m_vertices[tri.v0].p,
                                  m_mesh->m_vertices[tri.v1].p,
                                  m_mesh->m_vertices[tri.v2].p));
    }

    return dist;
}

bool Scene::Intersect(const Vec3f origin, Vec3f dir, float& t)
{
    t = std::numeric_limits<float>::max();
    float u, v;

    for (const auto& tri : m_mesh->m_triangles)
    {
        float cur_t;
        const bool hit = IntersectRayTri(
                             origin,
                             dir,
                             m_mesh->m_vertices[tri.v0].p,
                             m_mesh->m_vertices[tri.v1].p,
                             m_mesh->m_vertices[tri.v2].p,
                             cur_t,
                             u,
                             v);

        if (hit && cur_t < t)
            t = cur_t;
    }

    return t !=std::numeric_limits<float>::max(); 
}

