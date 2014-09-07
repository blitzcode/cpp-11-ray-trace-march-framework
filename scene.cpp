
#include "scene.h"

#include <limits>

#include "triangle.h"

Scene::Scene(std::unique_ptr<Mesh> mesh)
    : m_mesh(std::move(mesh))
{

}

float Scene::Distance(Vec3f pos)
{
    float dist = std::numeric_limits<float>::max();

    for (const auto& tri : m_mesh->m_triangles)
    {
        dist = std::min(dist, DistancePointTri(pos, tri.v0, tri.v1, tri.v2));
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
        const bool hit = IntersectRayTri(origin, dir, tri.v0, tri.v1, tri.v2, cur_t, u, v);

        if (hit && cur_t < t)
            t = cur_t;
    }

    return t !=std::numeric_limits<float>::max(); 
}

