
#include "scene.h"

#include <limits>

#include "triangle.h"
#include "cornell_box.h"

float Scene::Distance(Vec3f pos)
{
    float dist = std::numeric_limits<float>::max();

    for (uint tri=0; tri<32; tri++)
    {
        const Vec3f v0 = g_cornell_geom[tri * 3 + 0];
        const Vec3f v1 = g_cornell_geom[tri * 3 + 1];
        const Vec3f v2 = g_cornell_geom[tri * 3 + 2];
        dist = std::min(dist, DistancePointTri(pos, v0, v1, v2));
    }

    return dist;
}

bool Scene::Intersect(Vec3f origin, Vec3f dir, float& t)
{
    t = std::numeric_limits<float>::max();
    float u, v;

    for (uint tri=0; tri<32; tri++)
    {
        const Vec3f v0 = g_cornell_geom[tri * 3 + 0];
        const Vec3f v1 = g_cornell_geom[tri * 3 + 1];
        const Vec3f v2 = g_cornell_geom[tri * 3 + 2];

        float cur_t;
        const bool hit = IntersectRayTri(origin, dir, v0, v1, v2, cur_t, u, v);

        if (hit && cur_t < t)
            t = cur_t;
    }

    return t !=std::numeric_limits<float>::max(); 
}

