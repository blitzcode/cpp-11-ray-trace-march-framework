
#ifndef AABB_H
#define AABB_H

#include "types.h"
#include "lin_alg.h"
#include "aabb_tri_internal.h"

inline bool IntersectTriAABB(Vec3f v0, Vec3f v1, Vec3f v2, Vec3f aabb_min, Vec3f aabb_max)
{
    // Wrapper, see aabb_tri_internal.h

    double center[3] = { (aabb_min.x + aabb_max.x) * 0.5f,
                         (aabb_min.y + aabb_max.y) * 0.5f,
                         (aabb_min.z + aabb_max.z) * 0.5f };

    double half_size[3] = { (aabb_max.x - aabb_min.x) * 0.5f,
                            (aabb_max.y - aabb_min.y) * 0.5f,
                            (aabb_max.z - aabb_min.z) * 0.5f };

    double tri[3][3] = { { v0.x, v0.y, v0.z },
                         { v1.x, v1.y, v1.z },
                         { v2.x, v2.y, v2.z } };

    return triBoxOverlap(center, half_size, tri) == 1;
}

inline bool IntersectRayAABB(Vec3f origin, Vec3f dir, Vec3f aabb_min, Vec3f aabb_max, float& t)
{
    uint i, quadrant[3], which_plane;
    bool inside = true;
    float max_t[3], candidate_plane[3];

    for (i=0; i<3; i++)
    {
        if (origin.m_vec[i] < aabb_min.m_vec[i])
        {
            quadrant[i] = 1;
            candidate_plane[i] = aabb_min.m_vec[i];
            inside = false;
        }
        else if (origin.m_vec[i] > aabb_max.m_vec[i])
        {
            quadrant[i] = 1;
            candidate_plane[i] = aabb_max.m_vec[i];
            inside = false;
        }
        else
            quadrant[i] = 0;
    }

    if (inside)
    {
        t = 0.0f;
        return true;
    }

    for (i=0; i<3; i++)
    {
        if (quadrant[i] != 0 && std::abs(dir.m_vec[i]) > 0.00001f)
            max_t[i] = (candidate_plane[i] - origin.m_vec[i]) / dir.m_vec[i];
        else
            max_t[i] = -1.0f;
    }

    which_plane = 0;

    if (max_t[which_plane] < max_t[1])
        which_plane = 1;
    if (max_t[which_plane] < max_t[2])
        which_plane = 2;

    if (max_t[which_plane] < 0.0f)
        return false;

    for (i=0; i<3; i++)
    {
        if (which_plane != i)
        {
            origin.m_vec[i] = origin.m_vec[i] + max_t[which_plane] * dir.m_vec[i];

            if (origin.m_vec[i] < aabb_min.m_vec[i] - 0.00001f ||
                origin.m_vec[i] > aabb_max.m_vec[i] + 0.00001f)
            {
                return false;
            }
        }
        else 
            origin.m_vec[i] = candidate_plane[i];
    }

    t = max_t[which_plane];
    return true;
}

#endif // AABB_H

