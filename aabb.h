
#ifndef AABB_H
#define AABB_H

#include "types.h"
#include "lin_alg.h"
#include "aabb_tri_internal.h"

inline bool IntersectPointAABB(Vec3f p, Vec3f aabb_min, Vec3f aabb_max)
{
    return p.x >= aabb_min.x && p.y >= aabb_min.y && p.z >= aabb_min.z &&
           p.x <= aabb_max.x && p.y <= aabb_max.y && p.z <= aabb_max.z;
}

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

inline bool IntersectRayAABB(
    Vec3f origin,
    Vec3f dir,
    Vec3f aabb_min,
    Vec3f aabb_max,
    float& tmin,
    float& tmax)
{
    // Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
    // "An Efficient and Robust Ray-Box Intersection Algorithm"
    // Journal of graphics tools, 10(1):49-54, 2005
    //
    // http://www.cs.utah.edu/~awilliam/box/

    // TODO: Can be pre-computed
    const Vec3f inv_direction(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);
    const int sign[3] =
    {
        (inv_direction.x < 0.0f) ? 1 : 0,
        (inv_direction.y < 0.0f) ? 1 : 0,
        (inv_direction.z < 0.0f) ? 1 : 0
    };

    const Vec3f aabb[2] = { aabb_min, aabb_max };

    tmin = (aabb[    sign[0]].x - origin.x) * inv_direction.x;
    tmax = (aabb[1 - sign[0]].x - origin.x) * inv_direction.x;

    float tymin = (aabb[    sign[1]].y - origin.y) * inv_direction.y;
    float tymax = (aabb[1 - sign[1]].y - origin.y) * inv_direction.y;

    if ((tmin > tymax) || (tymin > tmax)) 
        return false;
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;

    float tzmin = (aabb[    sign[2]].z - origin.z) * inv_direction.z;
    float tzmax = (aabb[1 - sign[2]].z - origin.z) * inv_direction.z;

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;

    return true;
}

inline bool IntersectRayAABB_Old(
    Vec3f origin,
    Vec3f dir,
    Vec3f aabb_min,
    Vec3f aabb_max,
    float& t)
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

