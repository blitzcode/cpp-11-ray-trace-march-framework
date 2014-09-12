
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <algorithm>
#include <limits>

#include "lin_alg.h"

// Fast, Minimum Storage Ray-Triangle Intersection
//
// http://www.jcenligne.fr/download/little3d/
//     jgt%20Fast,%20Minumum%20Storage%20Ray-Triangle%20Intersection.htm
//
inline bool IntersectRayTri(Vec3f origin,
                            Vec3f dir,
                            Vec3f vert0,
                            Vec3f vert1,
                            Vec3f vert2,
                            float& t,
                            float& u,
                            float& v)
{
    // Epsilon smaller than in reference, misses intersections for
    // small (< 1.0^3) scenes / objects otherwise
    #define EPSILON 0.00000001f
    #define CROSS(dest,v1,v2) \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
    #define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
    #define SUB(dest,v1,v2) \
              dest[0]=v1[0]-v2[0]; \
              dest[1]=v1[1]-v2[1]; \
              dest[2]=v1[2]-v2[2];

    float edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
    float det,inv_det;

    /* find vectors for two edges sharing vert0 */
    SUB(edge1, vert1, vert0);
    SUB(edge2, vert2, vert0);

    /* begin calculating determinant - also used to calculate U parameter */
    CROSS(pvec, dir, edge2);

    /* if determinant is near zero, ray lies in plane of triangle */
    det = DOT(edge1, pvec);

#ifdef TEST_CULL           /* define TEST_CULL if culling is desired */
    if (det < EPSILON)
        return false;

    /* calculate distance from vert0 to ray origin */
    SUB(tvec, origin, vert0);

    /* calculate U parameter and test bounds */
    u = DOT(tvec, pvec);
    if (u < 0.0f || u > det)
        return false;

    /* prepare to test V parameter */
    CROSS(qvec, tvec, edge1);

    /* calculate V parameter and test bounds */
    v = DOT(dir, qvec);
    if (v < 0.0f || u + v > det)
        return false;

    /* calculate t, scale parameters, ray intersects triangle */
    t = DOT(edge2, qvec);
    inv_det = 1.0f / det;
    t *= inv_det;
    u *= inv_det;
    v *= inv_det;
#else                    /* the non-culling branch */
    if (det > -EPSILON && det < EPSILON)
        return false;
    inv_det = 1.0f / det;

    /* calculate distance from vert0 to ray origin */
    SUB(tvec, origin, vert0);

    /* calculate U parameter and test bounds */
    u = DOT(tvec, pvec) * inv_det;
    if (u < 0.0f || u > 1.0f)
        return false;

    /* prepare to test V parameter */
    CROSS(qvec, tvec, edge1);

    /* calculate V parameter and test bounds */
    v = DOT(dir, qvec) * inv_det;
    if (v < 0.0f || u + v > 1.0f)
        return false;

    /* calculate t, ray intersects triangle */
    t = DOT(edge2, qvec) * inv_det;
#endif

    return t >= 0.0f;

    #undef EPSILON
    #undef CROSS
    #undef DOT
    #undef SUB
}

template<typename T> Vector_t<T, 3> TriangleNormal(const Vector_t<T, 3>& v0,
                                                   const Vector_t<T, 3>& v1,
                                                   const Vector_t<T, 3>& v2)
{
    return Normalize(Cross(v1 - v0, v2 - v0));
}

template<typename T> void TriangleAABB(const Vector_t<T, 3>& v0,
                                       const Vector_t<T, 3>& v1,
                                       const Vector_t<T, 3>& v2,
                                       Vector_t<T, 3>& aabb_min,
                                       Vector_t<T, 3>& aabb_max)
{
    aabb_min = Vec3f(std::numeric_limits<float>::max());
    aabb_max = Vec3f(std::numeric_limits<float>::min());

    aabb_min = ComponentMin(aabb_min, v0);
    aabb_max = ComponentMax(aabb_max, v0);
    aabb_min = ComponentMin(aabb_min, v1);
    aabb_max = ComponentMax(aabb_max, v1);
    aabb_min = ComponentMin(aabb_min, v2);
    aabb_max = ComponentMax(aabb_max, v2);
}

inline bool ComputeBarycentric(Vec3f pos, Vec3f v0, Vec3f v1, Vec3f v2, float& u, float& v)
{
    // Compute the barycentric coordinates of a point, return if the point is inside
    // the triangle, or more accurate, inside its triangular prism
    //
    // Source: http://www.blackpawn.com/texts/pointinpoly/

    Vec3f e0 = v2 - v0;
    Vec3f e1 = v1 - v0;
    Vec3f e2 = pos - v0;

    const float dot00 = Dot(e0, e0);
    const float dot01 = Dot(e0, e1);
    const float dot02 = Dot(e0, e2);
    const float dot11 = Dot(e1, e1);
    const float dot12 = Dot(e1, e2);

    const float inv_denom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot00 * dot12 - dot01 * dot02) * inv_denom;
    v = (dot11 * dot02 - dot01 * dot12) * inv_denom;

    // Check if point is in triangle
    return (u >= 0) && (v >= 0) && (u + v < 1);
}

template<typename T> T BarycentricInterpolate(float u, float v, T a0, T a1, T a2)
{
    return a1 * u + a2 * v + a0 * (1 - u - v);
}

inline float LineSegMinDistSq(Vec3f a, Vec3f b, Vec3f p)
{
    // Squared distance to the closest point from p on the line segment a b
    Vec3f ab = b - a;
    const float len_sq = Dot(ab, ab);
    float t = Dot(p - a, ab) / len_sq;
    t = Clamp(t, 0.0f, 1.0f);
    const Vec3f proj = a + t * ab;
    return Dot(p - proj, p - proj);
}

inline float DistancePointTri(Vec3f pos, Vec3f v0, Vec3f v1, Vec3f v2)
{
    // Compute the distance between a point and a triangle. This is either the closest
    // point on the plane (if it is inside the triangle), or the closest point on any of
    // the three edges. Note that if we remove the 'inside triangle' case we get a DE for
    // the edges only, allowing us to produce a wireframe rendering
    //
    // TODO: Explore some other, potentially faster methods of computing this
    //       http://www-compsci.swan.ac.uk/~csmark/PDFS/dist.pdf
    //       http://www.ann.jussieu.fr/~frey/papers/divers/
    //           Jones%20M.W.,%203d%20distance%20fields,%20a%20survey.pdf

    float u, v;
    if (ComputeBarycentric(pos, v0, v1, v2, u, v))
    {
        const Vec3f point_on_plane = BarycentricInterpolate(u, v, v0, v1, v2);
        return Distance(pos, point_on_plane);
    }
    else
    {
        return std::sqrt(std::min(LineSegMinDistSq(v0, v1, pos),
                         std::min(LineSegMinDistSq(v0, v2, pos),
                                  LineSegMinDistSq(v1, v2, pos))));
    }
}

#endif // TRIANGLE_H

