
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "lin_alg.h"

// Fast, Minimum Storage Ray-Triangle Intersection
//
// http://www.jcenligne.fr/download/little3d/
//     jgt%20Fast,%20Minumum%20Storage%20Ray-Triangle%20Intersection.htm
//
inline bool IntersectRayTri(const Vec3f& origin,
                            const Vec3f& dir,
                            const Vec3f& vert0,
                            const Vec3f& vert1,
                            const Vec3f& vert2,
                            float& t,
                            float& u,
                            float& v)
{
    #define EPSILON 0.000001
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

    if (det > -EPSILON && det < EPSILON)
        return false;
    inv_det = 1.0 / det;

    /* calculate distance from vert0 to ray origin */
    SUB(tvec, origin, vert0);

    /* calculate U parameter and test bounds */
    u = DOT(tvec, pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* prepare to test V parameter */
    CROSS(qvec, tvec, edge1);

    /* calculate V parameter and test bounds */
    v = DOT(dir, qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* calculate t, ray intersects triangle */
    t = DOT(edge2, qvec) * inv_det;

    return true;

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

#endif // TRIANGLE_H

