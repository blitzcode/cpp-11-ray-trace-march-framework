
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "lin_alg.h"

// Fast, Minimum Storage Ray-Triangle Intersection
//
// http://www.jcenligne.fr/download/little3d/
//     jgt%20Fast,%20Minumum%20Storage%20Ray-Triangle%20Intersection.htm
//
bool IntersectRayTri(const Vec3f& origin,
                     const Vec3f& dir,
                     const Vec3f& v0,
                     const Vec3f& v1,
                     const Vec3f& v2,
                     float& t,
                     float& u,
                     float& v);

template<typename T> Vector_t<T, 3> TriangleNormal(const Vector_t<T, 3>& v0,
                                                   const Vector_t<T, 3>& v1,
                                                   const Vector_t<T, 3>& v2)
{
    return Normalize(Cross(v1 - v0, v2 - v0));
}

#endif // TRIANGLE_H

