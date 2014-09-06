
#ifndef RAY_TRI_H
#define RAY_TRI_H

#include "lin_alg.h"

// Fast, Minimum Storage Ray-Triangle Intersection
//
// http://www.jcenligne.fr/download/little3d/
//     jgt%20Fast,%20Minumum%20Storage%20Ray-Triangle%20Intersection.htm

bool IntersectRayTri(const Vec3f& origin,
                     const Vec3f& dir,
                     const Vec3f& v0,
                     const Vec3f& v1,
                     const Vec3f& v2,
                     float& t,
                     float& u,
                     float& v);

#endif // RAY_TRI_H

