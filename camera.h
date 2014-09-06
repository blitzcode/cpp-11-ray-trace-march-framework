
#ifndef CAMERA_H
#define CAMERA_H

#include "types.h"
#include "lin_alg.h"

void GenerateRay(const Matrix44f& camera, // Camera transform
                 Vec2ui pixel,            // Screen pixel
                 uint width,              // Screen width
                 uint height,             // ...and height
                 Vec2f sample_offs,       // Sample offset [-.5, +.5]
                 bool ortho,              // Orthographic or perspective camera?
                 float width_or_hfov,     // Width of ortho viewing volume or horizontal FOV degrees
                 Vec3f& origin,
                 Vec3f& dir);

#endif // CAMERA_H

