
#ifndef CAMERA_H
#define CAMERA_H

#include "types.h"
#include "lin_alg.h"

inline void GenerateRay(
    const Matrix44f& camera, // Camera transform
    Vec2ui pixel,            // Screen pixel
    uint width,              // Screen width
    uint height,             // ...and height
    Vec2f sample_offs,       // Sample offset [-.5, +.5]
    bool ortho,              // Orthographic or perspective camera?
    float width_or_hfov,     // Width of ortho viewing volume or horizontal FOV degrees
    Vec3f& origin,
    Vec3f& dir)
{
    // Convert fragment coordinates and sample offset to NDC [-1, 1]
    const Vec2f ndc((pixel.x + sample_offs.x) / float(width)  * 2.0f - 1.0f,
                    (pixel.y + sample_offs.y) / float(height) * 2.0f - 1.0f);

    // Generate ray from NDC and camera transform
    const float aspect = float(width) / float(height);
    if (ortho)
    {
        // Orthographic projection. Frame [-w/2, w/2] on X,
        // center interval on Y while keeping aspect
        const float width  = width_or_hfov;
        const float height = float(width) / aspect;
        camera.Transf4x4(Vec3f(ndc.x * (float(width ) / 2.0),
                               ndc.y * (float(height) / 2.0),
                               0.0),
                         origin);
        camera.Transf3x3(Vec3f(0.0f, 0.0f, -1.0f), dir);
    }
    else
    {
        // Perspective projection. Unlike the usual vertical FOV we deal with a horizontal
        // one, just like the orthographic camera defined by its width
        const float hfov   = DegToRad(width_or_hfov);
        const float fov_xs = tan(hfov / 2);
        camera.Transf4x4(Vec3f(0.0), origin);
        camera.Transf3x3
            (Normalize(Vec3f(ndc.x * fov_xs, ndc.y * fov_xs / aspect, -1.0f)), dir);
    }
}

#endif // CAMERA_H

