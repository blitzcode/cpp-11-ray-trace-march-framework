
#ifndef SCENE_H
#define SCENE_H

#include <memory>

#include "lin_alg.h"
#include "mesh.h"

class Scene
{
public:
    Scene(std::unique_ptr<Mesh> mesh, float fov, Matrix44f cam_mat);

    float Distance(Vec3f pos);
    bool Intersect(const Vec3f origin, Vec3f dir, float& t);

    void GetCameraParameters(float& fov, Matrix44f& cam_mat)
        { fov = m_fov; cam_mat = m_cam_mat; }

protected:
    std::unique_ptr<Mesh> m_mesh;

    // Camera parameters
    float     m_fov;
    Matrix44f m_cam_mat;
};

#endif // SCENE_H

