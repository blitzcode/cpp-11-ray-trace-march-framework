
#ifndef SCENE_H
#define SCENE_H

#include <memory>

#include "types.h"
#include "lin_alg.h"
#include "grid.h"

struct Mesh;

class Scene
{
public:
    Scene(std::unique_ptr<Mesh> mesh, float fov, Matrix44f cam_mat);
    void GetCameraParameters(float& fov, Matrix44f& cam_mat)
        { fov = m_fov; cam_mat = m_cam_mat; }
    inline const Grid * GetGrid() const { return &m_grid; }

protected:
    Grid m_grid;

    // Camera parameters
    float     m_fov;
    Matrix44f m_cam_mat;
};

#endif // SCENE_H

