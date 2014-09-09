
#include "scene.h"

#include "mesh.h"

Scene::Scene(std::unique_ptr<Mesh> mesh, float fov, Matrix44f cam_mat)
    : m_grid(std::move(mesh), 64), m_fov(fov), m_cam_mat(cam_mat)
{

}

