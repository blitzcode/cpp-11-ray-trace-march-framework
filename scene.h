
#ifndef SCENE_H
#define SCENE_H

#include <memory>

#include "lin_alg.h"
#include "mesh.h"

class Scene
{
public:
    Scene(std::unique_ptr<Mesh> mesh);

    float Distance(Vec3f pos);
    bool Intersect(const Vec3f origin, Vec3f dir, float& t);

protected:
    std::unique_ptr<Mesh> m_mesh;
};

#endif // SCENE_H

