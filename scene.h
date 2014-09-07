
#ifndef SCENE_H
#define SCENE_H

#include "lin_alg.h"

class Scene
{
public:
    Scene();

    float Distance(Vec3f pos);
    bool Intersect(const Vec3f origin, Vec3f dir, float& t);

protected:
};

#endif // SCENE_H

