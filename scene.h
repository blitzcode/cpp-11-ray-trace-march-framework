
#ifndef SCENE_H
#define SCENE_H

#include "lin_alg.h"

class Scene
{
public:
    float Distance(Vec3f pos);
    bool Intersect(Vec3f origin, Vec3f dir, float& t);

protected:
};

#endif // SCENE_H

