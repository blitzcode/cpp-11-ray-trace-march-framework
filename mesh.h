
#ifndef MESH_H
#define MESH_H

#include "types.h"
#include "lin_alg.h"

#include <vector>

struct Mesh
{
    struct Triangle
    {
        Vec3f v0;
        Vec3f v1;
        Vec3f v2;

        Vec3f n0;
        Vec3f n1;
        Vec3f n2;

        Vec3f n;
    };

    std::vector<Triangle> m_triangles;

    static void MkCornellBoxMesh(Mesh *mesh);
};

#endif // MESH_H

