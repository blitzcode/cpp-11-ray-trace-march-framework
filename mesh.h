
#ifndef MESH_H
#define MESH_H

#include "types.h"
#include "lin_alg.h"

#include <vector>

struct Mesh
{
    struct Triangle
    {
        uint32 v0;
        uint32 v1;
        uint32 v2;
        Vec3f  n;
    };

    struct Vertex
    {
        Vec3f p;
        Vec3f n;
    };

    std::vector<Triangle> m_triangles;
    std::vector<Vertex>   m_vertices;

    static void MkCornellBoxMesh(Mesh *mesh);
};

#endif // MESH_H

