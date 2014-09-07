
#include "mesh.h"
#include "triangle.h"
#include "cornell_box.h"

void Mesh::MkCornellBoxMesh(Mesh *mesh)
{
    mesh->m_triangles.clear();

    for (uint i=0; i<g_cornell_geom_num_tri; i++)
    {
        Triangle tri;
        tri.v0 = g_cornell_geom[i * 3 + 0];
        tri.v1 = g_cornell_geom[i * 3 + 1];
        tri.v2 = g_cornell_geom[i * 3 + 2];

        // Geometry only has face normals
        tri.n  = TriangleNormal(tri.v0, tri.v1, tri.v2);
        tri.n0 = tri.n;
        tri.n1 = tri.n;
        tri.n2 = tri.n;

        mesh->m_triangles.push_back(tri);
    }
}

