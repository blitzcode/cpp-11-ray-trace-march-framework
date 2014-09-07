
#include "mesh.h"
#include "triangle.h"
#include "cornell_box.h"

void Mesh::MkCornellBoxMesh(Mesh *mesh)
{
    mesh->m_triangles.clear();
    mesh->m_vertices.clear();

    for (uint i=0; i<g_cornell_geom_num_tri; i++)
    {
        Vec3f v0 = g_cornell_geom[i * 3 + 0];
        Vec3f v1 = g_cornell_geom[i * 3 + 1];
        Vec3f v2 = g_cornell_geom[i * 3 + 2];

        Triangle tri;
        tri.v0 = i * 3 + 0;
        tri.v1 = i * 3 + 1;
        tri.v2 = i * 3 + 2;
        tri.n  = TriangleNormal(v0, v1, v2);
        mesh->m_triangles.push_back(tri);

        // Geometry only has face normals
        mesh->m_vertices.push_back({ v0, tri.n });
        mesh->m_vertices.push_back({ v1, tri.n });
        mesh->m_vertices.push_back({ v2, tri.n });
    }
}

