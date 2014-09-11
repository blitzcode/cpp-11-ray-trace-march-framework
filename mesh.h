
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

    void Clear();
    void ComputeAABB(Vec3f& aabb_min, Vec3f& aabb_max) const;
    void Transform(Matrix44f mat);
	void AddQuad(const float *quad_vtx);
	void AddMesh(const Mesh& mesh);
    void NormalizeDimensions();
    bool Read(const char *filename, bool flip_winding = false);
    void CornellBox();
};

#endif // MESH_H

