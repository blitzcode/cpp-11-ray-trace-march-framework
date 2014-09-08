
#include "mesh.h"

#include <limits>
#include <cassert>
#include <cstdio>

#include "triangle.h"
#include "cornell_box.h"

void Mesh::MkCornellBoxMesh(Mesh *mesh)
{
    mesh->m_triangles.clear();
    mesh->m_vertices.clear();

    for (uint i=0; i<g_cornell_num_quads; i++)
        mesh->AddQuad(&g_cornell_quads[i * 4][0]);
}

void Mesh::AddQuad(
    const float *quad_vtx) // 4x3 floats making up the quad
{
    const Vec3f n = TriangleNormal(
        Vec3f(&quad_vtx[0]),
        Vec3f(&quad_vtx[3]),
        Vec3f(&quad_vtx[6]));

    const uint32 base_idx = uint32(m_vertices.size());
    m_vertices.push_back({ Vec3f(&quad_vtx[0]), n });
    m_vertices.push_back({ Vec3f(&quad_vtx[3]), n });
    m_vertices.push_back({ Vec3f(&quad_vtx[6]), n });
    m_vertices.push_back({ Vec3f(&quad_vtx[9]), n });

    Triangle tri;
    tri.n = n;

    tri.v0 = base_idx + 0;
    tri.v1 = base_idx + 1;
    tri.v2 = base_idx + 2;
    m_triangles.push_back(tri);

    tri.v0 = base_idx + 0;
    tri.v1 = base_idx + 2;
    tri.v2 = base_idx + 3;
    m_triangles.push_back(tri);
}

void Mesh::ComputeAABB(Vec3f& aabb_min, Vec3f& aabb_max)
{
    // Handle empty mesh
    if (m_triangles.empty())
    {
        aabb_min = aabb_max = Vec3f(0.0f);
        return;
    }

    aabb_min = Vec3f(std::numeric_limits<float>::max());
    aabb_max = Vec3f(std::numeric_limits<float>::min());

    // Only look at vertices referenced by actual triangles
    for (const auto& tri : m_triangles)
    {
        aabb_min = ComponentMin(aabb_min, m_vertices[tri.v0].p);
        aabb_max = ComponentMax(aabb_max, m_vertices[tri.v0].p);
        aabb_min = ComponentMin(aabb_min, m_vertices[tri.v1].p);
        aabb_max = ComponentMax(aabb_max, m_vertices[tri.v1].p);
        aabb_min = ComponentMin(aabb_min, m_vertices[tri.v2].p);
        aabb_max = ComponentMax(aabb_max, m_vertices[tri.v2].p);
    }
}

void Mesh::Transform(Matrix44f mat)
{
    // Transform normals with the inverse transpose
    Matrix44f inverse_transpose = mat;
    bool ret = inverse_transpose.Invert();
    assert(ret);
    inverse_transpose.Transpose4x4();

    // Face normals
    for (auto& tri : m_triangles)
    {
        inverse_transpose.Transf3x3(tri.n);
        tri.n = Normalize(tri.n);
    }

    // Vertex positions and normals
    for (auto& vtx : m_vertices)
    {
        mat.Transf4x4(vtx.p);
        inverse_transpose.Transf3x3(vtx.n);
        vtx.n = Normalize(vtx.n);
    }
}

void Mesh::NormalizeDimensions()
{
    // Center and scale mesh to [-.5, +.5] / 1 x 1 x 1 centered at (0, 0, 0)

    Vec3f aabb_min, aabb_max;
    ComputeAABB(aabb_min, aabb_max);
    const Vec3f center = (aabb_min + aabb_max) / 2.0f;
    const Vec3f extends = aabb_max - aabb_min;

    Matrix44f mat_trans;
    mat_trans.Translation(-center.x, -center.y, -center.z);

    Matrix44f mat_scale;
    mat_scale.Scaling(1.0f / (std::max(std::max(extends.x, extends.y), extends.z)));

    Transform(mat_trans * mat_scale);
}

