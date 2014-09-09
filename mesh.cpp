
#include "mesh.h"

#include <limits>
#include <cassert>
#include <cstdio>

#include "triangle.h"
#include "cornell_box.h"
#include "trace.h"

void Mesh::Clear()
{
    m_triangles.clear();
    m_vertices.clear();
}

void Mesh::CornellBox()
{
    Clear();

    for (uint i=0; i<g_cornell_num_quads; i++)
        AddQuad(&g_cornell_quads[i * 4][0]);
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

void Mesh::AddMesh(const Mesh& mesh)
{
    // Merge the passed mesh into this one

    const uint32 vtx_offs = uint32(m_vertices.size());

    for (const auto& vtx : mesh.m_vertices)
        m_vertices.push_back(vtx);

    for (const auto& tri : mesh.m_triangles)
    {
        m_triangles.push_back(tri);
        m_triangles.back().v0 += vtx_offs;
        m_triangles.back().v1 += vtx_offs;
        m_triangles.back().v2 += vtx_offs;
    }
}

void Mesh::ComputeAABB(Vec3f& aabb_min, Vec3f& aabb_max) const
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

bool Mesh::Read(const char *filename)
{
    // Read simple ASCII mesh format
    //
    // Variant A
    // ---------
    //
    // <tri 0 vtx 0>
    // <tri 0 vtx 1>
    // <tri 0 vtx 2>
    // <tri 1 vtx 0>
    // ...
    //
    // Variant B, Indexed
    // ------------------
    //
    // <num vertices>
    //
    // <vtx 0>
    // <vtx 1>
    // ...
    //
    // <num triangles>
    //
    // <idx 0> <idx 1> <idx 2> ...
    // ...
    //
    // Vertices
    // --------
    //
    // There are several supported vertex specifications
    //
    // x y z
    // x y z nx ny nz
    // x y z nx ny nz u v
    // x y z nx ny nz r g b

    Clear();

    std::FILE *file = std::fopen(filename, "r");
    if (file == nullptr)
    {
        Trace("Mesh::Read() - Can't open file");
        return false;
    }

    int ret;
    char buf[1024];

    #define ERR_EXIT(msg) \
        do { std::fclose(file); Trace("Mesh::Read() - " msg); return false; } while (false)

    // Determine if we got an index mesh by checking if the first line has space separated
    // values (vertex specification) or just a single string (vertex count)
    if (std::fgets(buf, sizeof(buf), file) == NULL)
        ERR_EXIT("Can't read 1st line");
    const bool indexed = std::strchr(buf, ' ') == NULL;
    std::rewind(file);

    // Read vertex count for indexed meshes
    uint num_vtx = 0;
    if (indexed)
    {
        // Vertex count
        ret = std::fscanf(file, "%i\n\n", &num_vtx);
        if (ret != 1)
            ERR_EXIT("Can't get vertex count");
        if (num_vtx < 3)
            ERR_EXIT("Invalid vertex count");
    }

    // Read first vertex to determine format
    {
        const long vtx_start = std::ftell(file);
        if (std::fgets(buf, sizeof(buf), file) == NULL)
            ERR_EXIT("Can't read 1st vertex");
        std::fseek(file, vtx_start, SEEK_SET);
    }

    // Determine vertex specification format
    enum VtxSpec
    {
        VSPos = 0,
        VSPosNormal,
        VSPosNormalUV,
        VSPosNormalRGB
    } vtx_spec;
    float x, y, z, nx, ny, nz, r_or_u, g_or_v, b;
    ret = std::sscanf(buf, "%f %f %f %f %f %f %f %f %f",
        &x, &y, &z, &nx, &ny, &nz, &r_or_u, &g_or_v, &b);
    switch (ret)
    {
        case 3 : vtx_spec = VSPos;          break;
        case 6 : vtx_spec = VSPosNormal;    break;
        case 8 : vtx_spec = VSPosNormalUV;  break;
        case 9 : vtx_spec = VSPosNormalRGB; break;
        default: ERR_EXIT("Invalid vertex spec");
    }

    // Read vertices
    uint vtx_read = 0;
    while (std::feof(file) == 0)
    {
        Vertex vtx;
        vtx.n = Vec3f(0.0f);

        // Position
        ret = std::fscanf(file, "%f %f %f ", &vtx.p.x, &vtx.p.y, &vtx.p.z);
        if (ret != 3)
            ERR_EXIT("Can't read position");

        if (vtx_spec != VSPos)
        {
            // Normal
            ret = std::fscanf(file, "%f %f %f ", &vtx.n.x, &vtx.n.y, &vtx.n.z);
            if (ret != 3)
                ERR_EXIT("Can't read normal");

            // Discard normal or UV
            if (vtx_spec == VSPosNormalUV)
            {
                ret = std::fscanf(file, "%f %f ", &r_or_u, &g_or_v);
                if (ret != 2)
                    ERR_EXIT("Can't read UV");
            }
            else if (vtx_spec == VSPosNormalRGB)
            {
                ret = std::fscanf(file, "%f %f %f ", &r_or_u, &g_or_v, &b);
                if (ret != 3)
                    ERR_EXIT("Can't read RGB");
            }
        }

        m_vertices.push_back(vtx);

        vtx_read++;
        if (indexed && vtx_read >= num_vtx)
            break;
    }

    // Did we read all vertices?
    if (vtx_read == 0 || (indexed && vtx_read != num_vtx))
        ERR_EXIT("Can't read all vertices");
    if (!indexed && (vtx_read % 3 != 0))
        ERR_EXIT("Invalid vertex count");

    // Triangle count
    uint num_tri;
    if (indexed)
    {
        // Index count
        uint num_idx;
        ret = std::fscanf(file, "%i\n\n", &num_idx);
        if (ret != 1)
            ERR_EXIT("Can't get index count");
        if (num_idx < 3 || num_idx % 3 != 0)
            ERR_EXIT("Invalid index count");
        num_tri = num_idx / 3;
    }
    else
        num_tri = vtx_read / 3;
    m_triangles.resize(num_tri);

    // Read / build triangles
    if (indexed)
    {
        for (uint i=0; i<num_tri; i++)
        {
            Triangle& tri = m_triangles[i];

            // Read indices
            ret = std::fscanf(file, "%i %i %i ", &tri.v0, &tri.v1, &tri.v2);
            if (ret != 3)
                ERR_EXIT("Can't read triangle indices");
            if (tri.v0 >= vtx_read || tri.v1 >= vtx_read || tri.v2 >= vtx_read)
                ERR_EXIT("Vertex index out of bounds");

            tri.n = TriangleNormal
                (m_vertices[tri.v0].p, m_vertices[tri.v1].p, m_vertices[tri.v2].p);
        }
    }
    else
    {
        for (uint i=0; i<num_tri; i++)
        {
            Triangle& tri = m_triangles[i];

            // Each set of three vertices forms its own triangle
            tri.v0 = i * 3 + 0;
            tri.v1 = i * 3 + 1;
            tri.v2 = i * 3 + 2;

            tri.n  = TriangleNormal
                (m_vertices[tri.v0].p, m_vertices[tri.v1].p, m_vertices[tri.v2].p);
        }
    }

    // Use face normals if no vertex normals are supplied
    if (vtx_spec == VSPos)
        for (const auto& tri : m_triangles)
            m_vertices[tri.v0].n = m_vertices[tri.v1].n = m_vertices[tri.v2].n = tri.n;

    // Trace some mesh information
    const char spec_to_str[][32] =
    {
        "VSPos",
        "VSPosNormal",
        "VSPosNormalUV",
        "VSPosNormalRGB"
    };
    Trace(
        "Loaded mesh '%s', NumVtx: %i, NumTri: %i, %s, %s",
        filename,
        vtx_read,
        num_tri,
        indexed ? "Indexed" : "Non-Indexed",
        spec_to_str[vtx_spec]);

    std::fclose(file);
    return true;

    #undef ERR_EXIT
}

