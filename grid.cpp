
#include "grid.h"

#include <algorithm>
#include <cassert>

#include "trace.h"
#include "aabb.h"
#include "triangle.h"
#include "timer.h"

Grid::Grid(std::unique_ptr<Mesh> mesh, uint grid_res)
    : m_mesh(std::move(mesh))
{
    assert(!m_mesh->m_vertices.empty() && !m_mesh->m_triangles.empty());
    assert(grid_res > 0);

    m_mesh->ComputeAABB(m_aabb_min, m_aabb_max);

    // We want to allocate grid_res cells along the longest axis of the mesh
    const Vec3f extends     = m_aabb_max - m_aabb_min;
    const float largest_dim = std::max(std::max(extends.x, extends.y), extends.z);
    m_cell_wdh              = largest_dim / float(grid_res);
    m_grid_dim[0]           = uint(std::ceil(extends.x / m_cell_wdh));
    m_grid_dim[1]           = uint(std::ceil(extends.y / m_cell_wdh));
    m_grid_dim[2]           = uint(std::ceil(extends.z / m_cell_wdh));
    const uint cell_cnt     = m_grid_dim[0] * m_grid_dim[1] * m_grid_dim[2];
    m_cells.resize(cell_cnt);

    Trace(
        "Building %ix%ix%i grid (%i total cells, %.3f cell width)"
        " at (%.3f, %.3f, %.3f) - (%.3f, %.3f, %.3f)"
        " for mesh with %i triangles and %i vertices",
        m_grid_dim[0],
        m_grid_dim[1],
        m_grid_dim[2],
        m_grid_dim[0] * m_grid_dim[1] * m_grid_dim[2],
        m_cell_wdh,
        m_aabb_min.x,
        m_aabb_min.y,
        m_aabb_min.z,
        m_aabb_max.x,
        m_aabb_max.y,
        m_aabb_max.z,
        m_mesh->m_triangles.size(),
        m_mesh->m_vertices.size());

    const double voxel_start_time = TimerGetTick();

    uint intersecting_total = 0;
    uint intersecting_max   = 0;
    for (uint tri_idx=0; tri_idx<m_mesh->m_triangles.size(); tri_idx++)
    {
        const Mesh::Triangle& tri = m_mesh->m_triangles[tri_idx];

        // AABB of current triangle in respect to the grid's origin
        Vec3f tri_aabb_min, tri_aabb_max;
        TriangleAABB(
            m_mesh->m_vertices[tri.v0].p,
            m_mesh->m_vertices[tri.v1].p,
            m_mesh->m_vertices[tri.v2].p,
            tri_aabb_min,
            tri_aabb_max);
        tri_aabb_min -= m_aabb_min;
        tri_aabb_max -= m_aabb_min;

        // Range of potentially intersecting grid cells
        const uint start[3] =
        {
            uint(std::floor(tri_aabb_min.x / m_cell_wdh)),
            uint(std::floor(tri_aabb_min.y / m_cell_wdh)),
            uint(std::floor(tri_aabb_min.z / m_cell_wdh))
        };
        const uint end[3] =
        {
            uint(std::floor(tri_aabb_max.x / m_cell_wdh)),
            uint(std::floor(tri_aabb_max.y / m_cell_wdh)),
            uint(std::floor(tri_aabb_max.z / m_cell_wdh))
        };

        // Determine actual intersections
        uint intersecting_cells = 0;
        for (uint x=start[0]; x<=end[0]; x++)
            for (uint y=start[1]; y<=end[1]; y++)
                for (uint z=start[2]; z<=end[2]; z++)
                {
                    // In world space
                    const Vec3f cell_aabb_min(
                        m_aabb_min.x + x * m_cell_wdh,
                        m_aabb_min.y + y * m_cell_wdh,
                        m_aabb_min.z + z * m_cell_wdh);
                    const Vec3f cell_aabb_max(
                        m_aabb_min.x + (x + 1) * m_cell_wdh,
                        m_aabb_min.y + (y + 1) * m_cell_wdh,
                        m_aabb_min.z + (z + 1) * m_cell_wdh);

                    const bool intersect = IntersectTriAABB(
                        m_mesh->m_vertices[tri.v0].p,
                        m_mesh->m_vertices[tri.v1].p,
                        m_mesh->m_vertices[tri.v2].p,
                        cell_aabb_min,
                        cell_aabb_max);

                    if (intersect)
                    {
                        intersecting_cells++;
                        const uint cell_idx = GridIdx(x, y, z);
                        m_cells[cell_idx].m_isect_tri_idx.push_back(tri_idx);
                    }
                }
        assert(intersecting_cells > 0); // Need to touch at least one cell

        intersecting_total += intersecting_cells;
        intersecting_max    = std::max(intersecting_max, intersecting_cells);
    }

    const double voxel_end_time = TimerGetTick();

    // Statistics
    uint empty_cells = 0, max_tri = 0, total_tri = 0;
    for (const auto& cell : m_cells)
    {
        if (cell.m_isect_tri_idx.empty())
            empty_cells++;
        max_tri = std::max(max_tri, uint(cell.m_isect_tri_idx.size()));
        total_tri += uint(cell.m_isect_tri_idx.size());
    }
    Trace(
        "Voxelization in %.3fs, each triangle intersecting %.2f cells avg. / %i cells max, "
        "%i empty cells (%.2f%%), triangles per cell %.3f avg. / %i max",
        voxel_end_time - voxel_start_time,
        float(intersecting_total) / float(m_mesh->m_triangles.size()),
        intersecting_max,
        empty_cells,
        float(empty_cells) / float(cell_cnt) * 100.0f,
        float(total_tri) / float(cell_cnt),
        max_tri);
}

