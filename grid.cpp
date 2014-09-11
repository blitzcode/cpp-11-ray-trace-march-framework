
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

    // Grow the bounding box a little. This fixes two problems:
    //
    // 1. Due to the center / half_size conversion in IntersectTriAABB() there can
    //    be small precision problems, missing intersections with triangles directly
    //    on the AABB faces, like when we're trying to build a grid for an actual AABB
    //
    // 2. Geometry lying exactly on the 'positive' faces on the cube would be classified
    //    as having a grid coordinate one beyond the limit
    //
    m_aabb_min -= Vec3f(0.0001f);
    m_aabb_max += Vec3f(0.0001f);

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
        cell_cnt,
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
                        assert(cell_idx < uint(m_cells.size()));
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

    // TODO: Distance field
}

void TraceVec(const char *desc, Vec3f v)
{
    //Trace("%s: %f %f %f", desc, v.x, v.y, v.z);
}
void TraceFloat3(const char *desc, const float *f)
{
    //Trace("%s: %f %f %f", desc, f[0], f[1], f[2]);
}
void TraceInt3(const char *desc, const int *i)
{
    //Trace("%s: %i %i %i", desc, i[0], i[1], i[2]);
}

bool Grid::Intersect(
    Vec3f origin,
    Vec3f dir,
    // Output
    float& t,
    float& u,
    float& v,
    uint32& tri_idx,
    bool debug) const
{
    //debug = false;

    // Traverse the grid and find the closest intersection for the given ray
    //
    // 3DDDA traversal algorithm from "A Fast Voxel Traversal Algorithm For Ray Tracing"
    // http://www.cs.yorku.ca/~amana/research/grid.pdf

    // Intersect with AABB for early rejection and finding the on-grid origin
    // TODO: Does this work if the origin is inside the AABB?
    float t_enter, t_leave;
    if (!IntersectRayAABB(origin, dir, m_aabb_min, m_aabb_max, t_enter, t_leave))
        return false;
    Vec3f origin_adj = origin + dir * t_enter;

    // The ray origin might be slightly outside of the AABB, causing negative
    // grid indices. We also might be right on the max faces, causing a grid
    // one too large
    // TODO: There's a better way of fixing this
    origin_adj = ComponentMin(origin_adj, m_aabb_max-0.0000001f);
    origin_adj = ComponentMax(origin_adj, m_aabb_min);
    t_enter = Dot(dir, origin_adj - origin);

    //if (debug) TraceVec("origin on grid", origin_adj);
    //if (debug) Trace("t_enter: %f", t_enter);

    // In grid units
    const Vec3f origin_grid = (origin_adj - m_aabb_min) / m_cell_wdh;

    if (debug) TraceVec("origin_grid", origin_grid);

    // Signed delta
    const Vec3f delta = ((origin_adj + dir) - m_aabb_min) / m_cell_wdh - origin_grid;

    if (debug) TraceVec("delta", delta);

    // Starting cell. Clamp to valid grid indices as precision problems with
    // IntersectRayAABB() can cause coordinates to be off
    int cell[3] =
    {
        /*
        Clamp(int(std::floor(origin_grid[0])), 0, int(m_grid_dim[0]) - 1),
        Clamp(int(std::floor(origin_grid[1])), 0, int(m_grid_dim[1]) - 1),
        Clamp(int(std::floor(origin_grid[2])), 0, int(m_grid_dim[2]) - 1),
        */
        /*
        int(std::floor(origin_grid[0])),
        int(std::floor(origin_grid[1])),
        int(std::floor(origin_grid[2])),
        */
        int(std::floor((origin_adj.x - m_aabb_min.x) / m_cell_wdh)),
        int(std::floor((origin_adj.y - m_aabb_min.y) / m_cell_wdh)),
        int(std::floor((origin_adj.z - m_aabb_min.z) / m_cell_wdh))
    };

    if (debug) TraceInt3("cell", cell);

    //if (debug)
    {
        if (cell[0] < 0 || cell[0] >= int(m_grid_dim[0]) ||
            cell[1] < 0 || cell[1] >= int(m_grid_dim[1]) ||
            cell[2] < 0 || cell[2] >= int(m_grid_dim[2]))
        {
            assert(false);
            Trace("err");
            //return false;
        }
    }

    // Direction of the ray
    const int step[3] =
    {
        (delta.x < 0.0f) ? -1 : 1,
        (delta.y < 0.0f) ? -1 : 1,
        (delta.z < 0.0f) ? -1 : 1,
    };

    if (debug) TraceInt3("step", step);

    // How much T on the ray equals to a one grid unit change on the axis
    const float delta_t[3] =
    {
        1.0f / std::abs(delta.x),
        1.0f / std::abs(delta.y),
        1.0f / std::abs(delta.z)
    };

    if (debug) TraceFloat3("delta_t", delta_t);

    // Determine the value of T on the ray for each axis where we advance to
    // the next cell on the axis
    float max[3];
    for (uint axis=0; axis<3; axis++)
    {
        if (delta[axis] == 0.0f)
        {
            // No change on this axis during traversal
            max[axis] = std::numeric_limits<float>::max();
            continue;
        }

        // Distance to one cell up or down?
        const float next_cell = (step[axis] == -1) ?
            std::floor(origin_grid[axis]) : std::ceil(origin_grid[axis]);

        // Distance on ray
        max[axis] = std::abs(next_cell - origin_grid[axis]) * delta_t[axis];

        // Full unit when we are on the cell boundary
        if (max[axis] == 0.0f)
            max[axis] = delta_t[axis];
    }

    max[0] += t_enter;
    max[1] += t_enter;
    max[2] += t_enter;

    // TODO
    /*
    // Basic reverse mailboxes, see
    // Ray-Triangle Intersection Algorithm for Modern CPU Architectures
    // http://www.graphicon.ru/2007/proceedings/Papers/Paper_46.pdf
    uint32 mailbox[8];
    for (uint i=0; i<8; i++)
        mailbox[i] = uint32(-1);
    uint mailbox_idx = 0;
    */

    // Step trough the grid until we found an intersection or left the grid
    t = std::numeric_limits<float>::max();
    unsigned int axis_advance;
    while (true)
    {
        // Find on which axis we are stepping next. We need to do this at this
        // early point so we can reject intersection outside of the current cell
        // TODO: PBRT has a branch free version of this
        if (max[0] < max[1])
            axis_advance = (max[0] < max[2]) ? 0 : 2;
        else
            axis_advance = (max[1] < max[2]) ? 1 : 2;

        // Test all triangles in cell
        /*
        if (cell[0] >= 0 && cell[0] < int(m_grid_dim[0]) &&
            cell[1] >= 0 && cell[1] < int(m_grid_dim[1]) &&
            cell[2] >= 0 && cell[2] < int(m_grid_dim[2]))
        */
        {
            // Get triangle index list
            const uint cell_idx = GridIdx(cell[0], cell[1], cell[2]);
            assert(cell_idx < m_cells.size());
            const std::vector<uint32>& tri_indices = m_cells[cell_idx].m_isect_tri_idx;

            for (auto cur_idx : tri_indices)
            {
                const Mesh::Triangle& tri = m_mesh->m_triangles[cur_idx];

                /*
                // Check mailbox
                if (mailbox[0] == cur_idx ||
                    mailbox[1] == cur_idx ||
                    mailbox[2] == cur_idx ||
                    mailbox[3] == cur_idx ||
                    mailbox[4] == cur_idx ||
                    mailbox[5] == cur_idx ||
                    mailbox[6] == cur_idx ||
                    mailbox[7] == cur_idx)
                {
                    continue;
                }
                */

                // Intersect
                float cur_t, cur_u, cur_v;
                const bool hit = IntersectRayTri(
                    origin,
                    dir,
                    m_mesh->m_vertices[tri.v0].p,
                    m_mesh->m_vertices[tri.v1].p,
                    m_mesh->m_vertices[tri.v2].p,
                    cur_t,
                    cur_u,
                    cur_v);

                if (hit)
                {
                    assert(cur_t >= 0.0f);

                    if (cur_t/*-t_enter*/ < max[axis_advance] && // Intersection in current cell?
                        cur_t < t)                   // Closer than any previous?
                    {
                        // Note that this is the T from the adjusted on-grid origin
                        t       = cur_t;
                        u       = cur_u;
                        v       = cur_v;
                        tri_idx = cur_idx;
                    }
                }
                /*
                else
                    mailbox[mailbox_idx++ % 8] = cur_idx;
                */
            }
        }

        // Did we find one or more intersections in this cell ?
        if (t != std::numeric_limits<float>::max())
        {
            // Add distance from origin to on-grid origin
            //t += t_enter;
            return true;
        }

        // Advance to the next cell, abort if outside of grid
        cell[axis_advance] += step[axis_advance];
        if (cell[axis_advance] < 0 || cell[axis_advance] > int(m_grid_dim[axis_advance]) - 1)
            break;
        max[axis_advance] += delta_t[axis_advance];
    }

    return false;
}
