
#ifndef GRID_H
#define GRID_H

#include <vector>
#include <memory>

#include "types.h"
#include "lin_alg.h"
#include "mesh.h"

class Grid
{
public:
    Grid(std::unique_ptr<Mesh> mesh, uint grid_res);
    inline const Mesh * GetMesh() const { return m_mesh.get(); }
    bool Intersect(
        Vec3f origin,
        Vec3f dir,
        float& t,
        float& u,
        float& v,
        uint32& tri_idx) const;

protected:
    std::unique_ptr<Mesh> m_mesh;

    uint m_grid_dim[3];
    float m_cell_wdh;
    Vec3f m_aabb_min;
    Vec3f m_aabb_max;
    struct Cell
    {
        // TODO: Probably want to store all this tightly packed in a single uint32 array
        //       to reduce memory consumption and cache misses
        std::vector<uint32> m_isect_tri_idx;
    };
    std::vector<Cell> m_cells;

    inline uint GridIdx(uint x, uint y, uint z) const
        { return x + z * m_grid_dim[0] + y * m_grid_dim[0] * m_grid_dim[2]; }

};

#endif // GRID_H

