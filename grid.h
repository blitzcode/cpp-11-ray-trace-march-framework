
#ifndef GRID_H
#define GRID_H

#include <vector>

#include "types.h"
#include "lin_alg.h"

struct Mesh;

class Grid
{
public:
    Grid(const Mesh& mesh, uint grid_res);

protected:
    uint m_grid_dim[3];
    float m_cell_wdh;
    Vec3f m_aabb_min;
    Vec3f m_aabb_max;
    struct Cell
    {
        std::vector<uint32> m_isect_tri_idx;
    };
    std::vector<Cell> m_cells;

    inline uint GridIdx(uint x, uint y, uint z) const
        { return x + z * m_grid_dim[0] + y * m_grid_dim[0] * m_grid_dim[1]; }

};

#endif // GRID_H

