#pragma once

#include <array>
#include <cmath>
#include <vector>

namespace ndg {

  /// \brief Grid class for a nonuniform grid.
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct Nonuniform_Grid {

    using scalar_type = T;
    using neighbor_index_type = long int;
    using neighbor_index_vector_type = std::array<neighbor_index_type, 2>;

    /// uniform constructor
    Nonuniform_Grid(const scalar_type & left_boundary_pos, const scalar_type & right_boundary_pos,
                    const std::size_t & num_cells)
        : number_of_cells(num_cells) {
      delta_x_vec.resize(number_of_cells, fabs(right_boundary_pos - left_boundary_pos) /
                                              (scalar_type) number_of_cells);
      grid_cell_boundaries.resize(number_of_cells + 1, left_boundary_pos);
      for (std::size_t i = 1; i <= number_of_cells; ++i) {
        grid_cell_boundaries[i] = grid_cell_boundaries[i - 1] + delta_x_vec[i];
      }
    };

    /// uniform constructor
    Nonuniform_Grid(const std::size_t & num_cells) : Nonuniform_Grid(-1.0, 1.0, num_cells){};

    /// nonuniform constructor
    Nonuniform_Grid(const scalar_type & left_boundary_pos,
                    const std::vector<scalar_type> & cell_widths)
        : number_of_cells(cell_widths.size()), delta_x_vec(cell_widths) {

      grid_cell_boundaries.resize(number_of_cells, left_boundary_pos);
      for (std::size_t i = 1; i < number_of_cells; ++i) {
        grid_cell_boundaries[i] = grid_cell_boundaries[i - 1] + delta_x_vec[i];
      }
    };

    /// get the number of cells in the grid
    inline std::size_t get_number_of_cells() const {
      return number_of_cells;
    }

    /// get the volume of the i-th cell
    inline scalar_type dx(const std::size_t & i) const {
      return delta_x_vec[i];
    }

    /// get the surface of the left boundary of the i-th cell
    inline static constexpr scalar_type ds([[maybe_unused]] const std::size_t & i) {
      return 1.0;
    }

    /// get the left boundary coordinate of the i-th cell
    inline scalar_type coord(const std::size_t & i) const {
      return grid_cell_boundaries[i];
    }

    /// return the indices of the neighboring cells
    inline neighbor_index_vector_type get_neighbor_indices(const std::size_t & i) const {
      neighbor_index_vector_type neigh_vec;
      neigh_vec[0] = i - 1;
      neigh_vec[1] = i + 1;
      if (i + 1 >= number_of_cells) {
        neigh_vec[1] = -1;
      }
      return neigh_vec;
    }

    /// transformation from the reference element [-1,1] to the i-th cell.
    inline scalar_type transform_ref_to_global(const scalar_type & x, const std::size_t & i) const {
      const scalar_type y = 0.5 * (coord(i + 1) + coord(i)) + 0.5 * dx(i) * x;
      return y;
    }

    /// Jacobi determinant of the transformation from the reference element [-1,1] to the i-th cell.
    inline scalar_type jacobi_det_transform_ref_to_global(const std::size_t & i) const {
      return 0.5 * dx(i);
    }

  private:

    std::size_t number_of_cells;
    /// width of each cell
    std::vector<scalar_type> delta_x_vec;
    /// vector of all cell boundaries
    std::vector<scalar_type> grid_cell_boundaries;
  };

} // namespace ndg
