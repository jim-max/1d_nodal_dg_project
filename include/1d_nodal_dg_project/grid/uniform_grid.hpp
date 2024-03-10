#pragma once

#include <array>
#include <cmath>
#include <stdlib.h>

namespace ndg {

  /// Grid class for a uniform grid.
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct Uniform_Grid {

    using scalar_type = T;
    using neighbor_index_type = long int;
    using neighbor_index_vector_type = std::array<neighbor_index_type, 2>;

    Uniform_Grid(const scalar_type & a, const scalar_type & b, const std::size_t & n)
        : number_of_cells(n), left_boundary_pos(a), right_boundary_pos(b),
          delta_x(fabs(left_boundary_pos - right_boundary_pos) / (scalar_type) number_of_cells){};

    Uniform_Grid(const std::size_t & n) : Uniform_Grid(-1.0, 1.0, n){};

    /// get the number of cells in the grid
    inline std::size_t get_number_of_cells() const {
      return number_of_cells;
    }

    /// g the volume of the idx-th cell
    inline scalar_type dx([[maybe_unused]] const std::size_t & idx) const {
      return delta_x;
    }

    /// get the surface of the left boundary of the idx-th cell
    inline static constexpr scalar_type ds([[maybe_unused]] const std::size_t & idx) {
      return 1.0;
    }

    /// get the left boundary coordinate of the idx-th cell
    inline scalar_type coord(const std::size_t & idx) const {
      return ((scalar_type) idx * delta_x + left_boundary_pos);
    }

    /// return the indices of the neighboring cells
    inline neighbor_index_vector_type get_neighbor_indices(const std::size_t & idx) const {
      neighbor_index_vector_type neigh_vec;
      neigh_vec[0] = idx - 1;
      neigh_vec[1] = idx + 1;
      if (idx + 1 >= number_of_cells) {
        neigh_vec[1] = -1;
      }
      return neigh_vec;
    }

    /// transformation from the reference element [-1,1] to the idx-th cell.
    inline scalar_type transform_ref_to_global(const scalar_type & x,
                                               const std::size_t & idx) const {
      const scalar_type y = 0.5 * (coord(idx + 1) + coord(idx)) + 0.5 * delta_x * x;
      return y;
    }

    /// Jacobi determinant of the transformation from the reference element [-1,1] to the idx-th cell.
    inline scalar_type
    jacobi_det_transform_ref_to_global([[maybe_unused]] const std::size_t & idx) const {
      return 0.5 * delta_x;
    }

  private:

    const std::size_t number_of_cells;
    const scalar_type left_boundary_pos;
    const scalar_type right_boundary_pos;
    /// the grid cell size
    const scalar_type delta_x;
  };

} // namespace ndg
