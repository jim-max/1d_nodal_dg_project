#pragma once

#include <cmath>
#include <cstddef>

namespace ndg {

  /// \brief Computes the L1 and L2 error of the numerical solution w.r.current_t. the exact solution.
  ///
  /// @tparam B the basis type
  /// @tparam V the vector type
  /// @tparam G the grid type
  /// @tparam M the model type
  template <typename B, typename V, typename G, typename M>
  static double compute_error(const V & solution_vector, const G & grid, const M & model,
                              const double current_t = 0) {
    using basis_type = B;
    using scalar_type = typename basis_type::scalar_type;
    using quad_type = typename basis_type::quad_type;

    static constexpr std::size_t n_basis_elements = basis_type::n_basis_elements;

    const std::size_t number_of_cells = grid.get_number_of_cells();

    scalar_type l1_error = 0.0;
    scalar_type l2_error = 0.0;

#pragma omp parallel for reduction(+ : l1_error, l2_error)
    for (std::size_t cell_idx = 0; cell_idx < number_of_cells; ++cell_idx) {
      for (std::size_t i = 0; i < n_basis_elements; ++i) {
        const std::size_t idx_i = i + cell_idx * n_basis_elements;
        // quadrature in the cell (exploit fact that we have a nodal basis)
        const scalar_type & weight = quad_type::get_weight(i);
        const scalar_type & x_loc = quad_type::get_node(i);
        const scalar_type x_glob = grid.transform_ref_to_global(x_loc, cell_idx);

        const scalar_type exact_val = model.exact_solution(x_glob, current_t);
        const scalar_type err_l1 = std::fabs(exact_val - solution_vector[idx_i]);
        const scalar_type err_l2 = err_l1 * err_l1;
        l1_error += weight * err_l1 * grid.jacobi_det_transform_ref_to_global(cell_idx);
        l2_error += weight * err_l2 * grid.jacobi_det_transform_ref_to_global(cell_idx);
      } // end basis element loop
    } // end cell loop
    l2_error = std::sqrt(l2_error);
    return l2_error;
  }

} // namespace ndg
