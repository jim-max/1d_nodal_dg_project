#pragma once

#include <1d_nodal_dg_project/config.hpp>

#ifdef OPENMP_FOUND
  #include <omp.h>
#endif

#include <Eigen/Sparse>

#include <1d_nodal_dg_project/utils.hpp>

namespace ndg {

  /// assembles the system matrix and right hand side.
  ///
  /// \todo model_type::advection_term() and
  /// model_type::time_derivative_term() are not yet considered.
  template <typename T>
  class System_Assembler {
  public:

    /// the type of the nodal basis
    using basis_type = T;

    /// the type for scalar variables
    using scalar_type = typename basis_type::scalar_type;
    /// the quadrature type
    using quad_type = typename basis_type::quad_type;

    /// the number of basis elements
    static constexpr std::size_t n_basis_elements = basis_type::n_basis_elements;

    /// the sparse matrix type
    using sparse_matrix_type = Eigen::SparseMatrix<scalar_type>;
    /// a vector type for the right hand side
    using rhs_vector_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

    /// Assemble the sparse system matrix
    ///
    /// \tparam U the grid type
    /// \tparam V the model type
    template <typename U, typename V>
    static sparse_matrix_type assemble_diffusion_matrix(const U & grid, const V & model,
                                                        const scalar_type & eta) {
      using grid_type = U;

      const scalar_type current_t = 0;
      const std::size_t number_of_cells = grid.get_number_of_cells();
      const std::size_t number_of_dofs = number_of_cells * n_basis_elements;

      using triplet_type = Eigen::Triplet<scalar_type>;
      // The std::vector of triplets might contain the elements in arbitrary order,
      // and might even contain duplicated elements that will be summed up by setFromTriplets().
      std::vector<triplet_type> triplet_list;
      triplet_list.reserve(3 * number_of_cells * n_basis_elements * n_basis_elements);

// parallel loop over all cells
#pragma omp parallel for
      for (std::size_t cell_idx = 0; cell_idx < number_of_cells; ++cell_idx) {

        std::vector<triplet_type> local_triplet_list;
        local_triplet_list.reserve(3 * n_basis_elements * n_basis_elements + 2);

        const scalar_type cell_center = grid.coord(cell_idx) + 0.5 * grid.dx(cell_idx);
        const scalar_type diffusion_kappa = model.diffusion_term(cell_center, current_t);

        // get cell neighbors
        const typename grid_type::neighbor_index_vector_type neighbors =
            grid.get_neighbor_indices(cell_idx);

        for (std::size_t i = 0; i < n_basis_elements; ++i) {
          const std::size_t idx_i = i + cell_idx * n_basis_elements;

          // Inner integrals
          for (std::size_t j = 0; j < n_basis_elements; ++j) {
            const std::size_t idx_j = j + cell_idx * n_basis_elements;

            // quadrature in the cell
            scalar_type val_ij = 0.0;
            for (std::size_t q_idx = 0; q_idx < quad_type::get_n_nodes(); ++q_idx) {
              const scalar_type & weight = quad_type::get_weight(q_idx);
              const scalar_type & x_loc = quad_type::get_node(q_idx);
              const scalar_type x_glob = grid.transform_ref_to_global(x_loc, cell_idx);

              const scalar_type & dphi_i = basis_type::get_x_derivative(i, q_idx);
              const scalar_type & dphi_j = basis_type::get_x_derivative(j, q_idx);

              val_ij += weight * diffusion_kappa * dphi_i * dphi_j;

              // TODO exploit fact that we have a nodal basis
              const scalar_type & phi_i = basis_type::get_value(i, q_idx);
              const scalar_type & phi_j = basis_type::get_value(j, q_idx);
              val_ij += weight * model.reaction_term(x_glob, current_t) * phi_i * phi_j;
            }
            // 1/squared det from the derivatives and 1 times from the transformation
            val_ij /= grid.jacobi_det_transform_ref_to_global(cell_idx);

            local_triplet_list.push_back(triplet_type(idx_i, idx_j, val_ij));
          }
        } // end basis element loop

        // edge integrals

        constexpr std::array<scalar_type, 2> normals = {{-1.0, 1.0}};
        constexpr std::array<std::size_t, 2> node_idx_map = {{0, n_basis_elements - 1}};
        constexpr std::array<std::size_t, 2> node_jdx_map = {{n_basis_elements - 1, 0}};

        for (std::size_t k = 0; k < neighbors.size(); ++k) {
          typename grid_type::neighbor_index_type neigh_idx = neighbors[k];
          if (neigh_idx >= 0) {
            // if we are inside the domain
            // (if the neighbor is at the boundary we would have neigh_idx = -1)

            // coordinate of neighbor cell
            const scalar_type x_neigh = grid.coord(neigh_idx) + 0.5 * grid.dx(neigh_idx);
            const scalar_type diffusion_kappa_neigh = model.diffusion_term(x_neigh, current_t);

            const scalar_type gamma = 2.0 * diffusion_kappa_neigh * diffusion_kappa /
                                      (diffusion_kappa_neigh + diffusion_kappa);
            // const scalar_type local_length_scale_hf = fmin(grid.dx(cell_idx), grid.dx(neigh_idx));
            const scalar_type local_length_scale_hf =
                0.5 * (grid.dx(cell_idx) + grid.dx(neigh_idx));

            const scalar_type normal = normals[k];
            const std::size_t node_idx = node_idx_map[k];
            const std::size_t node_jdx = node_jdx_map[k];

            {
              // penalty for jumps:
              const std::size_t i = node_idx;
              const std::size_t idx_i = i + cell_idx * n_basis_elements;
              const std::size_t j = node_jdx;
              const std::size_t idx_j = j + neigh_idx * n_basis_elements;
              scalar_type val_ij = eta * gamma / local_length_scale_hf;
              local_triplet_list.push_back(triplet_type(idx_i, idx_j, -val_ij));
              local_triplet_list.push_back(triplet_type(idx_i, idx_i, val_ij));
            }

            // harmonic mean terms
            for (std::size_t j = 0; j < n_basis_elements; ++j) {
              const std::size_t i = node_idx;
              const std::size_t idx_i = i + cell_idx * n_basis_elements;
              const std::size_t idx_j_1 = j + cell_idx * n_basis_elements;
              const std::size_t idx_j_2 = j + neigh_idx * n_basis_elements;

              const scalar_type dphi_1 = basis_type::get_x_derivative(j, node_idx) /
                                         grid.jacobi_det_transform_ref_to_global(cell_idx);
              const scalar_type dphi_2 = basis_type::get_x_derivative(j, node_jdx) /
                                         grid.jacobi_det_transform_ref_to_global(neigh_idx);

              scalar_type val_1 =
                  harmonic_mean(dphi_1, 0.0, diffusion_kappa, diffusion_kappa_neigh);
              val_1 *= normal;

              scalar_type val_2 =
                  harmonic_mean(0.0, dphi_2, diffusion_kappa, diffusion_kappa_neigh);
              val_2 *= normal;

              local_triplet_list.push_back(triplet_type(idx_i, idx_j_1, -val_1));
              local_triplet_list.push_back(triplet_type(idx_i, idx_j_2, -val_2));
            }

            for (std::size_t i = 0; i < n_basis_elements; ++i) {
              const std::size_t idx_i = i + cell_idx * n_basis_elements;
              const std::size_t idx_j_1 = node_idx + cell_idx * n_basis_elements;
              const std::size_t idx_j_2 = node_jdx + neigh_idx * n_basis_elements;

              const scalar_type & dphi_1 = basis_type::get_x_derivative(i, node_idx) /
                                           grid.jacobi_det_transform_ref_to_global(cell_idx);

              scalar_type val_ij =
                  harmonic_mean(dphi_1, 0.0, diffusion_kappa, diffusion_kappa_neigh);
              val_ij *= normal;
              local_triplet_list.push_back(triplet_type(idx_i, idx_j_1, -val_ij));
              local_triplet_list.push_back(triplet_type(idx_i, idx_j_2, val_ij));
            }

          } else {
            // At boundary
            const scalar_type gamma = diffusion_kappa;
            const scalar_type local_length_scale_hf = grid.dx(cell_idx);
            const scalar_type normal = normals[k];
            const std::size_t node_idx = node_idx_map[k];

            {
              // penalty for jumps:
              const std::size_t i = node_idx;
              const std::size_t idx_i = i + cell_idx * n_basis_elements;
              scalar_type val_ij = eta * gamma / local_length_scale_hf;
              local_triplet_list.push_back(triplet_type(idx_i, idx_i, val_ij));
            }

            // harmonic mean terms
            for (std::size_t j = 0; j < n_basis_elements; ++j) {
              const std::size_t i = node_idx;
              const std::size_t idx_i = i + cell_idx * n_basis_elements;
              const std::size_t idx_j = j + cell_idx * n_basis_elements;

              const scalar_type & dphi = basis_type::get_x_derivative(j, node_idx) /
                                         grid.jacobi_det_transform_ref_to_global(cell_idx);
              scalar_type val = -1.0 * dphi * normal;
              local_triplet_list.push_back(triplet_type(idx_i, idx_j, val));
            }
            for (std::size_t i = 0; i < n_basis_elements; ++i) {
              const std::size_t idx_i = i + cell_idx * n_basis_elements;
              const std::size_t idx_j = node_idx + cell_idx * n_basis_elements;
              const scalar_type & dphi = basis_type::get_x_derivative(i, node_idx) /
                                         grid.jacobi_det_transform_ref_to_global(cell_idx);
              scalar_type val = -1.0 * normal * dphi;
              local_triplet_list.push_back(triplet_type(idx_i, idx_j, val));
            }
          }
        }

#pragma omp critical
        { vec_append(triplet_list, local_triplet_list); }
      }

      // write triplet vector to sparse matrix
      sparse_matrix_type mat(number_of_dofs, number_of_dofs);
      mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
      return mat;
    }

    /// Assemble the right hand side of the system
    ///
    /// \tparam U the grid type
    /// \tparam V the model type
    template <typename U, typename V>
    static rhs_vector_type assemble_rhs(const U & grid, const V & model, const scalar_type & eta) {
      using grid_type = U;

      const scalar_type current_t = 0;

      const std::size_t number_of_cells = grid.get_number_of_cells();
      const std::size_t number_of_dofs = number_of_cells * n_basis_elements;

      rhs_vector_type rhs = rhs_vector_type::Zero(number_of_dofs);

#pragma omp parallel for
      for (std::size_t cell_idx = 0; cell_idx < number_of_cells; ++cell_idx) {

        const scalar_type cell_center = grid.coord(cell_idx) + 0.5 * grid.dx(cell_idx);
        const scalar_type diffusion_kappa = model.diffusion_term(cell_center, current_t);

        // get cell neighbors
        const typename grid_type::neighbor_index_vector_type neighbors =
            grid.get_neighbor_indices(cell_idx);

        for (std::size_t i = 0; i < n_basis_elements; ++i) {
          const std::size_t idx_i = i + cell_idx * n_basis_elements;

          // quadrature in the cell (exploit fact that we have a nodal basis)
          const scalar_type & weight = quad_type::get_weight(i);
          const scalar_type & x_loc = quad_type::get_node(i);
          const scalar_type x_glob = grid.transform_ref_to_global(x_loc, cell_idx);
          scalar_type val = weight * model.rhs_term(x_glob, current_t);
          val *= grid.jacobi_det_transform_ref_to_global(cell_idx);
          rhs[idx_i] += val;

        } // end basis element loop

        // boundary integrals
        constexpr std::array<scalar_type, 2> normals = {{-1.0, 1.0}};
        constexpr std::array<std::size_t, 2> bndry_node_idx_map = {{0, n_basis_elements - 1}};

        for (std::size_t k = 0; k < neighbors.size(); ++k) {
          typename grid_type::neighbor_index_type neigh_idx = neighbors[k];
          if (neigh_idx < 0) {
            // At boundary
            const scalar_type gamma = diffusion_kappa;
            const scalar_type local_length_scale_hf = grid.dx(cell_idx);
            const scalar_type normal = normals[k];
            const std::size_t node_idx = bndry_node_idx_map[k];

            const scalar_type & x_loc = quad_type::get_node(node_idx);
            const scalar_type x_glob = grid.transform_ref_to_global(x_loc, cell_idx);
            const scalar_type boundary_val = model.boundary_value(x_glob, current_t);

            {
              // penalty for jumps:
              const std::size_t i = node_idx;
              const std::size_t idx_i = i + cell_idx * n_basis_elements;
              scalar_type val = eta * gamma / local_length_scale_hf;
              rhs[idx_i] += val * boundary_val;
            }

            for (std::size_t i = 0; i < n_basis_elements; ++i) {
              const std::size_t idx_i = i + cell_idx * n_basis_elements;
              const scalar_type & dphi = basis_type::get_x_derivative(i, node_idx) /
                                         grid.jacobi_det_transform_ref_to_global(cell_idx);
              scalar_type val = normal * dphi * boundary_val;
              rhs[idx_i] -= val;
            }
          }
        } // end neighbor loop

      } // end cell loop
      return rhs;
    }

    template <typename U, typename V>
    static sparse_matrix_type assemble_matrix(const U & grid, const V & model,
                                              const scalar_type & eta,
                                              const scalar_type & delta_t) {
      sparse_matrix_type mat = assemble_diffusion_matrix(grid, model, eta);

      const std::size_t number_of_cells = grid.get_number_of_cells();
      const std::size_t number_of_dofs = number_of_cells * n_basis_elements;
      sparse_matrix_type iden(number_of_dofs, number_of_dofs);
      for (std::size_t i = 0; i < number_of_dofs; ++i) {
        iden.insert(i, i) = 1.0;
      }
      iden.makeCompressed();

      mat = iden + delta_t * mat;
      return mat;
    }

    /// Returns an array of the coordinates of all nodal points
    ///
    /// \tparam U the grid type
    template <typename U>
    static rhs_vector_type get_coordinate_vector(const U & grid) {

      const std::size_t number_of_cells = grid.get_number_of_cells();
      const std::size_t number_of_dofs = number_of_cells * n_basis_elements;

      rhs_vector_type coord_vector = rhs_vector_type::Zero(number_of_dofs);

      // #pragma omp parallel for
      for (std::size_t cell_idx = 0; cell_idx < number_of_cells; ++cell_idx) {
        for (std::size_t i = 0; i < n_basis_elements; ++i) {
          const std::size_t idx_i = i + cell_idx * n_basis_elements;
          const scalar_type & x_loc = quad_type::get_node(i);
          const scalar_type x_glob = grid.transform_ref_to_global(x_loc, cell_idx);
          coord_vector[idx_i] = x_glob;
        } // end basis element loop
      } // end cell loop
      return coord_vector;
    }
  };

} // namespace ndg
