#define CATCH_CONFIG_MAIN

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <type_traits>

#include <Eigen/Dense>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <unsupported/Eigen/SparseExtra>

#include <1d_nodal_dg_project/basis/nodal_basis.hpp>
#include <1d_nodal_dg_project/discretization/assemble_system.hpp>
#include <1d_nodal_dg_project/discretization/compute_error.hpp>
#include <1d_nodal_dg_project/grid/uniform_grid.hpp>
#include <1d_nodal_dg_project/quadrature/quadrature_gauss_lobatto.hpp>
#include <test_models/discontinuous_diffusion_test_model.hpp>
#include <test_models/linear_test_model.hpp>
#include <test_models/quadratic_test_model.hpp>

using namespace ndg;

TEMPLATE_TEST_CASE("approximation error within range for model", "[integration][model][template]",
                   Linear_Test_Model, Quadratic_Test_Model, Discontinuous_Diffusion_Test_Model) {

  using model_type = TestType;

  using scalar_type = double;
  constexpr unsigned pol_order = 2;

  using grid_type = Uniform_Grid<scalar_type>;
  using quad_type = Gauss_Lobatto_Quadrature<scalar_type, pol_order>;
  using basis_type = Nodal_Basis<quad_type>;
  using assembler_type = System_Assembler<basis_type>;
  using sparse_matrix_type = typename assembler_type::sparse_matrix_type;
  using rhs_vector_type = typename assembler_type::rhs_vector_type;

  model_type model;

  constexpr int n_cells = 200;

  // discontinuous Galerkin regularization parameter:
  constexpr scalar_type eta = 5.0;

  grid_type grid(-1, 1, n_cells);

  sparse_matrix_type mat = assembler_type::assemble_diffusion_matrix(grid, model, eta);
  rhs_vector_type rhs = assembler_type::assemble_rhs(grid, model, eta);
  rhs_vector_type coord_vec = assembler_type::get_coordinate_vector(grid);

  Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>> solver;
  solver.compute(mat);
  Eigen::VectorXd solution = solver.solve(rhs);

  Eigen::VectorXd res = mat * solution - rhs;
  INFO("||Ax-b||_2 = " << res.norm());

  const scalar_type l2_error = compute_error<basis_type>(solution, grid, model, 0);
  INFO("L2-Error = " << l2_error);

  constexpr scalar_type residual_tolerance = 1.0e-10;
  REQUIRE(res.norm() < residual_tolerance);

  // higher tolerance for the discontinuous test model
  constexpr scalar_type l2_error_tolerance =
      (std::is_same_v<model_type, Discontinuous_Diffusion_Test_Model>) ? 1.0e-1 : 1.0e-11;
  REQUIRE(l2_error < l2_error_tolerance);
}

TEMPLATE_TEST_CASE_SIG("approximation error within range for polynomial order",
                       "[integration][pol_order][template]",
                       ((unsigned test_pol_order), test_pol_order), 0, 1, 2, 3, 4) {
  using scalar_type = double;

  constexpr unsigned pol_order = test_pol_order;
  // tolerance values for each pol_order
  constexpr std::array<scalar_type, 5> l2_error_tolerances{0.5, 1.0e-5, 1.0e-11, 1.0e-9, 1.0e-9};
  constexpr scalar_type l2_error_tolerance = l2_error_tolerances[test_pol_order];

  using model_type = Quadratic_Test_Model;

  using grid_type = Uniform_Grid<scalar_type>;
  using quad_type = Gauss_Lobatto_Quadrature<scalar_type, pol_order>;
  using basis_type = Nodal_Basis<quad_type>;
  using assembler_type = System_Assembler<basis_type>;
  using sparse_matrix_type = typename assembler_type::sparse_matrix_type;
  using rhs_vector_type = typename assembler_type::rhs_vector_type;

  model_type model;

  constexpr int n_cells = 200;

  // discontinuous Galerkin regularization parameter:
  constexpr scalar_type eta = 5.0;

  grid_type grid(-1, 1, n_cells);

  sparse_matrix_type mat = assembler_type::assemble_diffusion_matrix(grid, model, eta);
  rhs_vector_type rhs = assembler_type::assemble_rhs(grid, model, eta);
  rhs_vector_type coord_vec = assembler_type::get_coordinate_vector(grid);

  Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>> solver;
  solver.compute(mat);
  Eigen::VectorXd solution = solver.solve(rhs);

  Eigen::VectorXd res = mat * solution - rhs;
  INFO("||Ax-b||_2 = " << res.norm());

  const scalar_type l2_error = compute_error<basis_type>(solution, grid, model, 0);
  INFO("L2-Error = " << l2_error);

  constexpr scalar_type residual_tolerance = 1.0e-10;
  REQUIRE(res.norm() < residual_tolerance);
  REQUIRE(l2_error < l2_error_tolerance);
}
