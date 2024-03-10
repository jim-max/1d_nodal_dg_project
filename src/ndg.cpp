
#include <1d_nodal_dg_project/config.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string>
#include <thread>

#include <boost/thread.hpp>
#include <fmt/core.h>
#include <unsupported/Eigen/SparseExtra>

#ifdef OPENMP_FOUND
  #include <omp.h>
#endif

#include <1d_nodal_dg_project/basis/nodal_basis.hpp>
#include <1d_nodal_dg_project/discretization/assemble_system.hpp>
#include <1d_nodal_dg_project/discretization/compute_error.hpp>
#include <1d_nodal_dg_project/grid/nonuniform_grid.hpp>
#include <1d_nodal_dg_project/grid/uniform_grid.hpp>
#include <1d_nodal_dg_project/io/plot_solution.hpp>
#include <1d_nodal_dg_project/models/diffusion_model.hpp>
#include <1d_nodal_dg_project/parameters.hpp>
#include <1d_nodal_dg_project/quadrature/quadrature_gauss_lobatto.hpp>

using namespace ndg;

auto main(int argc, char ** argv) -> int {

  auto opt_parameters = read_parameters(argc, argv);
  if (!opt_parameters.has_value()) {
    return EXIT_FAILURE;
  }
  const Parameters parameters = opt_parameters.value();

#ifdef OPENMP_FOUND
  {
    int max_number_of_threads = 1;
    if (parameters.number_of_threads > 0) {
      max_number_of_threads = parameters.number_of_threads;
    } else {
      // get physical cores (without hyperthreading)
      max_number_of_threads = static_cast<int>(boost::thread::physical_concurrency());
    }
    omp_set_num_threads(max_number_of_threads);
    fmt::println("OpenMP: Number of threads is set to {}", max_number_of_threads);
  }
#endif

  fmt::println("Run diffusion model simulation...");

  using scalar_type = double;

  // the number of basis elements
  constexpr unsigned number_of_nodes = 2;

  // Uniform_Grid or Nonuniform_Grid
  using grid_type = Uniform_Grid<scalar_type>;

  using model_type = Diffusion_Model;
  using quad_type = Gauss_Lobatto_Quadrature<scalar_type, number_of_nodes>;
  using basis_type = Nodal_Basis<quad_type>;
  using assembler_type = System_Assembler<basis_type>;

  using sparse_matrix_type = typename assembler_type::sparse_matrix_type;
  using rhs_vector_type = typename assembler_type::rhs_vector_type;

  // initialize model
  const scalar_type diffusion_coefficient_left = parameters.diffusion_coefficient_left;
  const scalar_type diffusion_coefficient_right = parameters.diffusion_coefficient_right;
  model_type model(diffusion_coefficient_left, diffusion_coefficient_right);

  // set up grid
  grid_type grid(parameters.grid_xmin, parameters.grid_xmax, parameters.number_of_cells);

  // discontinuous-Galerkin regularization parameter
  const scalar_type eta = parameters.dg_regularization_factor;

  // the time step
  const scalar_type delta_t = parameters.delta_t;
  // the end time
  const scalar_type t_end = parameters.t_end;

  // NOTE: in case the system is time-independent we can
  // compute the system matrices beforehand:
  // sparse_matrix_type mat = assembler_type::assemble_matrix(grid, model, eta, delta_t);
  // rhs_vector_type rhs = assembler_type::assemble_rhs(grid, model, eta);
  // Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>> solver;
  // solver.compute(mat);

  // save all coordinates in a vector
  rhs_vector_type coord_vec = assembler_type::get_coordinate_vector(grid);

  // initial condition:
  rhs_vector_type solution = 0.0 * coord_vec;
  for (int i = 0; i < coord_vec.size(); ++i) {
    if (coord_vec[i] < 0.0) {
      solution[i] = parameters.initial_value_left;
    }
    if (coord_vec[i] >= 0.0) {
      solution[i] = parameters.initial_value_right;
    }
  }

  Gnuplot_Plotter plotter;

  // main time loop
  scalar_type current_t = 0.0;
  while (current_t < t_end) {

    // assemble the system
    sparse_matrix_type mat = assembler_type::assemble_matrix(grid, model, eta, delta_t);
    rhs_vector_type rhs = assembler_type::assemble_rhs(grid, model, eta);
    Eigen::SparseLU<Eigen::SparseMatrix<scalar_type>> solver;
    solver.compute(mat);

    current_t += delta_t;
    rhs_vector_type vec_b = solution + delta_t * rhs;
    // solve the linear system
    solution = solver.solve(vec_b);

    plotter.plot_solution(solution, coord_vec, current_t);

    Eigen::VectorXd res = mat * solution - vec_b;

    fmt::print("t = {:.3f} / {:.2f}, ||Ax-b||_2 = {:.2e}     \r", current_t, t_end, res.norm());
    std::this_thread::sleep_for(std::chrono::milliseconds(5));

  } // end main time loop
  fmt::println("");

  // compute the error w.r.t. the exact solution
  const scalar_type l2_error = compute_error<basis_type>(solution, grid, model, 0);
  fmt::println("L2 error = {:e}", l2_error);

  return EXIT_SUCCESS;
}
