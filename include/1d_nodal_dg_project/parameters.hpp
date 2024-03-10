

#pragma once

#include <1d_nodal_dg_project/config.hpp>

#include <fstream>
#include <iostream>
#include <optional>

#include <boost/program_options.hpp>

namespace bpo = boost::program_options;

namespace ndg {

  /// Set of parameters for the simulation
  struct Parameters {

    /// the end time
    double t_end = 100.0;
    /// the time step
    double delta_t = 0.1;

    /// left boundary position of the grid
    double grid_xmin = -5.0;
    /// right boundary position of the grid
    double grid_xmax = 5.0;
    /// the number of cells in the grid
    std::size_t number_of_cells = 100;

    /// discontinuous-Galerkin regularization parameter \f$\eta\f$
    double dg_regularization_factor = 25.0;

    /// the initial value u(x) for x < 0
    double initial_value_left = 1.0;
    /// the initial value u(x) for x >= 0
    double initial_value_right = 0.0;
    /// the diffusion cofficient for x < 0
    double diffusion_coefficient_left = 1.0;
    /// the diffusion cofficient for x >= 0
    double diffusion_coefficient_right = 0.25;

    /// the number of threads used
    int number_of_threads = -1;
  };

  /// Reads consolidated CLI options and options read from a config file
  inline std::optional<Parameters> read_parameters(int argc, char ** argv) {
    bpo::options_description config_opts("Options");
    config_opts.add_options()("version", "show version")("help", "produce help message")(
        "config", bpo::value<std::string>(), "read the configuration file");

    Parameters parameters;

    bpo::options_description parameter_opts("Parameters");
    // clang-format off
    parameter_opts.add_options()
    ("t_end", bpo::value<double>(&parameters.t_end)->default_value(parameters.t_end), "")
    ("delta_t", bpo::value<double>(&parameters.delta_t)->default_value(parameters.delta_t), "")
    ("grid_xmin", bpo::value<double>(&parameters.grid_xmin)->default_value(parameters.grid_xmin), "")
    ("grid_xmax", bpo::value<double>(&parameters.grid_xmax)->default_value(parameters.grid_xmax), "")
    ("number_of_cells", bpo::value<std::size_t>(&parameters.number_of_cells)->default_value(parameters.number_of_cells), "")
    ("dg_regularization_factor", bpo::value<double>(&parameters.dg_regularization_factor)->default_value(parameters.dg_regularization_factor), "")
    ("initial_value_left", bpo::value<double>(&parameters.initial_value_left)->default_value(parameters.initial_value_left), "")
    ("initial_value_right", bpo::value<double>(&parameters.initial_value_right)->default_value(parameters.initial_value_right), "")
    ("diffusion_coefficient_left", bpo::value<double>(&parameters.diffusion_coefficient_left)->default_value(parameters.diffusion_coefficient_left), "")
    ("diffusion_coefficient_right", bpo::value<double>(&parameters.diffusion_coefficient_right)->default_value(parameters.diffusion_coefficient_right), "")
    ("number_of_threads", bpo::value<int>(&parameters.number_of_threads)->default_value(parameters.number_of_threads), "")
        // clang-format on
        ;

    bpo::options_description cmdline_opts;
    cmdline_opts.add(config_opts).add(parameter_opts);

    bpo::variables_map var_map;
    store(bpo::command_line_parser(argc, argv).options(cmdline_opts).run(), var_map);
    notify(var_map);

    if (var_map.count("help") > 0) {
      std::cout << cmdline_opts << std::endl;
      return std::nullopt;
    }

    if (var_map.count("config") > 0) {
      const std::string config_file = var_map["config"].as<std::string>();
      std::ifstream in_file_stream(config_file);
      if (!in_file_stream) {
        std::cout << "Error: could not open config file: " << config_file << "\n";
        return std::nullopt;
      }
      store(parse_config_file(in_file_stream, parameter_opts), var_map);
      notify(var_map);
    }

    return parameters;
  }

} // namespace ndg
