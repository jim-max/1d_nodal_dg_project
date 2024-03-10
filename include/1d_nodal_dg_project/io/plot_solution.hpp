#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <fmt/core.h>

namespace ndg {

  /// Class to support gnuplot plotting
  class Gnuplot_Plotter {
  public:

    Gnuplot_Plotter() {
      gnuplot_pipe = popen("gnuplot -persist", "w");
      if (gnuplot_pipe == nullptr) {
        fmt::println("Warning: could not open a gnuplot pipe!");
      }
      fmt::println(gnuplot_pipe, "set terminal x11");
      fmt::println(gnuplot_pipe, "set title \"1d nodal dg project\"");
    };

    ~Gnuplot_Plotter() {
      pclose(gnuplot_pipe);
    }

    /// \brief update the plot
    ///
    /// \param u vector of values
    /// \param x vector of value positions
    /// \param t current simulation time
    /// \tparam T vector type supporting []-access and .size()
    /// \tparam U vector type supporting []-access and .size()
    template <typename T, typename U>
    void plot_solution(const T & u, const U & x, const double & t = 0) {
      if (gnuplot_pipe == nullptr) {
        return;
      }

      const std::size_t n_dofs = x.size();
      const double xmin = x[0];
      const double xmax = x[n_dofs - 1];

      fmt::println(gnuplot_pipe, "set title \"t = {:f}\"", t);
      fmt::println(gnuplot_pipe, "set xrange [{:f}:{:f}]", xmin, xmax);

      fmt::println(gnuplot_pipe,
                   "plot '-' using 1:2 with lines lt rgb \"#8da0cb\" lw 2 title \"u(x)\"");

      for (std::size_t i = 0; i < n_dofs; ++i) {
        fmt::println(gnuplot_pipe, "{:e}\t{:e}", x[i], u[i]);
      }
      fmt::println(gnuplot_pipe, "e");

      fflush(gnuplot_pipe); // flush the gnuplot_pipe to update the plot
    }

  private:

    std::FILE * gnuplot_pipe = nullptr;
  };

} // namespace ndg
