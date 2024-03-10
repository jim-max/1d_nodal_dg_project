#pragma once

#include <cstdlib>
#include <tuple>
#include <type_traits>

#include <1d_nodal_dg_project/quadrature/gauss_lobatto_weights.hpp>

namespace ndg {

  /// Wrapper class for the Gauss--Lobatto quadrature
  template <typename T, unsigned M>
  class Gauss_Lobatto_Quadrature {
  public:

    using scalar_type = T;
    /// the number of Gauss--Lobatto nodes
    static constexpr unsigned order = M;

    using available_quadratures_type =
        std::tuple<GL_Quad_1<scalar_type>, GL_Quad_2<scalar_type>, GL_Quad_3<scalar_type>,
                   GL_Quad_4<scalar_type>, GL_Quad_5<scalar_type>, GL_Quad_6<scalar_type>>;
    static_assert((order >= 0) && (order < std::tuple_size_v<available_quadratures_type>),
                  "Invalid quadrature order!");
    using quad_type = typename std::tuple_element_t<order, available_quadratures_type>;

    using array_type = typename quad_type::array_type;

    /// returns the number of nodes
    inline static constexpr unsigned get_n_nodes() {
      return quad_type::n_nodes;
    }

    /// get all node coordinates
    inline static constexpr array_type get_nodes() {
      return quad_type::get_nodes();
    }

    /// get all node weights
    inline static constexpr array_type get_weights() {
      return quad_type::get_weights();
    }

    /// get the i-th node coordinate
    inline static constexpr scalar_type get_node(const std::size_t & idx) {
      return quad_type::get_nodes()[idx];
    }

    /// get the i-th weight
    inline static constexpr scalar_type get_weight(const std::size_t & idx) {
      return quad_type::get_weights()[idx];
    }
  };

} // namespace ndg
