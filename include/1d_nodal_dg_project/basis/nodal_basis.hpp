#pragma once

#include <array>
#include <tuple>
#include <type_traits>

#include <1d_nodal_dg_project/basis/nodal_basis_for_gauss_lobatto.hpp>
#include <1d_nodal_dg_project/quadrature/quadrature_gauss_lobatto.hpp>

namespace ndg {

  /// Checks whether template class is of type Gauss_Lobatto_Quadrature
  template <typename>
  struct Is_Gauss_Lobatto_Quad : std::false_type {};
  template <typename T, unsigned M>
  struct Is_Gauss_Lobatto_Quad<Gauss_Lobatto_Quadrature<T, M>> : std::true_type {};

  /// Wrapper class for the nodal basis
  ///
  /// @tparam T the quadrature type
  template <typename T>
  class Nodal_Basis {
  public:

    /// defines the quadrature type
    using quad_type = T;
    static_assert(Is_Gauss_Lobatto_Quad<quad_type>::value, "Invalid nodes / quadrature type!");

    /// the type of scalar variables
    using scalar_type = typename quad_type::scalar_type;

    /// the number of basis elements per cell
    static constexpr unsigned n_basis_elements = quad_type::get_n_nodes();

    using available_bases_type =
        std::tuple<Nodal_Basis_GL_1<scalar_type>, Nodal_Basis_GL_2<scalar_type>,
                   Nodal_Basis_GL_3<scalar_type>, Nodal_Basis_GL_4<scalar_type>,
                   Nodal_Basis_GL_5<scalar_type>>;

    static_assert((n_basis_elements - 1 >= 0) &&
                      (n_basis_elements - 1 < std::tuple_size_v<available_bases_type>),
                  "No derivatives available for higher order quadrature!");
    using basis_type = typename std::tuple_element_t<n_basis_elements - 1, available_bases_type>;

    /// returns the spatial derivative of the
    /// basis_element_index-th basis element at the node specified by node_index.
    inline static constexpr scalar_type get_x_derivative(const std::size_t & basis_element_index,
                                                         const std::size_t & node_index) {
      return basis_type::get_x_derivatives_at_nodes()[basis_element_index][node_index];
    }

    /// returns the value of the
    /// basis_element_index-th basis element at the node specified by node_index.
    inline static scalar_type get_value(const std::size_t & basis_element_index,
                                        const std::size_t & node_index) {
      return (basis_element_index == node_index);
    }
  };

} // namespace ndg
