#pragma once

#include <array>

namespace ndg {

  /// \brief The nodal basis for Gauss--Lobatto points. Order 1
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct Nodal_Basis_GL_1 {
    static constexpr unsigned n_basis_elements = 1;
    using array_type = std::array<T, n_basis_elements>;
    using matrix_type = std::array<array_type, n_basis_elements>;
    static constexpr matrix_type get_x_derivatives_at_nodes() {
      return {{{0.000000000e+00}}};
    }
  };

  /// \brief The nodal basis for Gauss--Lobatto points. Order 2
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct Nodal_Basis_GL_2 {
    static constexpr unsigned n_basis_elements = 2;
    using array_type = std::array<T, n_basis_elements>;
    using matrix_type = std::array<array_type, n_basis_elements>;
    static constexpr matrix_type get_x_derivatives_at_nodes() {
      return {{{-5.000000000e-01, -5.000000000e-01}, {5.000000000e-01, 5.000000000e-01}}};
    }
  };

  /// \brief The nodal basis for Gauss--Lobatto points. Order 3
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct Nodal_Basis_GL_3 {
    static constexpr unsigned n_basis_elements = 3;
    using array_type = std::array<T, n_basis_elements>;
    using matrix_type = std::array<array_type, n_basis_elements>;

    static constexpr matrix_type get_x_derivatives_at_nodes() {
      return {{
          {-1.500000000e+00, -5.000000000e-01, 5.000000000e-01},
          {2.000000000e+00, 0.0, -2.000000000e+00},
          {-5.000000000e-01, 5.000000000e-01, 1.500000000e+00},
      }};
    }
  };

  /// \brief The nodal basis for Gauss--Lobatto points. Order 4
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct Nodal_Basis_GL_4 {
    static constexpr unsigned n_basis_elements = 4;
    using array_type = std::array<T, n_basis_elements>;
    using matrix_type = std::array<array_type, n_basis_elements>;
    static constexpr matrix_type get_x_derivatives_at_nodes() {
      return {{{-3.000000000e+00, -8.090169944e-01, 3.090169944e-01, -5.000000000e-01},
               {4.045084972e+00, 0.0, -1.118033989e+00, 1.545084972e+00},
               {-1.545084972e+00, 1.118033989e+00, 0.0, -4.045084972e+00},
               {5.000000000e-01, -3.090169944e-01, 8.090169944e-01, 3.000000000e+00}}};
    }
  };

  /// \brief The nodal basis for Gauss--Lobatto points. Order 5
  template <typename T>
  ///
  /// @tparam T the scalar type
  struct Nodal_Basis_GL_5 {
    static constexpr unsigned n_basis_elements = 5;
    using array_type = std::array<T, n_basis_elements>;
    using matrix_type = std::array<array_type, n_basis_elements>;
    static constexpr matrix_type get_x_derivatives_at_nodes() {
      return {
          {{-5.000000000e+00, -1.240990253e+00, 3.750000000e-01, -2.590097470e-01, 5.000000000e-01},
           {6.756502489e+00, 0.0, -1.336584578e+00, 7.637626158e-01, -1.410164178e+00},
           {-2.666666667e+00, 1.745743122e+00, 0.0, -1.745743122e+00, 2.666666667e+00},
           {1.410164178e+00, -7.637626158e-01, 1.336584578e+00, 0.0, -6.756502489e+00},
           {-5.000000000e-01, 2.590097470e-01, -3.750000000e-01, 1.240990253e+00,
            5.000000000e+00}}};
    }
  };

} // namespace ndg
