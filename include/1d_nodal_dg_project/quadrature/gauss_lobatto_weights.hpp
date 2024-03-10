#pragma once

#include <array>

namespace ndg {

  /// \brief Gauss--Lobatto quadrature points and weights of order 1
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct GL_Quad_1 {
    static constexpr unsigned n_nodes = 1;
    using array_type = std::array<T, n_nodes>;

    constexpr static array_type get_nodes() {
      return {{0.0}};
    }
    constexpr static array_type get_weights() {
      return {{2.0}};
    }
  };

  /// \brief Gauss--Lobatto quadrature points and weights of order 2
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct GL_Quad_2 {
    static constexpr unsigned n_nodes = 2;
    using array_type = std::array<T, n_nodes>;

    constexpr static array_type get_nodes() {
      return {{-1.0, 1.0}};
    }
    constexpr static array_type get_weights() {
      return {{1.0, 1.0}};
    }
  };

  /// \brief Gauss--Lobatto quadrature points and weights of order 3
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct GL_Quad_3 {
    static constexpr unsigned n_nodes = 3;
    using array_type = std::array<T, n_nodes>;

    constexpr static array_type get_nodes() {
      return {{-1.0, 0.0, 1.0}};
    }
    constexpr static array_type get_weights() {
      return {{0.33333333333333, 1.33333333333333, 0.33333333333333}};
    }
  };

  /// \brief Gauss--Lobatto quadrature points and weights of order 4
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct GL_Quad_4 {
    static constexpr unsigned n_nodes = 4;
    using array_type = std::array<T, n_nodes>;

    constexpr static array_type get_nodes() {
      return {{-1.0, -0.4472135955, 0.4472135955, 1.0}};
    }
    constexpr static array_type get_weights() {
      return {{0.16666666666666, 0.83333333333333, 0.83333333333333, 0.16666666666666}};
    }
  };

  /// \brief Gauss--Lobatto quadrature points and weights of order 5
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct GL_Quad_5 {
    static constexpr unsigned n_nodes = 5;
    using array_type = std::array<T, n_nodes>;

    constexpr static array_type get_nodes() {
      return {{-1.0, -0.6546536707, 0, 0.6546536707, 1.0}};
    }
    constexpr static array_type get_weights() {
      return {{0.1, 0.54444444444444, 0.71111111111111, 0.54444444444444, 0.1}};
    }
  };

  /// \brief Gauss--Lobatto quadrature points and weights of order 6
  ///
  /// @tparam T the scalar type
  template <typename T>
  struct GL_Quad_6 {
    static constexpr unsigned n_nodes = 6;
    using array_type = std::array<T, n_nodes>;

    constexpr static array_type get_nodes() {
      return {{-1.0, -0.76505532, -0.28523152, 0.28523152, -0.76505532, 1.0}};
    }
    constexpr static array_type get_weights() {
      return {{0.066666666, 0.37847496, 0.55485838, 0.55485838, 0.37847496, 0.066666666}};
    }
  };

} // namespace ndg
