#pragma once

#include <cmath>

/// Test model with linear exact solution
///
/// See "Mathematical Aspects of Discontinuous Galerkin Methods" by Di Pietro & Ern, page 152
class Linear_Test_Model {
public:

  using scalar_type = double;

  static scalar_type diffusion_term(const scalar_type & x, const scalar_type & t) {
    return 1.0;
  }

  static scalar_type advection_term(const scalar_type & x, const scalar_type & t) {
    return 0.0;
  }

  static scalar_type reaction_term(const scalar_type & x, const scalar_type & t) {
    return 0.0;
  }

  static scalar_type time_derivative_term(const scalar_type & x, const scalar_type & t) {
    return 1.0;
  }

  static scalar_type rhs_term(const scalar_type & x, const scalar_type & t) {
    return 0.0;
  }

  static scalar_type boundary_value(const scalar_type & x, const scalar_type & t) {
    if (x <= 0) {
      return 1.0;
    }
    return -1.0;
  }

  /// the exact solution to this test problem
  static scalar_type exact_solution(const scalar_type & x, const scalar_type & t) {
    return -x;
  }
};
