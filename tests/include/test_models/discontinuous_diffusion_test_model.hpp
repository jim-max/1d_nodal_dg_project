#pragma once

#include <cmath>

/// Test model with discontinuous diffusion term
///
/// See "Mathematical Aspects of Discontinuous Galerkin Methods" by Di Pietro & Ern, page 152
class Discontinuous_Diffusion_Test_Model {
public:

  using scalar_type = double;

  static constexpr scalar_type alpha = 0.01;

  /// the diffusion coefficient function
  static scalar_type diffusion_term(const scalar_type & x, const scalar_type & t) {
    if (x <= 0) {
      return alpha;
    }
    return 1.0;
  }

  /// the advection coefficient
  static scalar_type advection_term(const scalar_type & x, const scalar_type & t) {
    return 0.0;
  }

  /// the reactive term, reaction_term() * u
  static scalar_type reaction_term(const scalar_type & x, const scalar_type & t) {
    return 0.0;
  }

  /// the term in front of the time derivative
  static scalar_type time_derivative_term(const scalar_type & x, const scalar_type & t) {
    return 1.0;
  }

  /// the right hand side of the equation
  static scalar_type rhs_term(const scalar_type & x, const scalar_type & t) {
    return 1.0;
  }

  /// specifies the boundary values
  static scalar_type boundary_value(const scalar_type & x, const scalar_type & t) {
    if (x <= 0) {
      return 0.0;
    }
    return 0.0;
  }

  /// returns the exact solution
  static scalar_type exact_solution(const scalar_type & x, const scalar_type & t) {
    if (x <= 0) {
      const scalar_type a1 = -0.5 / alpha;
      const scalar_type b1 = (1.0 + 3.0 * alpha) / (2.0 * alpha * (1.0 + alpha));
      return a1 * (1.0 + x) * (1.0 + x) + b1 * (1.0 + x);
    }
    const scalar_type a2 = -0.5;
    const scalar_type b2 = -(3.0 + alpha) / (2.0 * (1.0 + alpha));
    return a2 * (x - 1.0) * (x - 1.0) + b2 * (x - 1.0);
  }
};
