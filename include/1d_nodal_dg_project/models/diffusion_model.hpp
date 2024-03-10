#pragma once

#include <cmath>

namespace ndg {

  /// General diffusion model
  class Diffusion_Model {
  public:

    using scalar_type = double;

    /// default costructor
    Diffusion_Model() : alpha_l(1), alpha_r(0.1) {}

    /// constructor, setting the diffusion coefficients
    Diffusion_Model(const scalar_type & alpha_l_, const scalar_type & alpha_r_)
        : alpha_l(alpha_l_), alpha_r(alpha_r_) {}

    /// the diffusion coefficient function
    scalar_type diffusion_term([[maybe_unused]] const scalar_type & x,
                               [[maybe_unused]] const scalar_type & t) const {
      return (x <= 0) ? alpha_l : alpha_r;
    }

    /// the advection coefficient
    scalar_type advection_term([[maybe_unused]] const scalar_type & x,
                               [[maybe_unused]] const scalar_type & t) const {
      return 0.0;
    }

    /// the reactive term, reaction_term() * u
    scalar_type reaction_term([[maybe_unused]] const scalar_type & x,
                              [[maybe_unused]] const scalar_type & t) const {
      return 0.0;
    }

    /// the term in front of the time derivative
    scalar_type time_derivative_term([[maybe_unused]] const scalar_type & x,
                                     [[maybe_unused]] const scalar_type & t) const {
      return 1.0;
    }

    /// the right hand side of the equation
    scalar_type rhs_term([[maybe_unused]] const scalar_type & x,
                         [[maybe_unused]] const scalar_type & t) const {
      return 0.0;
    }

    /// specifies the boundary values
    scalar_type boundary_value([[maybe_unused]] const scalar_type & x,
                               [[maybe_unused]] const scalar_type & t) const {
      if (x <= 0) {
        return 1.0;
      }
      return 0.0;
    }

    /// returns the exact solution,
    /// here just a dummy function
    scalar_type exact_solution([[maybe_unused]] const scalar_type & x,
                               [[maybe_unused]] const scalar_type & t) const {
      return 0;
    }

  private:

    /// diffusion coefficient for x <= 0
    const scalar_type alpha_l;
    /// diffusion coefficient for x > 0
    const scalar_type alpha_r;
  };

} // namespace ndg
