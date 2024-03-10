#pragma once

namespace ndg {

  /// \brief append vector from_vecv to to_vec
  ///
  /// @param to_vec
  /// @param from_vec
  /// @tparam T the vector class, compatible to std::vector
  template <typename T>
  inline void vec_append(T & to_vec, const T & from_vec) {
    to_vec.reserve(to_vec.size() + from_vec.size());
    to_vec.insert(to_vec.end(), from_vec.begin(), from_vec.end());
  }

  /// \brief compute the weighted average of two values
  ///
  /// @tparam T the scalar type
  template <class T>
  inline T weighted_mean(const T & value_u, const T & value_v, const T & weight_u,
                         const T & weight_v) {
    return (weight_v * value_u + weight_u * value_v) / (weight_u + weight_v);
  }

  /// \brief compute the harmonic average of two values
  ///
  /// @tparam T the scalar type
  template <class T>
  inline T harmonic_mean(const T & value_u, const T & value_v, const T & weight_u,
                         const T & weight_v) {
    return weight_u * weight_v * (value_u + value_v) / (weight_u + weight_v);
  }

} // namespace ndg
