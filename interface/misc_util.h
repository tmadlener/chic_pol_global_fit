#ifndef GLOBALFIT_MISC_UTILS_H
#define GLOBALFIT_MISC_UTILS_H

#include <type_traits>
#include <vector>
#include <cmath>

namespace util {
inline double interpolate(double x, double x0, double x1, double y0, double y1) {
  return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

/**
 * Linearly interpolate in the logarithmic space in the y values
 */
inline double log_interpolate(double x, double x0, double x1, double y0, double y1) {
  const double logy0 = std::log(y0);
  const double logy1 = std::log(y1);

  return std::exp(logy0 + (x - x0) * (logy1 - logy0) / (x1 - x0));
}

/**
 * Find the bin index corresponding to the given x value
 */
inline int find_bin(double x, const std::vector<double>& points) {
  for (size_t i = 0; i < points.size(); ++i) {
    if (x >= points[i] && x <= points[i + 1]) { return i; }
  }
  return -1; // fail somewhere but not here
}

/**
 * Convert an enum (class) value to something that is usable as an array index
 */
template <typename E>
constexpr auto to_index(E e) noexcept {
  return static_cast<std::underlying_type_t<E>>(e);
}

template<typename E>
constexpr auto to_enum(int index) noexcept {
  return static_cast<E>(index);
}

} // namespace util

#endif
