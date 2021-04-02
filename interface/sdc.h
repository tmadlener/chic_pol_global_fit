#ifndef GLOBAL_FIT_SDC_H
#define GLOBAL_FIT_SDC_H

#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace sdc {
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
   * Check whether the two vectors describe the same support points
   */
  bool compare_supports(const std::vector<double>& lhs, const std::vector<double>& rhs) {
    if (lhs.size() != rhs.size()) {
      std::cerr << "WARNING: different support sizes: " << lhs.size() << " vs " << rhs.size() << std::endl;
      return false;
    }

    for (size_t i = 0; i < lhs.size(); ++i) {
      if (lhs[i] != rhs[i]) {
        std::cerr << "WARNING: suppport point " << i << " differs: " << lhs[i] << " vs " << rhs[i] << std::endl;
        return false;
      }
    }
    return true;
  }

  /**
   * Sum the sdcs at the given support point using the corresponding ldmes
   */
  double sum_sdcs(const std::vector<std::vector<double>>& sdcs, const std::vector<double>& ldmes, size_t bin) {
    double sum = 0;
    for (size_t i = 0; i < sdcs.size(); ++i) {
      sum += sdcs[i][bin] * ldmes[i];
    }

    return sum;
  }

} // namespace util

class SDC {
  friend class MultiSDC;

public:
  /**
   * The one (and only) constructor from two vectors: the support points and
   * the values at these support points
   */
  explicit SDC(std::vector<double> supports, std::vector<double> values) :
      m_supports(supports), m_values(values), m_size(supports.size()) {
    if (m_size != m_values.size()) {
      // TODO: Fail?
      std::cerr << "WARNING: different numer of support points and values: " << m_size << " vs " << m_values.size()
                << std::endl;
    }
  }

  // Rule of 5
  SDC() = delete;
  SDC(const SDC&) = default;
  SDC(SDC&&) = default;
  SDC& operator=(const SDC&) = default;
  SDC& operator=(SDC&&) = default;
  ~SDC() = default;

  /**
   * Evaluate at a given ptM value. Use the second argument to control how this
   * interpolation should be done, by default a simple linear interpolation
   * between the two adjacent support points and values is done.
   *
   * The second argument has to be something that can be used as a
   * std::function<double(double,double,double,double,double)>), i.e. a function
   * that takes 5 double arguments, and returns a double. The arguments are: x,
   * x0, x1, y0, y1. These are the values that can be used in the interpolation.
   */
  template <typename InterpolationFunc = decltype(util::interpolate)>
  double operator()(double ptM, InterpolationFunc&& interpFunc = util::interpolate) const {
    const int i = util::find_bin(ptM, m_supports);
    if (i >= 0 && (unsigned)i < m_size) {
      return interpFunc(ptM, m_supports[i], m_supports[i + 1], m_values[i], m_values[i + 1]);
    }
    return std::nan("");
  }

  /**
   * stream an SDC as two columns: pt value
   */
  friend std::ostream& operator<<(std::ostream& os, const SDC& sdc);

  // /**
  //  * Scale the SDC values by a constant factor
  //  */
  // friend SDC operator*(double factor, const SDC& sdc);

  /**
   * Many different functions involving a scalar can be implemented generically,
   * they only differ in the way the values are combined
   */
  template <typename T, typename BinaryValueOp>
  friend SDC scalar_op(const T val, const SDC& rhs, BinaryValueOp&& binOP);

  /**
   * Many of the functions above can be implemented generically, they only
   * differ in the way the values are combined at each support point
   */
  template <typename BinaryValueOP>
  friend SDC binary_op(const SDC& lhs, const SDC& rhs, BinaryValueOP&& binOP);

  /**
   * Scale the support points
   */
  void scale_support(double factor) {
    for (auto& s : m_supports)
      s *= factor;
  }

private:
  std::vector<double> m_supports; ///< points where SDC values are provided
  std::vector<double> m_values;   ///< SDC values at the given support points
  size_t m_size;                  ///< The number of points (minor optimization)
};

/**
 * A combination of multiple SDCs that can be evaluated at a given ptM value if
 * additionally the corresponding LDMEs are provided.
 */
class MultiSDC {
public:
  /**
   * One and only constructor taking a collection of SDCs that all have to
   * have the same supports
   */
  explicit MultiSDC(const std::vector<SDC>& sdcs) : m_supports(sdcs[0].m_supports), m_size(sdcs[0].m_size) {
    m_values.reserve(sdcs.size());
    for (const auto& s : sdcs) {
      if (!util::compare_supports(m_supports, s.m_supports)) continue;
      m_values.emplace_back(s.m_values);
    }
  }

  // Rule of 5
  MultiSDC() = delete;
  MultiSDC(const MultiSDC&) = default;
  MultiSDC(MultiSDC&&) = default;
  MultiSDC& operator=(const MultiSDC&) = default;
  MultiSDC& operator=(MultiSDC&&) = default;

  /**
   * Evaluate the combination of SDCs at the ptM value using the corresponding
   * LDMEs and an interpolation function to interpolate between the supports
   * where the SDCs are given. Default is linear interpolation.
   *
   * To evaluate the combination of SDCs first the values of all SDCs are
   * multiplied with the corresponding LDMEs. These values are then summed and
   * the interpolation is done using the supports and the summed SDCs at the
   * supports. In this way it is possible to deal with negative SDCs easily as,
   * at least for physically meaningful combinations of SDCs, their sum should
   * again be positive.
   */
  template <typename InterpolationFunc = decltype(util::interpolate)>
  double operator()(double ptM, const std::vector<double>& ldmes,
                    InterpolationFunc interpFunc = util::interpolate) const {
    const auto i = util::find_bin(ptM, m_supports);
    if (i >= 0 && (unsigned)i < m_size) {
      return interpFunc(ptM, m_supports[i], m_supports[i + 1], util::sum_sdcs(m_values, ldmes, i),
                        util::sum_sdcs(m_values, ldmes, i + 1));
    }
    return std::nan("");
  }

  /**
   * Scale the support points
   */
  void scale_support(double factor) {
    for (auto& s : m_supports)
      s *= factor;
  }

  /**
   * Stream the MultiSDC: pt SDC0 SDC1 ...
   */
  friend std::ostream& operator<<(std::ostream&, const MultiSDC&);

private:
  std::vector<double> m_supports;
  const size_t m_size;
  std::vector<std::vector<double>> m_values{};
};

/**
 * Add two SDCs (assuming they have the same supports!)
 */
SDC operator+(const SDC& lhs, const SDC& rhs) { return binary_op(lhs, rhs, std::plus<double>{}); }

/**
 * Subtract rhs SDC from lhs SDC (assuming they have the same support)
 */
SDC operator-(const SDC& lhs, const SDC& rhs) { return binary_op(lhs, rhs, std::minus<double>{}); }

/**
 * Multiply two SDC (values), assuming they have the same supports.
 *
 * NOTE: not strictly necessary for SDCs themselves, but allows to simplify some
 * of the manipulations we do at the start
 */
SDC operator*(const SDC& lhs, const SDC& rhs) { return binary_op(lhs, rhs, std::multiplies<double>{}); }

/**
 * Divide the numerator SDC (values) by the denominator SDC (values), given the
 * two have the same support
 */
SDC operator/(const SDC& num, const SDC& denom) { return binary_op(num, denom, std::divides<double>{}); }

/**
 * Scale an SDC by a constant factor
 */
SDC operator*(const double scale, const SDC& sdc) { return scalar_op(scale, sdc, std::multiplies<double>{}); }

/**
 * Shift an SDC by a constant factor
 */
SDC operator+(const double shift, const SDC& sdc) { return scalar_op(shift, sdc, std::plus<double>{}); }

/**
 * Subtract an SDC from a constant factor
 */
SDC operator-(const double shift, const SDC& sdc) { return scalar_op(shift, sdc, std::minus<double>{}); }

std::ostream& operator<<(std::ostream& os, const SDC& sdc) {
  for (size_t i = 0; i < sdc.m_size; ++i) {
    os << sdc.m_supports[i] << " " << sdc.m_values[i] << '\n';
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const MultiSDC& sdc) {
  for (size_t i = 0; i < sdc.m_size; ++i) {
    os << sdc.m_supports[i];
    for (const auto& vals : sdc.m_values) {
      os << " " << vals[i];
    }
    os << '\n';
  }
  return os;
}

template <typename T, typename BinaryValueOp>
SDC scalar_op(const T val, const SDC& sdc, BinaryValueOp&& binOP) {
  static_assert(std::is_arithmetic_v<T>, "Can only do scalar operations on SDCs with arithmetic types");

  auto values = sdc.m_values;
  for (auto& v : values)
    binOP(val, v);

  return SDC(sdc.m_supports, std::move(values));
}

template <typename BinaryValueOP>
SDC binary_op(const SDC& lhs, const SDC& rhs, BinaryValueOP&& binOP) {
  if (!util::compare_supports(lhs.m_supports, rhs.m_supports)) {
    std::cerr << "ERROR: cannot combine SDCs" << std::endl;
    return SDC({}, {});
  }

  auto values = lhs.m_values;
  for (size_t i = 0; i < lhs.m_size; ++i) {
    values[i] = binOP(values[i], rhs.m_values[i]);
  }

  return SDC(lhs.m_supports, std::move(values));
}

/**
 * Different types of SDCs that are availalbe from the input SDC files
 */
enum class SDCType {
  LP_NLO = 0, ///< LP + NLO is the first column of values (after pt)
  LO = 1,     ///< LO is the second column
  NLO = 2,    ///< NLO is the third column
};

constexpr std::array<const char*, 3> SDCTypeNames = {"LP+NLO", "LO", "NLO"};

namespace util {
  /**
   * Convert an enum (class) value to something that is usable as an array index
   */
  template <typename E>
  constexpr auto to_index(E e) noexcept {
    return static_cast<std::underlying_type_t<E>>(e);
  }
} // namespace util

/**
 * Read an SDC (LP + NLO) from a file
 */
SDC read_from_file(std::string filename, SDCType sdcType = SDCType::LP_NLO) {
  using util::to_index;

  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << "ERROR: Cannot open file \'" << filename << "\' to read SDC" << std::endl;
    std::exit(1); // TODO: Does this have to be fatal?
  }

  std::cout << "Reading SDC from file \'" << filename << "\' using type \'" << SDCTypeNames[to_index(sdcType)]
            << "\' ... ";

  std::string line;
  std::vector<double> points, values;
  while (std::getline(infile, line)) {
    // Ignore comments
    if (line[0] == '#') continue;

    std::stringstream linestr{line};
    // Assuming format: pT LP+NLP LO NLO
    double p;
    std::array<double, 3> vals;
    linestr >> p >> vals[to_index(SDCType::LP_NLO)] >> vals[to_index(SDCType::LO)] >> vals[to_index(SDCType::NLO)];

    points.push_back(p);
    values.push_back(vals[to_index(sdcType)]);
  }

  std::cout << "(" << points.size() << " value pairs)\n";

  return SDC(std::move(points), std::move(values));
}

/**
 * Helper struct for easier passing of SDCs
 */
struct StateSDCs {
  MultiSDC tot;
  MultiSDC lng; // Cannot name it long because that is a builtin
};
} // namespace sdc

#endif
