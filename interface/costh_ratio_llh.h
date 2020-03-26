#ifndef H_COSTH_RATIO_LLH__
#define H_COSTH_RATIO_LLH__

#include "data_structures.h"
#include "simple_compile_time_map.h"

#include <vector>
#include <array>

constexpr std::array<ParameterIndex, 5> COSTH_RATIO_PARAMS = {{
    {"lambda_1", 0},
    {"lambda_2", 1},
    {"norm_1", 2},
    {"norm_2", 3},
    {"norm_3", 4}
  }};

static constexpr int CPAR(const char* name) {
  return getParIdx(COSTH_RATIO_PARAMS, name);
}

class CosthRatioLikelihood {
public:
  CosthRatioLikelihood(const std::vector<CosthRatioMeasurement>& data) : m_costh_ratios(data)
  {
    defineStartParams();
  }

  double operator()(const double* p) const;

  size_t nPars() const { return COSTH_RATIO_PARAMS.size(); }

  int getParIdx(const std::string& name) const { return CPAR(name.c_str()); }

  ParamsSettings getStartParams() const { return m_startParams; }

private:
  void defineStartParams();

  template<class... Args>
  void setParam(const char* name, Args&&... args)
  { m_startParams[CPAR(name)] = ROOT::Fit::ParameterSettings(name, args...); }

  std::vector<CosthRatioMeasurement> m_costh_ratios;
  ParamsSettings m_startParams{mapSize(COSTH_RATIO_PARAMS)};
};


double costhRatio(double costh, double lamD, double lamN, double norm)
{
  return norm * (1 + lamN * costh * costh) / (1 + lamD * costh * costh);
}

double loglikeCosthRatio(const CosthRatioMeasurement& data,
                         const double lamD, const double lamN, const double norm)
{
  double loglike = 0;

  for (const auto& point : data) {
    const double rcosth = costhRatio(point.costh, lamD, lamN, norm);
    const double uncer = (rcosth > point.ratio) ? point.u_high : point.u_low;
    const double relDiff = (rcosth - point.ratio) / uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return -loglike;
}

double CosthRatioLikelihood::operator()(const double* p) const
{
  const double lambda1 = p[CPAR("lambda_1")];
  const double lambda2 = p[CPAR("lambda_2")];
  const double norm1 = p[CPAR("norm_1")];
  const double norm2 = p[CPAR("norm_2")];
  const double norm3 = p[CPAR("norm_3")];

  double loglike = 0;

  loglike += loglikeCosthRatio(m_costh_ratios[0], lambda1, lambda2, norm1);
  loglike += loglikeCosthRatio(m_costh_ratios[1], lambda1, lambda2, norm2);
  loglike += loglikeCosthRatio(m_costh_ratios[2], lambda1, lambda2, norm3);

  return loglike;
}

void CosthRatioLikelihood::defineStartParams()
{
  setParam("lambda_1", 0, 0.5);
  setParam("lambda_2", 0, 0.5);
  setParam("norm_1", 0.45, 0.2);
  setParam("norm_2", 0.45, 0.2);
  setParam("norm_3", 0.45, 0.2);
}


#endif


