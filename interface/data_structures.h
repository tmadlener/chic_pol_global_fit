#ifndef H_DATA_STRUCTURES__
#define H_DATA_STRUCTURES__

#include "constants.h"

#include "TGraphAsymmErrors.h"
#include "Fit/ParameterSettings.h"

#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <utility>

struct CrossSectionData {

  CrossSectionData(const double ptm, const double lo, const double hi) :
    ptM(ptm), xSec(0), uncer(0), K(0), low(lo), high(hi) {}

  CrossSectionData(std::string line) {
    std::stringstream sline(line);
    sline >> this->ptM;
    sline >> this->xSec;
    sline >> this->uncer;
    sline >> this->K;
    sline >> this->low;
    sline >> this->high;
  }

  CrossSectionData(const CrossSectionData&) = default;
  CrossSectionData(CrossSectionData&&) = default;
  CrossSectionData& operator=(const CrossSectionData&) = default;
  CrossSectionData& operator=(CrossSectionData&&) = default;
  ~CrossSectionData() = default;

  double ptM;
  double xSec;
  double uncer;
  double K;
  double low;
  double high;
};

struct PolarizationData {

  PolarizationData(std::string line) {
    std::stringstream sline(line);

    sline >> this->ptM;
    sline >> this->lth;
    sline >> this->uncer;
    sline >> this->low;
    sline >> this->high;
  }

  PolarizationData(const PolarizationData&) = default;
  PolarizationData(PolarizationData&&) = default;
  PolarizationData& operator=(const PolarizationData&) = default;
  PolarizationData& operator=(PolarizationData&&) = default;
  ~PolarizationData() = default;

  double ptM;
  double lth;
  double uncer;
  double low;
  double high;
};

using CrossSectionMeasurement = std::vector<CrossSectionData>;
using PolarizationMeasurement = std::vector<PolarizationData>;

template<typename DataType>
std::vector<DataType> readData(const std::string& filename)
{
  std::vector<DataType> data;
  std::ifstream file(filename);
  std::string line;

  if (!file) {
    std::cerr << "Could not open file: " << filename << std::endl;
    return data;
  }

  std::cout << "Reading data from: \'" << filename << "\' ... ";

  while(std::getline(file, line)) {
    data.emplace_back(line);
  }

  std::cout << "(" << data.size() << " data points)\n";

  return data;
}

template<typename DataType>
std::vector<DataType> readData(const std::string& filename, const double minPtM)
{
  auto data = readData<DataType>(filename);
  const auto sizePre = data.size();
  data.erase(std::remove_if(data.begin(), data.end(),
                            [minPtM](const DataType& p){ return p.ptM < minPtM; }),
             data.end());

  if (data.size() != sizePre) {
    std::cerr << "INFO: Removed " << sizePre - data.size() << " data points due to min pT/M ("
              << minPtM << ") requirement" << std::endl;
  }

  return data;
}

template<typename DataType>
double symError(const DataType& point) { return point.uncer; }

template<typename DataType>
double ptM(const DataType& point) { return point.ptM; }

template<typename DataType>
using PFunc = std::function<double(const DataType&)>;


/**
 * This is the all purpose function that does the heavy lifting of creating a
 * TGraph from the internal data structures. In principle it can process all
 * vectors of DataTypes where the DataType provides a low and high member that
 * gives the lower and upper bounds of the bin. All other function retrieval
 * values are parametrized and expect a function taking a const DataType& as
 * single argument and returing a double. The defaults are such that the
 * CrossSectionData and PolarizationData types can be processed only defining a
 * function that gives the value on the y-axis.
 */
template<typename DataType, typename ValF=PFunc<DataType>, typename XF=PFunc<DataType>,
         typename ErrLoF=PFunc<DataType>, typename ErrHiF=PFunc<DataType>>
TGraphAsymmErrors asTGraph(const std::vector<DataType>& data, ValF vFunc, XF xFunc=ptM<DataType>,
                           ErrLoF errLoF=symError<DataType>, ErrHiF errHiF=symError<DataType>)
{
  std::vector<double> ptm;
  std::vector<double> val;
  std::vector<double> uhi;
  std::vector<double> ulow;
  std::vector<double> low;
  std::vector<double> high;

  for (const auto& point : data) {
    const double xval = xFunc(point);
    ptm.push_back(xval);
    val.push_back(vFunc(point));
    ulow.push_back(errLoF(point));
    uhi.push_back(errHiF(point));
    low.push_back(xval - point.low);
    high.push_back(point.high - xval);
  }

  return TGraphAsymmErrors(data.size(), ptm.data(), val.data(), low.data(), high.data(), ulow.data(), uhi.data());
}


/**
 * Overload for CrossSectionData
 */
TGraphAsymmErrors asTGraph(const CrossSectionMeasurement& data)
{
  return asTGraph(data, [](const CrossSectionData& point) { return point.xSec; } );
}

/**
 * Overload for PolarizationData
 */
TGraphAsymmErrors asTGraph(const PolarizationMeasurement& data)
{
  return asTGraph(data, [](const PolarizationData& point) { return point.lth; } );
}

struct NuissanceParameter {
public:
  NuissanceParameter(const std::array<double, 2>& par) : central(par[0]), uncer(par[1]) {}
  NuissanceParameter(const double unc) : central(1.0), uncer(unc) {}

  double operator()(double val) const
  {
    const double sigmas = (val - 1) / (uncer / central);
    return -0.5 * sigmas * sigmas;
  }

  double central;
  double uncer;
};

struct ParameterScanSettings {
  std::vector<double> values;
  std::string name;
};

using ScanSettings = std::pair<ParameterScanSettings, ParameterScanSettings>;


struct CosthRatioData {
  CosthRatioData(const std::string& line) {
    std::stringstream sstr(line);
    sstr >> this->costh;
    sstr >> this->ratio;
    sstr >> this->low;
    sstr >> this->high;
    sstr >> this->u_low;
    sstr >> this->u_high;
  }

  CosthRatioData(const CosthRatioData&) = default;
  CosthRatioData(CosthRatioData&&) = default;
  CosthRatioData& operator=(const CosthRatioData&) = default;
  CosthRatioData& operator=(CosthRatioData&&) = default;
  ~CosthRatioData() = default;


  double costh;
  double low;
  double high;
  double u_low;
  double u_high;
  double ratio;
};

using CosthRatioMeasurement = std::vector<CosthRatioData>;

struct BinInfo {
  double ptM;
  double low;
  double high;
};

using PtCosthRatioMeasurement = std::pair<BinInfo, CosthRatioMeasurement>;

/**
 * Overload for CosthRatioData
 */
TGraphAsymmErrors asTGraph(const CosthRatioMeasurement& data)
{
  return asTGraph(data,
                  [] (const CosthRatioData& p) { return p.ratio; },
                  [] (const CosthRatioData& p) { return p.costh; },
                  [] (const CosthRatioData& p) { return p.u_low; },
                  [] (const CosthRatioData& p) { return p.u_high; });
}


using ParamsSettings = std::vector<ROOT::Fit::ParameterSettings>;

/**
 * Helper struct to hold all the global fit data, to make it easier to pass
 * around
 */
struct GlobalFitData {
  CrossSectionMeasurement psi2S_ATLAS_cs;
  CrossSectionMeasurement psi2S_CMS_cs;
  CrossSectionMeasurement jpsi_CMS_cs;
  CrossSectionMeasurement chic2_ATLAS_cs;
  CrossSectionMeasurement chic1_ATLAS_cs;
  CrossSectionMeasurement chic_ratio_CMS_cs;

  PolarizationMeasurement psi2S_CMS_pol;
  PolarizationMeasurement jpsi_CMS_pol;

  std::vector<PtCosthRatioMeasurement> chic_costh_ratios_CMS;
};


#endif
