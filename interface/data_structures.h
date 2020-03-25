#ifndef H_DATA_STRUCTURES__
#define H_DATA_STRUCTURES__

#include "constants.h"

#include "TGraphAsymmErrors.h"

#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <utility>

struct CrossSectionData {

  CrossSectionData(std::string line) {
    std::stringstream sline(line);
    sline >> this->ptM;
    sline >> this->xSec;
    sline >> this->uncer;
    sline >> this->K;
    sline >> this->low;
    sline >> this->high;
  }

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

  double ptM;
  double lth;
  double uncer;
  double low;
  double high;
};

using CrossSectionMeasurement = std::vector<CrossSectionData>;
using PolarizationMeasurement = std::vector<PolarizationData>;

template<typename DataType>
std::vector<DataType> readData(const std::string& filename, const double ptMMin=MIN_PTM)
{
  std::vector<DataType> data;
  std::ifstream file(filename);
  std::string line;

  if (!file) {
    std::cerr << "Could not open file: " << filename << std::endl;
    return data;
  }

  while(std::getline(file, line)) {
    DataType point(line);
    if (point.ptM > ptMMin) {
      data.push_back(point);
    }
  }

  return data;
}

template<typename DataType, typename ValFunc>
TGraphAsymmErrors asTGraph(const std::vector<DataType>& data, ValFunc vFunc)
{
  std::vector<double> ptm;
  std::vector<double> val;
  std::vector<double> uncer;
  std::vector<double> low;
  std::vector<double> high;

  for (const auto& point : data) {
    ptm.push_back(point.ptM);
    val.push_back(vFunc(point));
    uncer.push_back(point.uncer);
    low.push_back(point.ptM - point.low);
    high.push_back(point.high - point.ptM);
  }

  return TGraphAsymmErrors(data.size(), ptm.data(), val.data(), low.data(), high.data(), uncer.data(), uncer.data());
}


TGraphAsymmErrors asTGraph(const CrossSectionMeasurement& data)
{
  return asTGraph(data, [](const CrossSectionData& point) { return point.xSec;} );
}

TGraphAsymmErrors asTGraph(const PolarizationMeasurement& data)
{
  return asTGraph(data, [](const PolarizationData& point) { return point.lth;} );
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


#endif
