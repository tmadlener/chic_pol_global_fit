#ifndef H_FITTER__
#define H_FITTER__

#include "Math/Minimizer.h"
#include "Fit/Fitter.h"
#include "Math/MinimizerOptions.h"

class LikelihoodFitter {
public:
  LikelihoodFitter() {
    ROOT::Math::MinimizerOptions minOpt;
    minOpt.SetMinimizerType("Minuit");
    minOpt.SetPrintLevel(2);

    fitter.Config().SetMinimizerOptions(minOpt);
  }


  template<typename LLH>
  bool Fit(LLH& llh);

private:
  ROOT::Fit::Fitter fitter{};
};

template<typename LLH>
bool LikelihoodFitter::Fit(LLH &llh)
{
  const int nPar = llh.nPars();

  fitter.Config().SetParamsSettings(llh.getStartParams());

  const auto ret = fitter.FitFCN(nPar, llh); // default arguments are the ones
                                             // we want in this case

  fitter.GetMinimizer()->PrintResults();
  return ret;
}



#endif
