#ifndef H_FITTER__
#define H_FITTER__

#include "Math/Minimizer.h"
#include "Fit/Fitter.h"
#include "Math/MinimizerOptions.h"
#include "Fit/FitResult.h"

class LikelihoodFitter {
public:
  LikelihoodFitter() {
    ROOT::Math::MinimizerOptions minOpt;
    minOpt.SetMinimizerType("Minuit2");
    minOpt.SetPrintLevel(2);
    minOpt.SetDefaultErrorDef(0.5); // likelihood fit
    fitter.Config().SetMinimizerOptions(minOpt);
    // fitter.Config().SetMinosErrors(true);
  }


  template<typename LLH>
  bool Fit(LLH& llh);

  template<typename LLH>
  TGraph Contour(const LLH& llh, const char* par1, const char* par2, double errLvl=0.5, unsigned int nPoints=25) const;

  const ROOT::Fit::FitResult& Result() const { return fitter.Result(); }

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
