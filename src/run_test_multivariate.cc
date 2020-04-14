#include "multivariate_normal_distribution.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>


int main(int, char**)
{
  std::vector<double> a = {1,2};
  std::vector<double> c = {
    5, 1.4,
    1.4, 1
  };

  MultivariateNormalDistribution<> mvn{a, c};

  TFile* file = new TFile("test_multivar_norm_dist.root", "recreate");
  TTree* tree = new TTree("rand_vals", "random values sampled from multivariate normal distribution");

  tree->Branch("x", &a[0]);
  tree->Branch("y", &a[1]);

  for (size_t i = 0; i < 10000; ++i) {
    const auto vals = mvn();
    a.assign(vals.cbegin(), vals.cend());
    tree->Fill();
  }

  tree->Write();
  file->Close();

  return 0;
}
