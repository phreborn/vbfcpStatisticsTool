// Return equivalent significance (standard deviations) for single-sided pvalue.
//
// A.L. Read June 2012

Double_t Significance(Double_t pval) {
  return sqrt(ROOT::Math::chisquared_quantile_c(pval*2,1));
}
