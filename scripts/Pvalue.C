// Return single-sided pvalue for equivalent significance (standard deviations)
//
// A.L. Read June 2012

Double_t Pvalue(Double_t significance) {
  return ROOT::Math::chisquared_cdf_c(pow(significance,2),1)/2;
}
