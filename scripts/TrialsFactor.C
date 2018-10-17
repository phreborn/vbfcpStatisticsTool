// Implementation of Gross, E. & Vitells, "Trial factors for the look elsewhere effect in high energy physics",
// O. Eur. Phys. J. C (2010) 70: 525. https://doi.org/10.1140/epjc/s10052-010-1470-8
//
// A.L. Read June, 2012
//
// Print the local and global p-value and trials factor given the average number of up-crossings, the threshold for
// counting the up-crossings (number of standard deviations) and the observed local significance (number of standard
// deviations) for a single degree of freedom. See Equation 3 in the article.
//
// A graph showing the trials factor (minus 1) per up-crossing as a function of local significance for various
// thresholds for the up-crossings is also shown. For example, an observed local significance of 5 and 2 up-crossings
// at 0.5 sigma threshold gives a trials factor (minus 1) of ~30 (29.46).
  
void TrialsFactor(Double_t crossings=1, Double_t crossing_significance=0, Double_t local_significance=2.5){
  TF1 *f00 = new TF1("Trials factor","exp(-(pow(x,2)-pow(0.0,2))/2)/(ROOT::Math::chisquared_cdf_c(pow(x,2),1)/2)",2.0,8);
  TF1 *f05 = new TF1("Trials factor","exp(-(pow(x,2)-pow(0.5,2))/2)/(ROOT::Math::chisquared_cdf_c(pow(x,2),1)/2)",2.5,8);
  TF1 *f10 = new TF1("Trials factor","exp(-(pow(x,2)-pow(1.0,2))/2)/(ROOT::Math::chisquared_cdf_c(pow(x,2),1)/2)",3.0,8);
  TF1 *f15 = new TF1("Trials factor","exp(-(pow(x,2)-pow(1.5,2))/2)/(ROOT::Math::chisquared_cdf_c(pow(x,2),1)/2)",3.5,8);
 
  f00->SetMaximum(50);
  f00->SetMinimum(5);

  f00->SetLineColor(1);
  f05->SetLineColor(2);
  f10->SetLineColor(3);
  f15->SetLineColor(4);
 
  TCanvas *c1 = new TCanvas();
  c1->SetGrid();

  f00->Draw();
  f05->Draw("same");
  f10->Draw("same");
  f15->Draw("same");
 
  TLegend *legend = new TLegend(0.65,0.45,0.85,0.85);
  legend->AddEntry(f00,"Per 0.0 #sigma crossing","l");
  legend->AddEntry(f05,"Per 0.5 #sigma crossing","l");
  legend->AddEntry(f10,"Per 1.0 #sigma crossing","l");
  legend->AddEntry(f15,"Per 1.5 #sigma crossing","l");
  legend->Draw();

   f00->SetTitle("p_{0}^{global} = p_{0}^{local} + <N>*k, k=(Trials factor-1)/crossing;Local signficance (#sigma);(Trials factor-1)/crossing");
  
  Double_t local_p0 = ROOT::Math::chisquared_cdf_c(pow(local_significance,2),1)/2;
  Double_t global_p0 = local_p0 + crossings * exp(-(pow(local_significance,2)-pow(crossing_significance,2))/2);
  Double_t trials_factor = global_p0/local_p0;
  Double_t global_significance = sqrt(ROOT::Math::chisquared_quantile_c(global_p0*2,1));
  cout << "Local  p-value, significance: " <<  local_p0 << ", " << local_significance << " sigma" << endl;
  cout << "Trials factor: " << trials_factor << " (Trials factor minus 1: " << trials_factor-1 << ")" << endl;
  cout << "Global p-value, significance: " <<   global_p0 << ", " << global_significance << " sigma" << endl;
}
