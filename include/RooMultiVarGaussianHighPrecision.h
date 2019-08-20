/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_MULTI_VAR_GAUSSIAN_HP
#define ROO_MULTI_VAR_GAUSSIAN_HP

#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"

class RooRealVar;
class RooFitResult ;

#include <map>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/io.hpp> 
#include <boost/numeric/ublas/matrix_proxy.hpp> 
#include <boost/numeric/ublas/lu.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using boost::multiprecision::cpp_dec_float_50;

#if defined __GNUC__ || defined __APPLE__
#include <Eigen/Dense>
#else
#include <eigen3/Eigen/Dense>
#endif

typedef boost::multiprecision::cpp_dec_float<50> mp_backend;
typedef boost::multiprecision::number<mp_backend, boost::multiprecision::et_off> SuperFloat;
typedef Eigen::Matrix<cpp_dec_float_50, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixR;

class RooMultiVarGaussianHighPrecision : public RooAbsPdf {
public:

  RooMultiVarGaussianHighPrecision() {} ;
  RooMultiVarGaussianHighPrecision(const char *name, const char *title, const RooArgList& xvec, const RooArgList& mu, const TMatrixDSym& covMatrix) ;
  void setAnaIntZ(Double_t z) { _z = z ; }

  RooMultiVarGaussianHighPrecision(const RooMultiVarGaussianHighPrecision& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooMultiVarGaussianHighPrecision(*this,newname); }
  inline virtual ~RooMultiVarGaussianHighPrecision() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ; 
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ; 

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const; 
  void initGenerator(Int_t code) ;
  void generateEvent(Int_t code); 

  const TMatrixDSym& covarianceMatrix() const { return _cov ; }
  
  class AnaIntData {
  public:
    TMatrixD    S22bar ;
    cpp_dec_float_50    S22det ;
    std::vector<int> pmap ;
    Int_t       nint ;
  } ;

  class GenData {
  public:
    TMatrixD    UT ;
    std::vector<int> omap ;
    std::vector<int> pmap ;
    TVectorD    mu1 ;
    TVectorD    mu2 ;
    TMatrixD    S12S22I ;
  } ;

  class BitBlock {
  public:
    BitBlock() : b0(0), b1(0), b2(0), b3(0) {} ;

    void setBit(Int_t ibit) ;      
    Bool_t getBit(Int_t ibit) ;
    Bool_t operator==(const BitBlock& other) ;

    Int_t b0 ;
    Int_t b1 ;
    Int_t b2 ;
    Int_t b3 ;
  } ;

  static void blockDecompose(const TMatrixD& input, const std::vector<int>& map1, const std::vector<int>& map2, TMatrixDSym& S11, TMatrixD& S12, TMatrixD& S21, TMatrixDSym& S22) ;

protected:
  
  void decodeCode(Int_t code, std::vector<int>& map1, std::vector<int>& map2) const;
  AnaIntData& anaIntData(Int_t code) const ;
  GenData& genData(Int_t code) const ;

  mutable std::map<int,AnaIntData> _anaIntCache ; //!
  mutable std::map<int,GenData> _genCache ; //!

  mutable std::vector<BitBlock> _aicMap ; //!

  RooListProxy _x ;
  RooListProxy _mu ;
  TMatrixDSym _cov ;
  TMatrixDSym _covI ;
  cpp_dec_float_50    _det ; 
  Double_t    _z ; 

  void syncMuVec() const ;
  mutable TVectorD _muVec ; //! Do not persist

  Double_t evaluate() const ;

private:

  ClassDef(RooMultiVarGaussianHighPrecision,1) // Multivariate Gaussian PDF with correlations
};

#endif
