#ifndef _Angles_H
#define _Angles_H


#include "TLorentzVector.h"

// c++ libraries
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

class Angles {
 public:
  Angles();
  virtual ~Angles();
  float getCosThetaStar_CS(TLorentzVector h1, TLorentzVector h2, float ebeam = 3500);
  void SetZGAngles(const TLorentzVector& l1, const TLorentzVector& l2, const TLorentzVector& g);
  void GetZGAngles(double& cos1, double& cos2, double& phi, double& cos3);
  void GetZGAngles(const TLorentzVector& l1, const TLorentzVector& l2, const TLorentzVector& g, double& cos1, double& cos2, double& phi, double& cos3);
  double GetCos1();
  double GetCos2();
  double GetCosTheta();
  double GetPhi();
    
 private:
  double _costheta_lplus, _costheta_lminus, _phi, _cosTheta;
};


#endif
