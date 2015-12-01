#include "../interface/Angles.h"

// Constructor, set default values
Angles::Angles():
  _costheta_lplus(-99),
  _costheta_lminus(-88),
  _phi(-77),
  _cosTheta(-66)
{}
Angles::~Angles(){}


float Angles::getCosThetaStar_CS(TLorentzVector h1, TLorentzVector h2, float ebeam) {
  // cos theta star angle in the Collins Soper frame
  // Copied directly from here: https://github.com/ResonantHbbHgg/Selection/blob/master/selection.h#L3367-L3385
  TLorentzVector p1, p2;
  p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
  p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

  TLorentzVector hh;
  hh = h1 + h2;

  TVector3 boost = - hh.BoostVector();
  p1.Boost(boost);
  p2.Boost(boost);
  h1.Boost(boost);

  TVector3 CSaxis = p1.Vect().Unit() - p2.Vect().Unit();
  CSaxis.Unit();

  return cos(   CSaxis.Angle( h1.Vect().Unit() )    );
}


void Angles::SetZGAngles(const TLorentzVector& lplus, const TLorentzVector& lminus, const TLorentzVector& gamma)
{
  //Set the angles (they are private variables of the class)

  TLorentzVector z = lplus + lminus;
  TLorentzVector h = z + gamma;

  //----------------------//
  //Boosting into CM frame//
  //----------------------//
  TLorentzVector beamAxis(0,0,1,1);
  TVector3 b1  = -1*h.BoostVector();  //this is the boost vector
  TLorentzVector h_CM(h);
  h_CM.Boost(b1);
  TLorentzVector z_CM(z);
  z_CM.Boost(b1);

  // ANGLE:
  _cosTheta = z_CM.Vect().Unit().Dot(h.Vect().Unit());

  TLorentzVector beamAxis_CM(beamAxis);
  beamAxis_CM.Boost(b1);


  //------------------------//
  //- Rotating the CM frame://
  //------------------------//
  TVector3 axis_z_CM = z_CM.Vect().Unit();
  TVector3 axis_y_CM = beamAxis_CM.Vect().Cross(z_CM.Vect()).Unit();
  TVector3 axis_x_CM = axis_y_CM.Cross(axis_z_CM).Unit();
  TRotation rotation;
  rotation = rotation.RotateAxes(axis_x_CM, axis_y_CM, axis_z_CM).Inverse();

  TLorentzVector lplus_CM(lplus), lminus_CM(lminus);
  lplus_CM.Boost(b1);
  lplus_CM.Transform(rotation);
  lminus_CM.Boost(b1);
  lminus_CM.Transform(rotation);


  //cout<<"\t\t z_CM before rotation  x = "<<z_CM.X()<<"  y="<<z_CM.Y()<<"  z="<<z_CM.Z()<<"  t="<<z_CM.T()<<endl;
  z_CM.Transform(rotation);
  //cout<<"\t\t z_CM after rotation  x = "<<z_CM.X()<<"  y="<<z_CM.Y()<<"  z="<<z_CM.Z()<<"  t="<<z_CM.T()<<endl;

  //This is checked are is gives the same results:
  //TLorentzVector hRot(h);
  //hRot.Transform(rotation);
  //Float_t c3 = hRot.CosTheta();
  //Float_t c2 = cos(z_CM.Angle(hRot.Vect()));

  //cout<<"higgs not boostedd  x = "<<h.X()<<"  y="<<h.Y()<<"  z="<<h.Z()<<"  t="<<h.T()<<endl;
  //cout<<"higgs boostedd   x = "<<h_CM.X()<<"  y="<<h_CM.Y()<<"  z="<<h_CM.Z()<<"  t="<<h_CM.T()<<endl;
  //cout<<"\t\t z boostedd  x = "<<z_CM.X()<<"  y="<<z_CM.Y()<<"  z="<<z_CM.Z()<<"  t="<<z_CM.T()<<endl;

  //cout<<"Debug \n"<<"\t\t c1 = "<<_cosTheta<<"\n \t\t c2 = "<<c2<<"\n \t\t c3 = "<<c3<<endl;

  if (fabs(_cosTheta)>1 || std::isnan(_cosTheta)){
    //cout<<"Cosine is greater than 1 or nan,  wha?"<<endl;
    exit(0);
  //throw cms::Exception("This is exceptional! ")<<"Cosine is greater than 1 or nan,  wha?"<<endl;
  }

  //Boosting to Z frame from the CM frame!
  TVector3 b2 = -1*z_CM.BoostVector();

  TLorentzVector lplus_inZFrame(lplus_CM);
  TLorentzVector lminus_inZFrame(lminus_CM);
  lminus_inZFrame.Boost(b2);
  lplus_inZFrame.Boost(b2);

  //Angle between the lepton(+) and a Z direction in the Z rest frame
  // ANGLES:
  _costheta_lplus  = lplus_inZFrame.CosTheta();
  _costheta_lminus = lminus_inZFrame.CosTheta();
  //_phi             = lplus_inZFrame.Phi();
  _phi             = lminus_inZFrame.Phi();
  //Filled!
  //cout<<"c1 ="<<_costheta_lplus<<"  c2 ="<<_costheta_lminus<<"  phi ="<<_phi<<"c1 ="<<_cosTheta<<endl;

  //cout<<"lplus in Z frame: px = "<<lplus_inZFrame.Px()<<"  py = "<<lplus_inZFrame.Py()<<endl;
  //cout<<"lminus in Z frame: px = "<<lminus_inZFrame.Px()<<"  py = "<<lminus_inZFrame.Py()<<endl;



  //ATT:
  //There are some problems with other definitions (c2, and c3). They are not always the same
  //It has to be investigated later.

  //double c2 = (lplus_CM.E() - lminus_CM.E())/ (lplus_CM.Vect() + lminus_CM.Vect()).Mag(); //from a formula
  //if (abs(costheta1)>1 || isnan(costheta1)
  // throw cms::Exception("This is exceptional! ")<<"Cosine is greater than 1,  wha?"<<endl;

  //Boosting directly to Z
  // TVector3 b3  = -1*z.BoostVector();
  //TLorentzVector lplus_ZFrame(lplus);
  //TLorentzVector lminus_ZFrame(lminus);
  //lplus_ZFrame.Boost(b3);
  //lminus_ZFrame.Boost(b3);

  //cout<<"\t lplus in Z frame: px = "<<lplus_ZFrame.Px()<<"  py = "<<lplus_ZFrame.Py()<<endl;
  //cout<<"\t lminus in Z frame: px = "<<lminus_ZFrame.Px()<<"  py = "<<lminus_ZFrame.Py()<<endl;

  //c3 = lplus_ZFrame.Vect().Unit().Dot(b3.Unit());
  //double  c3 = cos(TMath::Pi() - lplus_ZFrame.Angle(b3));
  //cout<<"Debug \n"<<"\t\t cos_lp = "<<_costheta_lplus<<"\t cos_lm = "<<_costheta_lminus<<"\n \t\t c2 = "<<c2<<"\n \t\t c3 = "<<c3<<"\n \t\t c4 = "<<endl;
  //costheta1 = c3;



}
void Angles::GetZGAngles(double& costheta_lplus, double& costheta_lminus, double& phi, double& cosTheta)
{
  //Just a simple getter to acces private angles
  costheta_lplus  = _costheta_lplus;
  costheta_lminus = _costheta_lminus;
  phi = _phi;
  cosTheta = _cosTheta;
}


void Angles::GetZGAngles(const TLorentzVector& lplus, const TLorentzVector& lminus, const TLorentzVector& g,
			 double& costheta_lplus, double& costheta_lminus, double& phi, double& cosTheta)
{
  SetZGAngles(lminus, lplus, g);

  costheta_lplus = _costheta_lplus;
  costheta_lminus = _costheta_lminus;
  phi = _phi;
  cosTheta = _cosTheta;
}

double Angles::GetCos1(){return _costheta_lplus;}
double Angles::GetCos2(){return _costheta_lminus;}
double Angles::GetCosTheta(){return _cosTheta;}
double Angles::GetPhi(){return _phi;}
