#include "GFField.h"

#include <iostream>

#include <CLHEP/Units/SystemOfUnits.h>

namespace SQGenFit
{
GFField* GFField::_p_field = nullptr;
GFField* GFField::instance(const PHField* field)
{
  if(field != nullptr)
  {
    _p_field = new GFField(field);
  }

  return _p_field;
}

GFField::GFField(const PHField* field): 
  _field(field),
  _scale(1.),
  _zoffset(0.),
  _disable(false)
{
  if(field == nullptr)
  {
    std::cerr << "PHField not intialted for the GFField !!" << std::endl; 
  }
}

TVector3 GFField::get(const TVector3& pos) const
{
  double x = pos.X();
  double y = pos.Y();
  double z = pos.Z();

  double Bx, By, Bz;
  get(x, y, z, Bx, By, Bz);

  return TVector3(Bx, By, Bz);
}

void GFField::get(const double& x, const double& y, const double& z, double& Bx, double& By, double& Bz) const
{
  if(_disable)
  {
    Bx = 0.;
    By = 0.;
    Bz = 0.;
    return;
  }

  double Point[] = {x*CLHEP::cm, y*CLHEP::cm, z*CLHEP::cm, 0.};
  if(z < 650.) 
  {
    Point[2] -= (_zoffset*CLHEP::cm); //we assume positive zoffset means moving magnet downstream, which is equivalent to moving track upstream
    //std::cout << z << "  " << Point[2] << "  " << _zoffset << std::endl;
  }

  double Bfield[6];
  for(int i = 0; i < 6; ++i)
  {
    Bfield[i] = 0.;
  }

  _field->GetFieldValue(Point, Bfield);
  Bx = _scale*Bfield[0]/CLHEP::kilogauss;
  By = _scale*Bfield[1]/CLHEP::kilogauss;
  Bz = _scale*Bfield[2]/CLHEP::kilogauss;
}

}