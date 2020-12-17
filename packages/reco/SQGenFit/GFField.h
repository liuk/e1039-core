#ifndef _GFFIELD_H
#define _GFFIELD_H

#include <GenFit/AbsBField.h>
#include <TVector3.h>

#include <phfield/PHField.h>

namespace SQGenFit
{
class GFField: public genfit::AbsBField
{
public:
  static GFField* instance(const PHField* field = nullptr);

  GFField(const PHField* field);
  virtual ~GFField() {}

  TVector3 get(const TVector3& pos) const;
  void get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const;

  void setScale(double scale) { _scale = scale; }
  void setZOffset(double offset) { _zoffset = offset; }
  void disable() { _disable = true; }

private:
  const PHField* _field;
  double _scale;
  double _zoffset;
  bool   _disable;

  static GFField* _p_field;
};
}

#endif