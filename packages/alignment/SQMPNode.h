#ifndef _SQMPNODE_H
#define _SQMPNODE_H

#include <iostream>
#include <phool/PHObject.h>

class SignedHit;
class Tracklet;

class SQMPNode: public PHObject
{
public:
  SQMPNode();
  SQMPNode(int detector_index);
  SQMPNode(SignedHit& hit_input, Tracklet& trk);

  void      identify(std::ostream& os = std::cout) const;
  void      Reset() { flag = -1; detectorID = -1; }
  int       isValid() const { return flag; }
  SQMPNode* Clone()   const { return (new SQMPNode(*this)); }

  //Calculate derivatives
  void setDerivatives(double z, double cosphi, double sinphi, double tx, double ty, double x0, double y0);

  //Flag indicating whether this node is valid or not
  int flag;

  //detector ID
  int detectorID;
  int elementID;

  //Left/right
  int sign;

  //Measurement (actually residual due to the definition of millepede)
  double meas;

  //other auxilary info
  int charge;
  double tdctime;
  double pos;
  double drift;

  //Resolution of the measurement
  double sigma;

  //Derivatives w.r.t global parameters
  double dwdz;
  double dwdphi;
  double dwdw;

  //Derivarives w.r.t local parameters
  double dwdx;
  double dwdy;
  double dwdtx;
  double dwdty;

  //Local parameters
  double x0;            // x position
  double y0;            // y position
  double tx;            // dxdz, i.e. px/pz
  double ty;            // dydz, i.e. py/pz

  //z position of the node
  double z;

  //Overiden comparison operator
  bool operator<(const SQMPNode& elem) const { return detectorID < elem.detectorID; }

  ClassDef(SQMPNode, 1);
};


#endif
