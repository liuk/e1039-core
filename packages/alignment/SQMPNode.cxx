#include "SQMPNode.h"

#include <ktracker/FastTracklet.h>
#include <geom_svc/GeomSvc.h>
#include <iostream>

ClassImp(SQMPNode)

namespace
{
  //static flag to indicate the initialized has been done
  static bool inited = false;

  //static flag of kmag on/off
	static bool KMAG_ON;

  //initialize global variables
  void initGlobalVariables()
  {
    if(!inited) 
    {
      inited = true;

      recoConsts* rc = recoConsts::instance();
      KMAG_ON = rc->get_BoolFlag("KMAG_ON");
    }
  }
};

SQMPNode::SQMPNode():
  detectorID(-1),
  flag(0)
{
  initGlobalVariables();
}

SQMPNode::SQMPNode(int detector_index):
  detectorID(detector_index),
  flag(0)
{
  initGlobalVariables();
}

SQMPNode::SQMPNode(SignedHit& hit_signed, Tracklet& trk):
  detectorID(hit_signed.hit.detectorID),
  elementID(hit_signed.hit.elementID),
  sign(hit_signed.sign),
  charge(trk.getCharge()),
  tdctime(hit_signed.hit.tdcTime),
  drift(hit_signed.hit.driftDistance),
  pos(hit_signed.hit.pos),
  flag(1),
  tx(trk.tx),
  ty(trk.ty),
  x0(trk.x0),
  y0(trk.y0)
{
  initGlobalVariables();

  if(detectorID < 1 || detectorID > nChamberPlanes)
  {
    flag = 0;
    return;
  }

  //Measurements
  meas = trk.residual[detectorID-1];
  // if(fabs(meas) > 1.5 )
  // {
  //     flag = false;
  //     return;
  // }

  //Local parameters
  GeomSvc* p_geomSvc = GeomSvc::instance();
  z = p_geomSvc->getPlanePosition(detectorID);
  if(KMAG_ON && detectorID <= 12) trk.getXZInfoInSt1(tx, x0);

  //Plane tilt angles
  double cosphi = p_geomSvc->getCostheta(detectorID);
  double sinphi = p_geomSvc->getSintheta(detectorID);
  sigma = p_geomSvc->getPlaneResolution(detectorID);

  //Fill direivatives
  setDerivatives(z, cosphi, sinphi, tx, ty, x0, y0);
}

void SQMPNode::setDerivatives(double z, double cosphi, double sinphi, double tx, double ty, double x0, double y0)
{
  //Global derivatives
  dwdz = (tx*cosphi + ty*sinphi);
  dwdphi = (-x0*sinphi - tx*z*sinphi + y0*cosphi + ty*z*cosphi);
  dwdw = -1.;

  //Local derivatives
  dwdx = cosphi;
  dwdy = sinphi;
  dwdtx = z*cosphi;
  dwdty = z*sinphi;
}

void SQMPNode::identify(std::ostream& os) const
{
  using namespace std;

  os << "========= Alignment node of detector: " << detectorID << " =========" << endl;
  if(!isValid())
  {
    os << "No measurement on this node." << endl;
    return;
  }

  os << "Measurements: " << meas << " +/- " << sigma << endl;

  os << "Local track parameters: " << endl;
  os << "z      x0     y0     dx/dz     dy/dz" << endl;
  os << z << "    " << x0 << "    " << y0 << "    " << tx << "    " << ty << endl;

  os << "Global derivative:" << endl;
  os << "dw/dz       dw/dphi       dw/dw" << endl;
  os << dwdz << "   " << dwdphi << "    " << dwdw << endl;

  os << "Local derivative:" << endl;
  os << "dw/dx     dw/dy      dw/dtx      dw/dty" << endl;
  os << dwdx << "    " << dwdy << "    " << dwdtx << "    " << dwdty << endl;
}