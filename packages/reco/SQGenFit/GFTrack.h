#ifndef _GFTRACK_H
#define _GFTRACK_H

#include <vector>
#include <set>
#include <memory>

#include <GenFit/AbsTrackRep.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/Track.h>
#include <GenFit/Tools.h>

#include <TString.h>
#include <TVector3.h>

#include "SRecEvent.h"
#include "FastTracklet.h"
#include "GFMeasurement.h"

namespace SQGenFit
{
class GFTrack
{
public:
  typedef std::set<GFMeasurement*, GFMeasurementComp> HitList_t;

  GFTrack();
  GFTrack(SRecTrack& recTrack);
  GFTrack(Tracklet& tracklet);
  ~GFTrack();

  //! create a clone of the current track - either by track candidate or by simple copy
  GFTrack* clone(bool clean = false) const;

  //! convert the genfit track to SRecTrack to be saved in the DST tree
  SRecTrack getSRecTrack();

  void setVerbosity(unsigned int v);

  //! initialize the track by a tracklet, default reference point is in front of the station-1
  void setTracklet(const Tracklet& tracklet, double z_reference = 590., bool wildseedcov = false);
  void addMeasurements(HitList_t& measurements);

  //! id = -1 means insert hits after the last point, -2 means before the last point
  //! id = 0 means before the first point, 1 means before the second point
  void addMeasurement(GFMeasurement* measurement, const int id = -1);

  //! add hit to the track, default insert position is to the front, assuming we always add hits from downstream to upstream
  void addMeasurement(const SignedHit& hit, const int id = 0);

  // void addVertex(double zvtx);

  //! try adding a hit and calculate/return the incremental chi2
  double testOneHitKalman(const SignedHit& hit);

  //! remove all the fitted state, optionally update the seed
  void reset(bool updateSeed = true);

  //! get the projected x/y and w/dw on the detector plane
  void getProjection(int detID, double& x, double& y, double& w, double& dw);

  //! get the track fit quality parameters
  double getChi2();
  double getNDF();
  int getCharge();

  //! The extrapolation is implemented for line, plane and point, but the update/filter is only implemeted 
  //! for line and plane, user needs to be careful and pass correct measurement and cov to get sensible result
  double extrapolateToLine(TVector3& endPoint1, TVector3& endPoint2, const int startPtID = 0);
  double extrapolateToPlane(TVector3& pO, TVector3& pU, TVector3& pV, const int startPtID = 0);
  double extrapolateToPoint(TVector3& point, const int startPtID = 0);
  void getExtrapPosMomCov(TVector3* pos, TVector3* mom, TMatrixDSym* cov);

  //! perform one step of filter for the _propState, returns the chi2/ndf
  double updatePropState(const TVectorD& meas, const TMatrixDSym& V);
  
  //! propagate the track to a plane at the expected z_vertex location, the x and y are obtained from recoConsts
  //! set bias to false to do the extrapolation only
  double swimToVertex(double z, TVector3* pos = nullptr, TVector3* mom = nullptr, TMatrixDSym* cov = nullptr, bool biased = true);

  //! load the extrapolation cache to calculate the z location of the point of the closest approach, needs to be called after extrapolateToXXX
  double getPOCA(SRecTrack* strack = nullptr); 

  //! load the extrapolation cache to calculate the total distance of the track traveled inside the beam dump
  double getDumpPathLen();

  //! throw exception if the internal object onwership is messed up, implemented in GenFit
  void checkConsistency();

  //! this function will restore the information from a GenFit track, including: trkrep, measurements, pdg code, etc.
  void restoreFromGenFitTrack();

  //! update all the hits uinfo (residual, chi2, etc.) after fit
  void postFitUpdate(bool updateMeasurements = true);
  
  //! debug output
  void print(unsigned int debugLvl = 0);

  genfit::Track* getGenFitTrack() { return _track; }
  genfit::AbsTrackRep* getGenFitTrkRep() { return _trkrep; }

private:
  //! Auxilary function to initialize the extrapolation/projection, should not be seen by external user
  bool setInitialStateForExtrap(const int startPtID = 0);

  //! interface with genfit operations
  genfit::AbsTrackRep* _trkrep; 
  genfit::Track* _track;
  HitList_t _measurements;   // --- TODO: should I make this a list instead?

  //! container of the MSOP obtained from the SRecTrack
  std::vector<genfit::MeasuredStateOnPlane> _fitstates;

  //! Used for the propagation and hypothesis test
  std::unique_ptr<genfit::MeasuredStateOnPlane> _propState;
  std::unique_ptr<genfit::AbsMeasurement> _virtMeas;

  //! Store the original track candidate information from track finding
  Tracklet* _trkcand;
  int _pdg;

};
}

#endif