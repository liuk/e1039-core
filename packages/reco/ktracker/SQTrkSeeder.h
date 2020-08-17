#ifndef _SQTRKSEEDER_H
#define _SQTRKSEEDER_H

#include <set>
#include <list>
#include <unordered_map>
#include <string>

#include <Math/Factory.h>
#include <Math/Minimizer.h>

#include <GlobalConsts.h>
#include <phool/PHTimer.h>
#include <geom_svc/GeomSvc.h>

#include "SRawEvent.h"
#include "SRawEventSvc.h"
#include "FastTracklet.h"

class TrackingStationInfo
{
public:
  enum StID_t {DC0, DC1, DC2, DC3p, DC3m};
  StID_t stationID;
  std::string stationName();

  int xPairID;
  double xZ;

  int uPairID;
  double uZ;
  double uCostheta;
  double uSintheta;
  double uCellWidth;
  double uWinSize;

  int vPairID;
  double vZ;
  double vWinSize;
  double vWinSizeX;
  double vWinSizeY;
};

class SQTrkSeeder
{
public:
  SQTrkSeeder();
  ~SQTrkSeeder();

  void registerSeedingStation(TrackingStationInfo::StID_t stID);

  void processEvent();

  void buildSeedInSt(TrackingStationInfo::StID_t stID);

  int fitSeed(Tracklet& seed);

  bool acceptSeed(Tracklet& seed);

  std::list<Tracklet>& getSeeds() { return seeds; }

  void printTimer();

private:
  //! stationIDs used for the seeder
  std::set<TrackingStationInfo::StID_t> seedStIDs;

  //! RawHit query tool
  SRawEventSvc* p_rawEvtSvc;

  GeomSvc* p_geomSvc;

  //! container of cached tracking station info
  std::unordered_map<int, TrackingStationInfo> stationInfo;
  std::unordered_map<int, PHTimer> stationTimer;

  //! list of the seeds created
  std::list<Tracklet> seeds;

  //! Least chi square fitter and functor
  ROOT::Math::Minimizer* minimizer;

  //! cache for chamber z position
  double z_plane[nChamberPlanes+1];
};

#endif
