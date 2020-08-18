#include "SQTrkSeeder.h"

#include <Math/Functor.h>
#include <bitset>

namespace 
{
  //static flag to indicate the initialized has been done
  static bool inited = false;

  //Track quality cuts
  static double TX_MAX;
  static double TY_MAX;
  static double X0_MAX;
  static double Y0_MAX;

  //initialize global variables
  void initGlobalVariables()
  {
    if(!inited) 
    {
      inited = true;

      recoConsts* rc = recoConsts::instance();

      TX_MAX = rc->get_DoubleFlag("TX_MAX");
      TY_MAX = rc->get_DoubleFlag("TY_MAX");
      X0_MAX = rc->get_DoubleFlag("X0_MAX");
      Y0_MAX = rc->get_DoubleFlag("Y0_MAX");
    }
  }
}

std::string TrackingStationInfo::stationName()
{
  switch(stationID)
  {
    case DC0: return "D0";
    case DC1: return "D1";
    case DC2: return "D2";
    case DC3p: return "D3p";
    case DC3m: return "D3m";
    default: return "DC";
  }
}

SQTrkSeeder::SQTrkSeeder(): p_geomSvc(GeomSvc::instance()), p_rawEvtSvc(SRawEventSvc::instance())
{
  initGlobalVariables();

  //Minimize ROOT output
  extern Int_t gErrorIgnoreLevel;
  gErrorIgnoreLevel = 9999;

  minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  minimizer->SetMaxFunctionCalls(10000);
  minimizer->SetMaxIterations(20);
  minimizer->SetTolerance(1E-2);
  minimizer->SetPrintLevel(0);

  for(int i = 1; i <= nChamberPlanes; ++i) z_plane[i] = p_geomSvc->getPlanePosition(i);
}

SQTrkSeeder::~SQTrkSeeder()
{
  delete minimizer;
}

void SQTrkSeeder::registerSeedingStation(TrackingStationInfo::StID_t stID)
{
  TrackingStationInfo stInfo;
  stInfo.stationID = stID;

  stInfo.xPairID = (p_geomSvc->getDetectorID(stInfo.stationName() + "X") + 1)/2;
  stInfo.xZ = 0.5*(z_plane[2*stInfo.xPairID-1] + z_plane[2*stInfo.xPairID]);

  stInfo.uPairID = (p_geomSvc->getDetectorID(stInfo.stationName() + "U") + 1)/2;
  stInfo.uZ = 0.5*(z_plane[2*stInfo.uPairID-1] + z_plane[2*stInfo.uPairID]);
  stInfo.uCostheta = p_geomSvc->getCostheta(2*stInfo.uPairID-1);
  stInfo.uSintheta = p_geomSvc->getSintheta(2*stInfo.uPairID-1);
  stInfo.uCellWidth = p_geomSvc->getCellWidth(2*stInfo.uPairID-1);
  double uScaleY = p_geomSvc->getPlaneScaleY(2*stInfo.uPairID-1);
  stInfo.uWinSize = fabs(0.5*uScaleY*stInfo.uSintheta) + TX_MAX*fabs((stInfo.uZ- stInfo.xZ)*stInfo.uCostheta) + TY_MAX*fabs((stInfo.uZ - stInfo.xZ)*stInfo.uSintheta) + 2.*stInfo.uSintheta + 10.;  //TODO: 10 is something we need to optimize

  stInfo.vPairID = (p_geomSvc->getDetectorID(stInfo.stationName() + "V") + 1)/2;
  stInfo.vZ = 0.5*(z_plane[2*stInfo.vPairID-1] + z_plane[2*stInfo.vPairID]);
  stInfo.vWinSize = 2.*stInfo.uCellWidth*stInfo.uCostheta;
  stInfo.vWinSizeX = stInfo.uCostheta*TX_MAX;
  stInfo.vWinSizeY = stInfo.uSintheta*TY_MAX;

  seedStIDs.insert(stID);
  stationInfo[stID] = stInfo;
  stationTimer[stID] = PHTimer(stInfo.stationName());
  stationTimer[stID].reset();

  LogDebug(stID << "  " << stInfo.stationName() << "  " << stInfo.xPairID << "  " << stInfo.uPairID << "  " << stInfo.vPairID);
}

void SQTrkSeeder::processEvent()
{
  seeds.clear();
  for(auto it = seedStIDs.begin(); it != seedStIDs.end(); ++it)
  {
    stationTimer[*it].restart();
    buildSeedInSt(*it);
    stationTimer[*it].stop();
  }

  seeds.sort();
  if(seeds.size() > 200) seeds.resize(200);
}

bool SQTrkSeeder::acceptSeed(Tracklet& seed)
{
  if(seed.isValid() == 0) return false;

  // if(!(p_rawEvtSvc->hodoMasked(39, seed.getExpPositionW(39), 1) && p_rawEvtSvc->hodoMasked(37, seed.getExpPositionW(37), 2)) && 
  //    !(p_rawEvtSvc->hodoMasked(40, seed.getExpPositionW(40), 1) && p_rawEvtSvc->hodoMasked(38, seed.getExpPositionW(38), 2)))
  //    return false;
  return true;
}


void SQTrkSeeder::buildSeedInSt(TrackingStationInfo::StID_t stID)
{  
  auto pairs_X = p_rawEvtSvc->getHitPairsInDetectorPair(stationInfo[stID].xPairID);
  if(pairs_X.empty()) return;
  auto pairs_U = p_rawEvtSvc->getHitPairsInDetectorPair(stationInfo[stID].uPairID);
  if(pairs_U.empty()) return;
  auto pairs_V = p_rawEvtSvc->getHitPairsInDetectorPair(stationInfo[stID].vPairID);
  if(pairs_V.empty()) return;
  LogDebug(pairs_X.size() << " x pairs, " << pairs_U.size() << " u pairs, " << pairs_V.size() << " v pairs.");

  //X-U combination first, then add V pairs
  std::vector<Hit>& hits = p_rawEvtSvc->hits();
  for(auto xiter = pairs_X.begin(); xiter != pairs_X.end(); ++xiter)
  {
    //U projections from X plane
    double x_pos = xiter->second >= 0 ? 0.5*(hits[xiter->first].pos + hits[xiter->second].pos) : hits[xiter->first].pos;
    double u_min = x_pos*stationInfo[stID].uCostheta - stationInfo[stID].uWinSize;
    double u_max = u_min + 2.*stationInfo[stID].uWinSize;


    LogDebug("Trying X hits " << xiter->first << "  " << xiter->second << "  " << hits[xiter->first].elementID << " at " << x_pos);
    LogDebug("U plane window:" << u_min << "  " << u_max);
    for(auto uiter = pairs_U.begin(); uiter != pairs_U.end(); ++uiter)
    {
      double u_pos = uiter->second >= 0 ? 0.5*(hits[uiter->first].pos + hits[uiter->second].pos) : hits[uiter->first].pos;
      LogDebug("Trying U hits " << uiter->first << "  " << uiter->second << "  " << hits[uiter->first].elementID << " at " << u_pos);
      
      if(u_pos < u_min || u_pos > u_max) continue;

      //V projections from X and U plane
      double z_x = xiter->second >= 0 ? stationInfo[stID].xZ : z_plane[hits[xiter->first].detectorID];
      double z_u = uiter->second >= 0 ? stationInfo[stID].uZ : z_plane[hits[uiter->first].detectorID];
      double z_v = stationInfo[stID].vZ;
      double v_win1 = stationInfo[stID].vWinSize;
      double v_win2 = fabs((z_u + z_v - 2.*z_x)*stationInfo[stID].vWinSizeX);
      double v_win3 = fabs((z_v - z_u)*stationInfo[stID].vWinSizeY);
      double v_win = v_win1 + v_win2 + v_win3 + 2.*stationInfo[stID].uCellWidth;
      double v_min = 2.*x_pos*stationInfo[stID].uCostheta - u_pos - v_win;
      double v_max = v_min + 2.*v_win;
      LogDebug("V plane window:" << v_min << "  " << v_max);

      for(auto viter = pairs_V.begin(); viter != pairs_V.end(); ++viter)
      {
        double v_pos = viter->second >= 0 ? 0.5*(hits[viter->first].pos + hits[viter->second].pos) : hits[viter->first].pos;
        LogDebug("Trying V hits " << viter->first << "  " << viter->second << "  " << hits[viter->first].elementID << " at " << v_pos);

        if(v_pos < v_min || v_pos > v_max) continue;
        if(xiter->second < 0 && uiter->second < 0 && viter->second < 0) continue;

        //Now create a seed with large momentum to remove the multiple scattering effect
        int LR1 = 0;
        int LR2 = 0;
        Tracklet seed;
        seed.stationID = stID;
        seed.invP = 0.01;

        if(xiter->first >= 0)
        {
          seed.hits.push_back(SignedHit(hits[xiter->first], LR1));
          ++seed.nXHits;
        }
        if(xiter->second >= 0)
        {
          seed.hits.push_back(SignedHit(hits[xiter->second], LR2));
          ++seed.nXHits;
        }

        if(uiter->first >= 0)
        {
          seed.hits.push_back(SignedHit(hits[uiter->first], LR1));
          ++seed.nUHits;
        }
        if(uiter->second >= 0)
        {
          seed.hits.push_back(SignedHit(hits[uiter->second], LR2));
          ++seed.nUHits;
        }

        if(viter->first >= 0)
        {
          seed.hits.push_back(SignedHit(hits[viter->first], LR1));
          ++seed.nVHits;
        }
        if(viter->second >= 0)
        {
          seed.hits.push_back(SignedHit(hits[viter->second], LR2));
          ++seed.nVHits;
        }

        fitSeed(seed);
#ifdef _DEBUG_ON
        seed.print();
#endif
        if(acceptSeed(seed))
        {
          LogDebug("Accepted!");
          seeds.push_back(seed);
        }
        else
        {
          LogDebug("Rejected!");
        }
      }
    }
  }

  //TODO: maybe set a cap to each seeding station separately? create local list and move it using std::move or sth similar
}

int SQTrkSeeder::fitSeed(Tracklet& seed)
{
  ROOT::Math::Functor fcn(&seed, &Tracklet::Eval4, 4);
  minimizer->SetFunction(fcn);

  minimizer->SetLimitedVariable(0, "tx", seed.tx, 0.001, -TX_MAX, TX_MAX);
  minimizer->SetLimitedVariable(1, "ty", seed.ty, 0.001, -TY_MAX, TY_MAX);
  minimizer->SetLimitedVariable(2, "x0", seed.x0, 0.1,   -X0_MAX, X0_MAX);
  minimizer->SetLimitedVariable(3, "y0", seed.y0, 0.1,   -Y0_MAX, Y0_MAX);
  minimizer->Minimize();

  seed.tx = minimizer->X()[0];
  seed.ty = minimizer->X()[1];
  seed.x0 = minimizer->X()[2];
  seed.y0 = minimizer->X()[3];

  seed.err_tx = minimizer->Errors()[0];
  seed.err_ty = minimizer->Errors()[1];
  seed.err_x0 = minimizer->Errors()[2];
  seed.err_y0 = minimizer->Errors()[3];
  seed.chisq  = minimizer->MinValue();

  int status = minimizer->Status();
  return status;
}

void SQTrkSeeder::printTimer()
{
  std::cout << "============================== SQTrkSeeder timer info ====================================" << std::endl;
  for(auto it = stationTimer.begin(); it != stationTimer.end(); ++it)
  {
    it->second.print_stat();
  }
}
