#include "SQTrkProp.h"

#include <phool/recoConsts.h>
#include <phfield/PHFieldConfig_v3.h>
#include <phfield/PHFieldUtility.h>
#include <phgeom/PHGeomUtility.h>

#include <GlobalConsts.h>

SQTrkProp::SQTrkProp(const PHField* field): 
  p_geomSvc(GeomSvc::instance()), 
  p_rawEvtSvc(SRawEventSvc::instance()),
  timer(PHTimer("TrkProp")),
  fitter(nullptr)
{
  sequence = std::vector<int>{9, 8, 7, 3, 2, 1};

  //init fitter
  gfield = new SQGenFit::GFField(field);
  fitter = new SQGenFit::GFFitter();
  fitter->init(gfield, "KalmanFitterRefTrack");

  timer.reset();
}

SQTrkProp::~SQTrkProp()
{
  if(fitter != nullptr) delete fitter;
}

void SQTrkProp::processSeed(Tracklet& seed)
{
  timer.restart();
#ifdef _DEBUG_ON
  seed.identify();
#endif
  SQGenFit::GFTrackPtr seedtrk(new SQGenFit::GFTrack());
  seedtrk->setTracklet(seed, 1900.);

  int status = fitter->processTrack(*seedtrk);
  if(status != 0) 
  {
    LogDebug("Seed track failed fitting: " << status);
    return;
  }

  tracks.clear();
  tracks.push_back(seedtrk);
  for(auto pairID = sequence.begin(); pairID != sequence.end(); ++pairID)
  {
    LogDebug("Processing detector pair: " << *pairID);
    LogDebug(tracks.size() << " base tracks to process ");

    std::list<SQGenFit::GFTrackPtr> newTracks;
    for(auto track = tracks.begin(); track != tracks.end(); ++track)
    {
#ifdef _DEBUG_ON
      LogDebug("Process this track: ");
      (*track)->print();
#endif

      propagateTo(track->get(), *pairID, newTracks);
      LogDebug("Now we have " << newTracks.size() << " new tracks");
    }

    trimTracklist(newTracks);
    if(newTracks.empty()) break;

    tracks.assign(newTracks.begin(), newTracks.end());
  }
  tracks.sort(SQGenFit::GFTrackPtrComp());
  timer.stop();
}

void SQTrkProp::propagateTo(SQGenFit::GFTrack* btrk, int detPairID, std::list<SQGenFit::GFTrackPtr>& tracklist)
{
  double x, y, w, dw;
  btrk->getProjection(2*detPairID, x, y, w, dw);
  LogDebug("Projected position and windows on detPair " << detPairID << ": " << x << " " << y << "  " << w << " +/- " << dw);
  
  std::list<SRawEventSvc::HitPair_t> hitlist = p_rawEvtSvc->getHitPairsInDetectorPair(detPairID, w, 3.*dw);
  if(hitlist.empty()) return;

  LogDebug(hitlist.size() << " hit pairs to add");
  for(auto it = hitlist.begin(); it != hitlist.end(); ++it)
  {
    LogDebug("Try adding hit " << it->first << " and " << it->second << " at p = " << p_rawEvtSvc->hit(it->first).pos);
#ifdef _DEBUG_ON
    if(it->first  >= 0) p_rawEvtSvc->hit(it->first ).identify();
    if(it->second >= 0) p_rawEvtSvc->hit(it->second).identify();
#endif
    /* 
    Two tings we could do to potential speed up:
    1. use GFTrack::testOneHitKalman() to get a quick one step fit and cut on chi2
    2. use GFTrack::testOneHitKalman() to get the chi2 for all possible candidates, sort by chi2/nhits, select only the top n pairs 
       for a full fit
    */
    SQGenFit::GFTrackPtr newTrack(btrk->clone());
    if(it->first  >= 0) newTrack->addMeasurement(SignedHit(p_rawEvtSvc->hit(it->first),  0));
    if(it->second >= 0) newTrack->addMeasurement(SignedHit(p_rawEvtSvc->hit(it->second), 0));
    
    int status = fitter->processTrack(*newTrack);
    if(status == 0) 
    {
      tracklist.push_back(newTrack);
      LogDebug("Fit successful for this pair");
    }
    else
      LogDebug("Fit failed for this pair: " << status);
  }
}

void SQTrkProp::trimTracklist(std::list<SQGenFit::GFTrackPtr>& tracklist)
{
  for(auto track = tracklist.begin(); track != tracklist.end(); )
  {
    if(!acceptTrack(track->get()))
    {
      track = tracklist.erase(track);
    }
    else
    {
      ++track;
    }
  }
}

bool SQTrkProp::acceptTrack(SQGenFit::GFTrack* track)
{
  if(track->getQuality() > 10.) return false;

  TVector3 pos, mom;
  track->getFittedPosMom(pos, mom);
  if(mom.Mag() < 10.) return false;

  return true;
}

