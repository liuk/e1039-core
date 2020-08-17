#ifndef _SRAWEVENTSVC_H
#define _SRAWEVENTSVC_H

#include "SRawEvent.h"

#include <GlobalConsts.h>
#include <geom_svc/GeomSvc.h>

class SRawEventSvc
{
public:
  typedef std::pair<int, int> HitPair_t;

  static SRawEventSvc* instance();

  SRawEventSvc();

  bool setRawEvent(SRawEvent* r, bool buildAuxInfo = true);
  SRawEvent* getRawEvent() { return rawEvent; }

  //Hit list querying tools
  std::list<int> getHitsIndexInDetector(short detectorID, double p_exp = 0., double win = 99999.);
  std::list<int> getHitsIndexInDetectorPair(short pairID);  //pairID is the id of prime/unprime pairs
  std::list<int> getHitsIndexInDetectors(std::vector<int>& detectorIDs);
  std::list<HitPair_t> getHitPairsInDetectorPair(short pairID, double p_exp = 0., double win = 99999.);

  bool acceptEvent(SRawEvent* r);

  //list of hits from the SRawEvent with public access
  std::vector<Hit>& hits() { return *hitsPtr; }

  //access to single hit reference
  Hit& hit(int idx) { return hitsPtr->at(idx); }

  //hodoscope masking
  bool hodoMasked(int detID, int eleID, int win);
  bool hodoMasked(int detID, double pos, int win);

  //get the hodo patter
  unsigned int hodoPattern(int detID) { return hodoPatterns[detID-nChamberPlanes-1]; }

private:
  SRawEvent* rawEvent;
  std::vector<Hit>* hitsPtr;

  //additional cache for hit index locations, idx=0 not used, idx=max used for the total number of hits
  unsigned int idxFirst[nChamberPlanes+nHodoPlanes+nPropPlanes+nDarkPhotonPlanes+1+1];
  unsigned int hodoPatterns[nHodoPlanes];

  GeomSvc* p_geomSvc;

  //something that is larger than 1/2 of the cellWidth
  double spacing_limit[(nChamberPlanes+nHodoPlanes+nPropPlanes+nDarkPhotonPlanes)/2+1];

  static SRawEventSvc* p_rawEventSvc;
};

#endif