#include "SRawEventSvc.h"

#include <bitset>

//#define _DEBUG_ON
#ifdef _DEBUG_ON
#  define LogDebug(exp) std::cout << "DEBUG: " << typeid(*this).name() << " " << __FUNCTION__ << " " << __LINE__ << " :: " << exp << std::endl
#else
#  define LogDebug(exp)
#endif

SRawEventSvc* SRawEventSvc::p_rawEventSvc = nullptr;
SRawEventSvc* SRawEventSvc::instance()
{
  if(p_rawEventSvc == nullptr)
  {
    p_rawEventSvc = new SRawEventSvc();
  }
  return p_rawEventSvc;
}

SRawEventSvc::SRawEventSvc():
  rawEvent(nullptr),
  hitsPtr(nullptr),
  p_geomSvc(GeomSvc::instance())
{
  spacing_limit[0] = 0.;
  for(int i = 1; i <= (nChamberPlanes+nHodoPlanes+nPropPlanes+nDarkPhotonPlanes)/2; ++i)
  {
    spacing_limit[i] = p_geomSvc->getPlaneSpacing(2*i-1);
  }

  idxFirst[0] = 0;
}

bool SRawEventSvc::setRawEvent(SRawEvent* r, bool buildAuxInfo)
{
  if(!acceptEvent(r)) return false;

  rawEvent = r;
  hitsPtr  = &(r->getAllHits());
  if(!buildAuxInfo) return true;

  std::vector<Hit>& hits = *hitsPtr;

  int detectorID_prev = 0;
  for(int i = 0; i < nHodoPlanes; ++i) hodoPatterns[i] = 0;
  for(unsigned int i = 0; i < hits.size(); ++i)
  {
    if(hits[i].detectorID > detectorID_prev)
    {
      //idxFirst[hits[i].detectorID] = i;
      for(int j = detectorID_prev+1; j <= hits[i].detectorID; ++j) idxFirst[j] = i;
      detectorID_prev = hits[i].detectorID;
    }

    if(hits[i].detectorID > nChamberPlanes && hits[i].detectorID <= nChamberPlanes+nHodoPlanes)
    {
      hodoPatterns[hits[i].detectorID-nChamberPlanes-1] |= (1 << hits[i].elementID);
    }
  }
  for(int i = detectorID_prev+1; i < nChamberPlanes+nHodoPlanes+nPropPlanes+nDarkPhotonPlanes+1+1; ++i) idxFirst[i] = rawEvent->fNHits[0];

  for(int i = 0; i < nChamberPlanes+nHodoPlanes+nPropPlanes+nDarkPhotonPlanes+1; ++i) LogDebug(i << " - " << idxFirst[i]);
  for(int i = 0; i < nHodoPlanes; ++i) LogDebug(i << "  " << std::bitset<32>(hodoPatterns[i]));
  return true;
}

bool SRawEventSvc::acceptEvent(SRawEvent* r)
{
  return true;
}

bool SRawEventSvc::hodoMasked(int detID, int eleID, int win)
{
  if(eleID < 1 || eleID > p_geomSvc->getPlaneNElements(detID)) return false;
  if(win < 0) return false;

  unsigned int pattern = hodoPattern(detID);
  if(pattern == 0) return false;

  unsigned int core = (1 << (2*win+1)) - 1;
  if(eleID > win)
  {
    core = core << (eleID - win);
  }
  else
  {
    core = core >> (win - eleID);
  }
  
  LogDebug("Matching core = " << std::bitset<32>(core) << ", hit pattern = " << std::bitset<32>(pattern));
  return (pattern & core) != 0;
}

bool SRawEventSvc::hodoMasked(int detID, double pos, int win)
{
  int elementID = p_geomSvc->getExpElementID(detID, pos);
  LogDebug("Expected elementID on " << detID << ": " << elementID);
  return hodoMasked(detID, elementID, win);
}

std::list<int> SRawEventSvc::getHitsIndexInDetector(short detectorID, double p_exp, double win)
{
  std::list<int> hit_list;
  for(int i = idxFirst[detectorID]; i < idxFirst[detectorID+1]; ++i)
  {
    if(win < 998. && fabs(p_geomSvc->getMeasurement(hit(i).detectorID, hit(i).elementID) - p_exp) > win) continue;
    hit_list.push_back(i);
  }

  return hit_list;
}

std::list<int> SRawEventSvc::getHitsIndexInDetectorPair(short pairID)
{
  std::list<int> hit_list;
  for(int i = idxFirst[2*pairID-1]; i < idxFirst[2*pairID+1]; ++i)
  {
    hit_list.push_back(i);
  }

  return hit_list;
}

std::list<int> SRawEventSvc::getHitsIndexInDetectors(std::vector<int>& detectorIDs)
{
  std::list<int> hit_list;
  unsigned int nDetectors = detectorIDs.size();
  for(unsigned int i = 0; i < nDetectors; ++i)
  {
    for(int j = idxFirst[detectorIDs[i]]; j < idxFirst[detectorIDs[i]+1]; ++j) hit_list.push_back(j);
  }

  return hit_list;
}

std::list<SRawEventSvc::HitPair_t> SRawEventSvc::getHitPairsInDetectorPair(short pairID, double p_exp, double win)
{
  std::list<HitPair_t> hitpairs;
  std::list<int> hitlist1 = getHitsIndexInDetector(2*pairID,   p_exp, win);
  std::list<int> hitlist2 = getHitsIndexInDetector(2*pairID-1, p_exp, win+spacing_limit[pairID]);

  std::vector<int> hitflag1(hitlist1.size(), -1);   //-1 means out-of-the window, 0 means in the window but not paired, 1 means used in pair(s)
  std::vector<int> hitflag2(hitlist2.size(), -1);  

  LogInfo(hitlist1.size() << "  " << hitlist2.size());

  int index1 = -1;
  int index2 = -1;
  for(auto iter = hitlist1.begin(); iter != hitlist1.end(); ++iter)
  {
    ++index1;
    double pos1 = p_geomSvc->getMeasurement(hit(*iter).detectorID, hit(*iter).elementID);
    //if(win < 998. && fabs(pos1 - p_exp) > win) continue;

    hitflag1[index1] = 0;
    index2 = -1;
    for(auto jter = hitlist2.begin(); jter != hitlist2.end(); ++jter)
    {
      ++index2;
      double pos2 = p_geomSvc->getMeasurement(hit(*jter).detectorID, hit(*jter).elementID);
      //if(win < 998. && fabs(pos2 - p_exp) > win) continue;
      
      hitflag2[index2] = 0;
      if(fabs(pos1 - pos2) > spacing_limit[pairID]) continue;

      hitpairs.push_back(std::make_pair(*iter, *jter));
      hitflag1[index1] = 1;
      hitflag2[index2] = 1;
    }
  }

  index1 = 0;
  for(auto iter = hitlist1.begin(); iter != hitlist1.end(); ++iter)
  {
    if(hitflag1[index1] == 0) hitpairs.push_back(std::make_pair(*iter, -1));
    ++index1;
  }

  index2 = 0;
  for(auto iter = hitlist2.begin(); iter != hitlist2.end(); ++iter)
  {
    if(hitflag2[index2] == 0) hitpairs.push_back(std::make_pair(*iter, -1));
    ++index2;
  }

  return hitpairs;
}