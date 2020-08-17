#ifndef _SQTRKPROP_H
#define _SQTRKPROP_H

#include <list>
#include <initializer_list>

#include <geom_svc/GeomSvc.h>

#include "SRawEvent.h"
#include "SRawEventSvc.h"

#include "FastTracklet.h"
#include "SQTrkSeeder.h"

#include "GFFitter.h"
#include "GFField.h"
#include "GFTrack.h"

class SQTrkProp
{
public:
  SQTrkProp(const PHField* field);
  ~SQTrkProp();

  void setPropSequence(std::initializer_list<int> seq) { sequence = seq; }

  void processSeed(Tracklet& seed);

  void propagateTo(SQGenFit::GFTrack* btrk, int detPairID, std::list<SQGenFit::GFTrackPtr>& tracklist);

  void trimTracklist(std::list<SQGenFit::GFTrackPtr>& tracklist);

private:
  SRawEventSvc* p_rawEvtSvc;
  GeomSvc*      p_geomSvc;

  SQGenFit::GFFitter* fitter;
  SQGenFit::GFField* gfield;
  std::vector<int> sequence;
};

#endif