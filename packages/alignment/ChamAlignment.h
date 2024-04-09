#ifndef _CHAMALIGNMENT_H
#define _CHAMALIGNMENT_H

#include <TString.h>
#include <fun4all/SubsysReco.h>
#include "SQMillepede.h"


class TrackletVector;

class ChamAlignment: public SubsysReco
{
public:
  ChamAlignment();
  virtual ~ChamAlignment();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  //! interfaces to millepede
  void enableEval(TString evalName) { mp->enableEval(evalName); }

private:
  int GetNodes(PHCompositeNode* topNode);

  bool acceptTrack(Tracklet* tracklet);

  TrackletVector* trackletVec;

  SQMillepede* mp;
};

#endif
