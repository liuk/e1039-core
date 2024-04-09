#include "ChamAlignment.h"

#include <iostream>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <ktracker/FastTracklet.h>

ChamAlignment::ChamAlignment()
{
  mp = new SQMillepede();
}

ChamAlignment::~ChamAlignment()
{
  if(mp != nullptr) delete mp;
}

int ChamAlignment::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int ChamAlignment::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  //Initialize millepede fitter
  //Set initial values here
  for(int i = 1; i <= nChamberPlanes; ++i)
  {
    mp->setDetParaError(i, 0, 0.1);
    mp->setDetParaError(i, 1, 0.005);
    mp->setDetParaError(i, 2, 0.05);
  }

  //Fix D2X to 

  return Fun4AllReturnCodes::EVENT_OK;
}

int ChamAlignment::process_event(PHCompositeNode* topNode)
{
  int nTracklets = trackletVec->size();
  for(int i = 0; i < nTracklets; ++i)
  {
    Tracklet* tracklet = trackletVec->at(i);
    if(!acceptTrack(tracklet)) return;

    mp->addTrack(tracklet);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int ChamAlignment::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

bool ChamAlignment::acceptTrack(Tracklet* tracklet)
{
  return true;
}

int ChamAlignment::GetNodes(PHCompositeNode* topNode)
{
  trackletVec = findNode::getClass<TrackletVector>(topNode, "TrackletVector");
  if(!trackletVec) 
  {
    std::cout << "ERROR: cannot get tracklet vector needed for chamber alignment, will quit now" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
