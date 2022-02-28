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

  return Fun4AllReturnCodes::EVENT_OK;
}

int ChamAlignment::process_event(PHCompositeNode* topNode)
{


  return Fun4AllReturnCodes::EVENT_OK;
}

int ChamAlignment::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
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
