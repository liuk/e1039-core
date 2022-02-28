#include "PropAlignment.h"
#include <fun4all/Fun4AllReturnCodes.h>

PropAlignment::PropAlignment()
{}

PropAlignment::~PropAlignment()
{}

int PropAlignment::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PropAlignment::InitRun(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PropAlignment::process_event(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PropAlignment::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PropAlignment::GetNodes(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
