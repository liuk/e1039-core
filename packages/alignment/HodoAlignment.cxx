#include "HodoAlignment.h"
#include <fun4all/Fun4AllReturnCodes.h>

HodoAlignment::HodoAlignment()
{}

HodoAlignment::~HodoAlignment()
{}

int HodoAlignment::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int HodoAlignment::InitRun(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int HodoAlignment::process_event(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int HodoAlignment::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int HodoAlignment::GetNodes(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
