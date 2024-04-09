#include "HodoAlignment.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>

#include <geom_svc/GeomSvc.h>
#include <interface_main/SQHitVector_v1.h>
#include <ktracker/FastTracklet.h>


HodoAlignment::HodoAlignment(const std::string& name):
  SubsysReco(name),
  evalEnabled(false),
  p_geomSvc(nullptr),
  hitVector(nullptr),
  trackletVec(nullptr)
{}

HodoAlignment::~HodoAlignment()
{
}

int HodoAlignment::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int HodoAlignment::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

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
  trackletVec = findNode::getClass<TrackletVector>(topNode, "TrackletVector");
  if(!trackletVec) 
  {
    std::cout << "ERROR: cannot get tracklet vector needed for hodo alignment, will quit now" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  hitVector = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if(!hitVector)
  {
    std::cout << "ERROR: cannot get hit vector needed for chamber alignment, will quit now" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void HodoAlignment::registerHodo(TString hodoName)
{
  alignData.push_back(HodoData(hodoName));
}

void HodoAlignment::enableEval(TString evalName)
{
  evalEnabled = true;
  evalPrefix = evalName;
}

HodoAlignment::HodoData::HodoData(TString name)
{
  GeomSvc* p_geomSvc = GeomSvc::instance();

  hodoName = name;
  hodoID = p_geomSvc->getDetectorID(hodoName);
  nElements = p_geomSvc->getPlaneNElements(hodoID);
  width = p_geomSvc->getCellWidth(hodoID);

  histAll = new TH1D(hodoName, hodoName, 200, -width, width);
  for(int i = 0; i < nElements; ++i)
  {
    TH1D* htemp = new TH1D(TString::Format("%s_%i", hodoName, i+1), TString::Format("%s_%i", hodoName, i+1), 200, -width, width);
    elemHists.push_back(htemp);
  }

  nValidEntries = 0;
  offsetAll = 0.;
  offsetByElem = 0.;
}