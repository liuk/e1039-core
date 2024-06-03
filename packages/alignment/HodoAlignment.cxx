#include "HodoAlignment.h"

#include <memory>
#include <fstream>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TText.h>

#include <geom_svc/GeomSvc.h>
#include <interface_main/SQHitVector_v1.h>
#include <UtilAna/UtilSQHit.h>
#include <ktracker/FastTracklet.h>


HodoAlignment::HodoAlignment(const std::string& name):
  SubsysReco(name),
  textFileName("alignment_hodo.txt"),
  evalEnabled(false),
  evalFile(nullptr),
  evalTree(nullptr),
  p_geomSvc(nullptr),
  hitVector(nullptr),
  trackletVec(nullptr)
{}

HodoAlignment::~HodoAlignment()
{
}

int HodoAlignment::Init(PHCompositeNode* topNode)
{
  registerHodo("H1B");
  registerHodo("H1T");
  registerHodo("H1L");
  registerHodo("H1R");
  registerHodo("H2B");
  registerHodo("H2T");
  registerHodo("H2L");
  registerHodo("H2R");
  registerHodo("H3B");
  registerHodo("H3T");
  registerHodo("H4Y1L");
  registerHodo("H4Y1R");
  registerHodo("H4Y2L");
  registerHodo("H4Y2R");
  registerHodo("H4B");
  registerHodo("H4T");

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
  int nTracklets = trackletVec->size();
  if(nTracklets < 1) return Fun4AllReturnCodes::EVENT_OK;

  Tracklet* tracklet = trackletVec->at(0);   //only use the first tracklet for simplicity
  for(unsigned int i = 0; i < alignData.size(); ++i)
  {
    double x_exp = tracklet->getExpPositionX(alignData[i].zhodo);
    double y_exp = tracklet->getExpPositionY(alignData[i].zhodo);
    double p_exp = p_geomSvc->getUinStereoPlane(alignData[i].hodoID, x_exp, y_exp);
    if(!p_geomSvc->isInPlane(alignData[i].hodoID, x_exp, y_exp)) continue;

    std::shared_ptr<SQHitVector> hodoHits(UtilSQHit::FindHits(hitVector, alignData[i].hodoID));

    SQHit* hit1 = nullptr;
    SQHit* hit2 = nullptr;
    if(hodoHits->size() == 1)
    {
      hit1 = hodoHits->at(0);
    }
    else if(hodoHits->size() == 2)
    {
      hit1 = hodoHits->at(0);
      hit2 = hodoHits->at(1);
      if(abs(hit1->get_element_id() - hit2->get_element_id()) != 1) continue;
    }
    else
    {
      continue;
    }

    alignData[i].fillHit(hit1, p_exp);
    if(hit2 != nullptr) alignData[i].fillHit(hit2, p_exp);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int HodoAlignment::End(PHCompositeNode* topNode)
{
  std::fstream fout(textFileName.Data(), std::ios::out);
  if(evalEnabled)
  {
    evalFile = new TFile(evalFileName, "recreate");
  }

  for(unsigned int i = 0; i < alignData.size(); ++i)
  {
    alignData[i].calcOffset();
    fout << alignData[i].getOffset() << "    " << alignData[i].hodoID << std::endl;

    TCanvas* c = alignData[i].plot();
    c->Write();
  }
  fout.close();
  if(evalEnabled) evalFile->Close();

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
  evalFileName = evalName;
}

HodoAlignment::HodoData::HodoData(TString name)
{
  GeomSvc* p_geomSvc = GeomSvc::instance();

  hodoName = name;
  hodoID = p_geomSvc->getDetectorID(hodoName.Data());
  nElements = p_geomSvc->getPlaneNElements(hodoID);
  width = p_geomSvc->getCellWidth(hodoID);
  zhodo = p_geomSvc->getPlanePosition(hodoID);

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

void HodoAlignment::HodoData::fillHit(SQHit* hit, double p_exp)
{
  int elementID = hit->get_element_id();
  double p = hit->get_pos();
  if(fabs(p_exp - p) < width)
  {
    histAll->Fill(p_exp - p);
    elemHists[elementID-1]->Fill(p_exp - p);
  }
}

double HodoAlignment::HodoData::getOffset()
{
  return offsetAll;
}

HodoAlignment::HodoData::~HodoData()
{
  delete histAll;
  for(int i = 0; i < nElements; ++i)
  {
    delete elemHists[i];
  }
}

void HodoAlignment::HodoData::calcOffset()
{
  //calculate offset by all elements combined
  offsetAll = findHistCenter(histAll);

  //calculate offset by each individual elements
  offsetByElem = 0.;
  int nValidEntries = 0.;
  for(unsigned int j = 0; j < nElements; ++j)
  {
    int nEntries = elemHists[j]->GetEntries();
    if(nEntries < 500)
    {
      elemHists[j]->Rebin(2);
    }

    double offset = findHistCenter(elemHists[j]);
    if(offset < 100.)
    {
      offsetByElem  += offset*nEntries;
      nValidEntries += nEntries;
    }
    offsetByElem = offsetByElem/nValidEntries;
  }
}

TCanvas* HodoAlignment::HodoData::plot()
{
  TCanvas* c = new TCanvas();
  histAll->Draw();

  double y_max = histAll->GetBinContent(histAll->GetMaximumBin());
  TArrow* ar_center = new TArrow(offsetAll, y_max, offsetAll, 0., 0.01, ">");
  TArrow* ar_left   = new TArrow(offsetAll - 0.5*width, y_max, offsetAll - 0.5*width, 0., 0.01, ">");
  TArrow* ar_right  = new TArrow(offsetAll + 0.5*width, y_max, offsetAll + 0.5*width, 0., 0.01, ">");

  ar_center->SetLineWidth(2);
  ar_left->SetLineWidth(2);
  ar_right->SetLineWidth(2);

  ar_center->SetLineColor(kRed);
  ar_right->SetLineColor(kBlue);
  ar_left->SetLineColor(kBlue);

  ar_center->Draw();
  ar_left->Draw();
  ar_right->Draw();

  TText* text = new TText();
  text->DrawTextNDC(0.1, 0.92, TString::Format("Offsets: %f cm", offsetAll));

  return c;
}

double HodoAlignment::HodoData::findHistCenter(TH1D* hist)
{
  //Basic quality cut
  if(hist->GetEntries() < 150) return 9999.;
  if(hist->GetRMS() < width/5. && hist->GetEntries() < 1000) return 9999.;

  int nBin = hist->GetNbinsX();
  int nBinInSize = nBin/2;
  double binWidth = hist->GetBinWidth(1);

  //Search from left side
  int nEvt_max = 0;
  int index_max_left = 0;
  for(int i = 1; i <= nBinInSize; ++i)
  {
    int nEvt_curr = int(hist->Integral(i, i + nBinInSize));
    //cout << i << " : " << hist->GetBinCenter(i) << " <===> " << hist->GetBinCenter(i + nBinInSize) << " : " << nEvt_curr << " === " << nEvt_max << endl;
    if(nEvt_curr > nEvt_max)
    {
        nEvt_max = nEvt_curr;
        index_max_left = i;
    }
  }

  //Search from right side
  nEvt_max = 0;
  int index_max_right = nBin;
  for(int i = nBin; i >= nBinInSize; --i)
  {
    int nEvt_curr = int(hist->Integral(i - nBinInSize, i));
    if(nEvt_curr > nEvt_max)
    {
        nEvt_max = nEvt_curr;
        index_max_right = i;
    }
  }

  //Left-right difference should not be too large
  double left_center = hist->GetBinCenter(index_max_left + nBinInSize/2);
  double right_center = hist->GetBinCenter(index_max_right - nBinInSize/2);
  
  if(fabs(left_center - right_center) > 5.*binWidth) return 9999.;
  return (left_center + right_center)/2.;
}
