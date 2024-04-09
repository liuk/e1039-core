#include "SQMillepede.h"

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include <ktracker/FastTracklet.h>
#include <ktracker/SQTrackletFitter.h>
#include <geom_svc/GeomSvc.h>

#include "SQMPNode.h"

SQMillepede::SQMillepede():
  evalEnabled(false),
  evalNode(nullptr),
  evalFile(nullptr),
  evalTree(nullptr),
  refit_enable(false)
{
  for(int i = 0; i < nChamberPlanes; ++i) evalHist[i] = nullptr;
  trackletFitter = new SQTrackletFitter(false);
}

SQMillepede::~SQMillepede()
{
  closeEvaluation();
  delete trackletFitter;
}

void SQMillepede::initMillepede()
{
  for(int i = 0; i < NPARS; ++i)
  {
    par_align[i] = 0.;
    err_align[i] = 0.;
  }
}

void SQMillepede::globalFit()
{
  fitglo_(par_align);
  for(int i = 0; i < NPARS; ++i)
  {
    err_align[i] = errpar_(&i);
  }
}

void SQMillepede::addTrack(Tracklet* track)
{
  nodes.clear();
  for(auto iter = track->hits.begin(); iter != track->hits.end(); ++iter)
  {
    if(iter->hit.index < 0) continue;

    if(!refit_enable)
    {
      SQMPNode node(*iter, *track);
      nodes.push_back(node);
    }
    else
    {
      // create a temporary track object and disable the current hit
      Tracklet trkfit = *track;
      // std::cout << "----------------------------------------------" << std::endl;
      // trkfit.print();
      // trackletFitter->fit(trkfit);
      // trkfit.print();

      auto jter = trkfit.hits.begin();
      for(; jter != trkfit.hits.end(); ++jter)
      {
        if(iter->hit.index == jter->hit.index)
        {
          jter->hit.index = -iter->hit.index;
          break;
        }
      }
      
      trackletFitter->fit(trkfit);

      jter->hit.index = iter->hit.index;
      trkfit.calcChisq();

      SQMPNode node(*jter, trkfit);
      nodes.push_back(node);
    }
  }

  addMilleTrack();
  fillEvaluation();
}

void SQMillepede::addMilleTrack()
{
  float dergb[NPARS];
  float derlc[NPARTRK];
  float meas, sigma;

  for(auto node = nodes.begin(); node != nodes.end(); ++node)
  {
    if(node->isValid() == 0) continue;

    //Initialize all parameters
    zerloc_(dergb, derlc);

    //Get measurements and derivatives
    int idx = node->detectorID - 1;
    meas = node->meas;
    sigma = node->sigma;

    derlc[0] = node->dwdx;
    derlc[1] = node->dwdy;
    derlc[2] = node->dwdtx;
    derlc[3] = node->dwdty;

    dergb[NPARPLAN*idx + 0] = node->dwdz;
    dergb[NPARPLAN*idx + 1] = node->dwdphi;
    dergb[NPARPLAN*idx + 2] = node->dwdw;

    //Fill the numbers to millepede
    equloc_(dergb, derlc, &meas, &sigma);
  }

  //Perform a local fit
  fitloc_();
}

void SQMillepede::addConstraint(float par[], float val)
{
  float rhs = val;
  float dercs[NPARS];
  for(int i = 0; i < NPARS; ++i) dercs[i] = par[i];

  constf_(dercs, &rhs);
}

void SQMillepede::setDetParameter(int detectorID, int parID, float val)
{
  par_align[globalId(detectorID, parID)] = val;
}

void SQMillepede::setDetParaError(int detectorID, int parID, float val)
{
  err_align[globalId(detectorID, parID)] = val;
}

void SQMillepede::fixDetParameter(int detectorID, int parID, float val)
{
  setDetParameter(detectorID, parID, val);
  setDetParaError(detectorID, parID, 0.);

  //Add constraint to ensure the parameter won't change
  float dercs[NPARS];
  float rhs = val;

  for(int i = 0; i < NPARS; ++i) dercs[i] = 0.;
  dercs[globalId(detectorID, parID)] = 1.;

  //pass constraints to millepede
  constf_(dercs, &rhs);
}

void SQMillepede::fixDetectorPair(int detID1, int detID2, int parID)
{
  float rhs = 0.;
  float dercs[NPARS];
  for(int i = 0; i < NPARS; ++i) dercs[i] = 0.;

  dercs[globalId(detID1, parID)] = 1.;
  dercs[globalId(detID2, parID)] = -1.;

  constf_(dercs, &rhs);
}

void SQMillepede::enableEval(TString evalFileName)
{
  evalEnabled = true;

  evalNode = new SQMPNode(0);
  evalFile = new TFile(evalFileName, "recreate");
  evalTree = new TTree("eval", "eval");
  evalTree->Branch("MPNode", &evalNode, 256000, 99);

  GeomSvc* p_geomSvc = GeomSvc::instance();
  for(int i = 0; i < nChamberPlanes; ++i)
  {
    evalHist[i] = new TH1D(p_geomSvc->getDetectorName(i+1).c_str(), p_geomSvc->getDetectorName(i+1).c_str(), 100, -0.5*p_geomSvc->getPlaneSpacing(i+1), 0.5*p_geomSvc->getPlaneSpacing(i+1));
  }
}

void SQMillepede::fillEvaluation()
{
  if(!evalEnabled) return;

  for(auto node = nodes.begin(); node != nodes.end(); ++node)
  {
    if(node->isValid() == 0) continue;

    *evalNode = *node;
    evalTree->Fill();

    if(node->detectorID > 0) evalHist[node->detectorID-1]->Fill(node->meas);
  }
}

void SQMillepede::closeEvaluation()
{
  if(!evalEnabled) return;

  evalFile->cd();
  evalTree->Write();
  for(int i = 0; i < nChamberPlanes; ++i) evalHist[i]->Write();
  evalFile->Close();

  delete evalNode;
}
