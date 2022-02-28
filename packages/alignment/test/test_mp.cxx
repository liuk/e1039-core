#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <GlobalConsts.h>
#include <ktracker/FastTracklet.h>
#include <geom_svc/GeomSvc.h>
#include <phool/recoConsts.h>
#include "SQMillepede.h"

using namespace std;

int main(int argc, char* argv[])
{
  recoConsts* rc = recoConsts::instance();
  rc->set_DoubleFlag("FMAGSTR", -1.054);
  rc->set_DoubleFlag("KMAGSTR", 0.);
  // rc->set_CharFlag("AlignmentMille", "$E1039_RESOURCE/alignment/run3/align_mille.txt");
  // rc->set_CharFlag("AlignmentHodo", "$E1039_RESOURCE/alignment/run3/alignment_hodo.txt");
  // rc->set_CharFlag("AlignmentProp", "$E1039_RESOURCE/alignment/run3/alignment_prop.txt");
  // rc->set_CharFlag("Calibration", "$E1039_RESOURCE/alignment/run3/calibration.txt");
  // rc->set_CharFlag("Geometry", "geometry_G17_run3");
  rc->set_CharFlag("MySQLURL", "mysql://localhost");
  rc->set_BoolFlag("KMAG_ON", false);
  rc->Print();

  GeomSvc::UseDbSvc(true);
  GeomSvc* p_geomSvc = GeomSvc::instance();
  p_geomSvc->printTable();

  TClonesArray* tracklets = new TClonesArray("Tracklet");
  tracklets->Clear();

  TFile* dataFile = new TFile(argv[1], "read");
  TTree* dataTree = (TTree*)dataFile->Get("eval");

  dataTree->SetBranchAddress("tracklets", &tracklets);

  SQMillepede* mp = new SQMillepede();
  mp->initMillepede();
  mp->enableEval("alignment_eval.root");

  int nEvents = dataTree->GetEntries();
  for(int i = 0; i < nEvents; ++i)
  {
    dataTree->GetEntry(i);

    int nTracklets = tracklets->GetEntries();
    for(int j = 0; j < nTracklets; ++j)
    {
      cout << i << "   " << j << endl;
      Tracklet* tracklet = (Tracklet*)tracklets->At(j);
      tracklet->print();
      mp->addTrack(tracklet);
    }

    tracklets->Clear();
  }

  // mp->globalFit();
  for(int i = 1; i <= nChamberPlanes; ++i)
  {
    cout << i << "  ";
    for(int j = 0; j < 3; ++j) cout << mp->getDetParameter(i, j) << " +/- " << mp->getDetParaError(i, j) << "  ";
    cout << endl;
  }

  delete mp;
  return 0;
}