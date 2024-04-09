#ifndef SQMILLEPEDE_H
#define SQMILLEPEDE_H

#include <vector>
#include <GlobalConsts.h>

class TFile;
class TTree;
class TH1D;
class TString;
class SQMPNode;
class Tracklet;
class SQTrackletFitter;

///Interface for cross-compilation with fortran millepede program
extern "C"
{
  //Initialize global alignment parameter arrays
  extern void initgl_(int*, int*, int*, int*);

  //Optional: define sigma for single parameter
  extern void parsig_(int*, float*);

  //Optional: unit for iterations
  extern void initun_(int*, float*);

  //Optional: set initial value of global parameters
  extern void parglo_(float*);

  //Optional: constraints
  extern void constf_(float*, float*);

  //Set all local arrays to zero
  extern void zerloc_(float*, float*);

  //Equations for local fit
  extern void equloc_(float*, float*, float*, float*);

  //Local parameter fit (+entry KILLOC)
  extern void fitloc_();

  //Global parameter fit
  extern void fitglo_(float*);

  //Retrieve error for a given parameter
  extern float errpar_(int*);
}

// number of parameters per track
#define NPARTRK  4

// number of alignment parameters per detector plane
#define NPARPLAN 3

// number of maximum parameters Millepede will handle
#define NPARS NPARPLAN*nChamberPlanes

class SQMillepede
{
public:
  SQMillepede();
  ~SQMillepede();

  void initMillepede();

  void globalFit();

  void addTrack(Tracklet* track);
  void addMilleTrack();

  void addConstraint(float par[], float val);

  void setDetParameter(int detectorID, int parID, float val);
  void setDetParaError(int detectorID, int parID, float val);
  void fixDetParameter(int detectorID, int parID, float val);
  void fixDetectorPair(int detID1, int detID2, int parID);

  float getDetParameter(int detectorID, int parID) { return par_align[globalId(detectorID, parID)]; }
  float getDetParaError(int detectorID, int parID) { return err_align[globalId(detectorID, parID)]; }

  void enableEval(TString evalFileName);
  void fillEvaluation();
  void closeEvaluation();

  void enableRefit() { refit_enable = true; }

  TH1D* getEvalHist(int detectorID) { return evalHist[detectorID-1]; }

private:
  int globalId(int detectorID, int parID) { return (detectorID-1)*NPARPLAN + parID; }

private:
  // array of alignment parameters
  float par_align[NPARS];
  float err_align[NPARS];

  // intermediate container
  std::vector<SQMPNode> nodes;

  // whether we remove the hit of interest
  bool refit_enable;

  // tracklet fitter
  SQTrackletFitter* trackletFitter;

  // evaluation output
  bool evalEnabled;
  SQMPNode* evalNode;
  TFile* evalFile;
  TTree* evalTree;
  TH1D*  evalHist[nChamberPlanes];
};

#endif
