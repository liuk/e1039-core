#ifndef _HODOALIGNMENT_H
#define _HODOALIGNMENT_H

#include <string>
#include <vector>
#include <TString.h>
#include <fun4all/SubsysReco.h>

class TH1D;
class TString;
class TFile;
class TTree;
class TCanvas;

class GeomSvc;
class SQHit;
class SQHitVector;
class TrackletVector;

class HodoAlignment: public SubsysReco
{
public:
  class HodoData;

public:
  HodoAlignment(const std::string& name = "HodoAlignment");
  virtual ~HodoAlignment();

  int Init(PHCompositeNode* topNode);
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  void registerHodo(TString hodoName);
  void enableEval(TString evalName);
  void setOutputName(TString outName) { textFileName = outName; }

private:
  int GetNodes(PHCompositeNode* topNode);

  bool evalEnabled;
  TString evalFileName;
  TString textFileName;
  std::vector<HodoData> alignData;
  
  GeomSvc* p_geomSvc;
  SQHitVector* hitVector;
  TrackletVector* trackletVec;

  TFile* evalFile;
  TTree* evalTree;
};

class HodoAlignment::HodoData
{
public:
  HodoData(TString name);
  ~HodoData();

  void fillHit(SQHit* hit, double p_exp);
  void calcOffset();
  TCanvas* plot();
  double getOffset();
  double findHistCenter(TH1D* hist);

public:
  TString hodoName;
  int hodoID;
  int nElements;
  double width;
  double zhodo;

  TH1D* histAll;
  std::vector<TH1D*> elemHists;

  int nValidEntries;
  double offsetAll;
  double offsetByElem;  //ofset calculated by the weighted sum of the all elements
};

#endif
