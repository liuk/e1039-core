#ifndef SQTriggerAnalyzer_H
#define SQTriggerAnalyzer_H

#include <TString.h>
#include <TFormula.h>

#include <fun4all/SubsysReco.h>
#include <interface_main/SQEvent_v1.h>

#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <list>
#include <map>
#include <string>
#include <initializer_list>

class PHG4Hit;
class SQRun;
class SQSpillMap;
class SQEvent;
class SQHit;
class SQHitMap;
class SQHitVector;
class PHG4HitContainer;

class GeomSvc;

class SQTriggerRoad
{
public:
  SQTriggerRoad(const unsigned short N = 4);
  SQTriggerRoad(const std::list<int>& path);

  //!set the trigger element uniqueID
  void setTrID(const int i, const int uID) { uniqueTrIDs[i] = uID; }
  void setTrIDs(const int uIDs[]) { for(unsigned int i = 0; i < nTrHits; ++i) uniqueTrIDs[i] = uIDs[i]; }

  //!Get the sign of LR or TB
  int getTB() const;

  //!flip the LR or TB
  void flipTB();

  //! get the TriggerPrimitive name, i.e. top plus, bottom minus, etc
  TString getPrimitiveName() const;

  //!Other gets
  //@{
  int getRoadID() const { return roadID; }
  int getCharge() const { return charge; }
  double getSigWeight() const { return sigWeight; }
  double getBkgRate() const { return bkgRate; }
  double getPxMin() const { return pXmin; }
  int getTrID(unsigned int i) const { return i < uniqueTrIDs.size() ? uniqueTrIDs[i] : 0; }
  int getTrDetectorID(unsigned int i) const { return getTrID(i)/1000; }
  int getTrElementID(unsigned int i) const { return getTrID(i) % 1000; }
  TString getStringID() const;
  //@}

  //!Sets
  //@{
  void setRoadID(int id) { roadID = id; }
  void setCharge(int ch) { charge = ch; }
  void setSigWeight(double weight) { sigWeight = weight; }
  void setBkgRate(double rate) { bkgRate = rate; }
  void setPxMin(double pxmin) { pXmin = pxmin; }
  //@}

  //!comparison
  //@{
  bool operator == (const SQTriggerRoad& elem) const;
  bool operator <  (const SQTriggerRoad& elem) const;
  //@}

  //!printer
  friend std::ostream& operator << (std::ostream& os, const SQTriggerRoad& road);

private:
  //!unique road ID
  int roadID;

  //!charge
  int charge;

  //!total signal weight
  double sigWeight;

  //!total background occurence
  double bkgRate;

  //!Minimum Px
  double pXmin;

  //!unique detector element IDs: = 1000*detectorID + elementID
  std::vector<int> uniqueTrIDs;

  //!number of hits per road
  const unsigned short nTrHits;
};

class SQTriggerAnalyzer: public SubsysReco
{
public:
  //!Forward declaration of MatrixNode format -- TODO will be eventuall private
  class MatrixNode;
  class TriggerPrimitive;
  class TriggerCondition;

public:
  SQTriggerAnalyzer(const std::string& name = "SQTriggerAnalyzer", const unsigned short N = 4, const int verbosity = 0);
  virtual ~SQTriggerAnalyzer();

#ifndef __CINT__
  int Init(PHCompositeNode* topNode);
#endif
    
  //! fun4all interfaces
  int InitRun(PHCompositeNode* topNode);
  int process_event(PHCompositeNode* topNode);
  int End(PHCompositeNode* topNode);

  //! register trigger conditions
  void registerTrigger(const SQEvent::TriggerMask bit, const std::string& expr);

  //! Initalize default trigger configuration
  void initDefault();
  void initDefaultMatrix();
  void initDefaultNIM();

  //! update the hit flags for hodos, return number of hits updated
  int updateHitFlags();

  void set_matrix_file(const std::string& filename) { matrixFileName = filename; }

  //!Build the trigger matrix by the input roads list
  bool loadMatrixFile();
  void buildMatrix();

  //!create the trigger stations hit pattern, return false if one or more plane is missing
  bool buildHitPattern();

  //!search for possible roads
  void searchMatrix(MatrixNode* node, int level, int index);

  //!Tree deletion
  void deleteMatrix(MatrixNode* node);

  //!Helper function to retrieve the found road list
  std::list<SQTriggerRoad>& getRoadsFound(int index) { return roads_found[index]; }

  //!Helper functions to print various things
  void printHitPattern();
  void printPath();
  void print(std::ostream& os = std::cout) const;

private:
  //! road set input file name
  std::string matrixFileName;

  //!Internal hit pattern structure
  typedef std::vector<std::set<int> > TrHitPattern;
  TrHitPattern data;

  //!the trigger matrix, 0 for mu+, 1 for mu-
  //@{
  MatrixNode* matrix[2];
  std::map<TString, SQTriggerRoad> roads[2];
  //@}

  //!container of the roads found for +/-
  std::list<SQTriggerRoad> roads_found[2];

  //!temporary container of traversal path
  std::list<int> path;

  //! FPGA trigger detector level/depth mapping
  std::map<int, int> matrixLevelByID;

  //!flag of successful matrix loading from file
  bool MatrixAvailable;

  //! IO nodes
  SQEvent*     _event_header;
  SQHitVector* _hit_vector;

  //!pointer to the digitizer, or geometry for that matter
  GeomSvc* p_geomSvc;

  //!number of hits needed per TriggerRoad
  const unsigned short nTrHits;

  //!Trigger primitive maps
  std::map<std::string, TriggerPrimitive*> trPrimByName;

  //!Trigger conditions
  std::map<SQEvent::TriggerMask, TriggerCondition> trigConds;
};

class SQTriggerAnalyzer::MatrixNode
{
public:
  MatrixNode(int uID);

  //!add a child
  void add(MatrixNode* child);

  //!get detectorID
  int detectorID() { return uniqueID/1000; }

  //!get elementID
  int elementID() { return uniqueID % 1000; }

public:
  //! UniqueID = 1000*detectorID + elementID
  int uniqueID;
  std::map<int, MatrixNode*> children;
};

class SQTriggerAnalyzer::TriggerPrimitive
{
public:
  TriggerPrimitive(TString rawName);

public:
  std::string name;
  int lo;
  int hi;
  int counter;
};

class SQTriggerAnalyzer::TriggerCondition
{
public:
  bool eval();
  void reset();

public:
  SQEvent::TriggerMask trigBit;
  TFormula trigExpr;

  std::vector<TriggerPrimitive*> trigPrims;
};

#endif
