#include "SQTriggerAnalyzer.h"

#include <interface_main/SQHit.h>
#include <interface_main/SQHit_v1.h>
#include <interface_main/SQHitMap_v1.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQEvent_v1.h>
#include <interface_main/SQSpillMap_v1.h>
#include <geom_svc/GeomSvc.h>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <TObjArray.h>
#include <TObjString.h>

#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>

SQTriggerRoad::SQTriggerRoad(const unsigned short N): 
  roadID(0), 
  sigWeight(0.), 
  bkgRate(0.), 
  pXmin(0.),
  nTrHits(N)
{
  uniqueTrIDs.clear();
  for(unsigned int i = 0; i < nTrHits; ++i) uniqueTrIDs.push_back(0);
}

SQTriggerRoad::SQTriggerRoad(const std::list<int>& path): 
  roadID(0), 
  sigWeight(0.), 
  bkgRate(0.), 
  pXmin(0.),
  nTrHits(path.size())
{
  uniqueTrIDs.clear();
  for(std::list<int>::const_iterator iter = path.begin(); iter != path.end(); ++iter)
  {
    if(*iter < 0) continue;
    uniqueTrIDs.push_back(*iter);
  }
}

int SQTriggerRoad::getTB() const
{
  int TB = 0;
  for(unsigned int i = 0; i < nTrHits; ++i)
  {
    int detectorID = getTrDetectorID(i);
    TB = TB + (((detectorID & 1) == 0) ? 1 : -1);
  }
  
  if(TB == 4)  return 1;
  if(TB == -4) return -1;

  return 0;
}

void SQTriggerRoad::flipTB()
{
  int corr = getTB();
  for(unsigned int i = 0; i < nTrHits; ++i)
  {
    int detectorID = getTrDetectorID(i);
    int elementID  = getTrElementID(i);
    
    detectorID -= corr;
    uniqueTrIDs[i] = detectorID*1000 + elementID;
  }
}

TString SQTriggerRoad::getPrimitiveName() const
{
  return TString(Form("%s%s", getTB() == 1 ? "T" : "B", charge > 0 ? "P" : "M"));
}

bool SQTriggerRoad::operator == (const SQTriggerRoad& elem) const
{
  for(unsigned int i = 0; i < nTrHits; ++i)
  {
    if(uniqueTrIDs[i] != elem.uniqueTrIDs[i]) return false;
  }
  
  return true;
}

bool SQTriggerRoad::operator < (const SQTriggerRoad& elem) const
{
  return sigWeight < elem.sigWeight;
}

TString SQTriggerRoad::getStringID() const
{
  TString sid;
  for(unsigned int i = 0; i < nTrHits; ++i)
  {
    sid = sid + Form("%06d", uniqueTrIDs[i]);
  }
  
  return sid;
}

std::ostream& operator << (std::ostream& os, const SQTriggerRoad& road)
{
  os << "Trigger Road ID = " << road.roadID << ", signal = " << road.sigWeight << ", bkg = " << road.bkgRate << "\n   ";
  os << road.getStringID();
 
  return os;
}

SQTriggerAnalyzer::SQTriggerAnalyzer(const std::string& name, const unsigned short N, const int verbosity):
  SubsysReco(name),
  nTrHits(N),
  MatrixAvailable(false),
  matrixFileName(""),
  _event_header(nullptr),
  _hit_vector(nullptr),
  p_geomSvc(nullptr)
{
  Verbosity(verbosity);
}

SQTriggerAnalyzer::~SQTriggerAnalyzer()
{
  if(MatrixAvailable)
  {
    deleteMatrix(matrix[0]);
    deleteMatrix(matrix[1]);
  }

  for(auto iter = trPrimByName.begin(); iter != trPrimByName.end(); ++iter)
  {
    delete iter->second;
  }
}

int SQTriggerAnalyzer::Init(PHCompositeNode* topNode) 
{
  p_geomSvc = GeomSvc::instance();

  //Load the trigger roads
  MatrixAvailable = loadMatrixFile();
  if(!MatrixAvailable)
  {
    std::cout << Name() << ": no matrix roads loaded, will only emulate NIM triggers" << std::endl;
  }

  //build the search matrix
  if(MatrixAvailable) buildMatrix();

  //init default trigger configuration, if user has initialized something it won't be overwritten
  initDefault();

  if(Verbosity() > Fun4AllBase::VERBOSITY_SOME) print();

  return Fun4AllReturnCodes::EVENT_OK;
}

int SQTriggerAnalyzer::InitRun(PHCompositeNode* topNode) 
{
  //book output nodes
  PHNodeIterator iter(topNode);

  PHCompositeNode* eventNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if(!eventNode) 
  {
    std::cout << Name() << ": no DST node, create one" << std::endl;
    eventNode = new PHCompositeNode("DST");
    topNode->addNode(eventNode);
  }

  _event_header = new SQEvent_v1();
  PHIODataNode<PHObject>* eventHeaderNode = new PHIODataNode<PHObject>(_event_header, "SQEvent", "PHObject");
  eventNode->addNode(eventHeaderNode);
  if(verbosity >= Fun4AllBase::VERBOSITY_SOME) std::cout << Name() << ": DST/SQEvent Added" << std::endl;

  //look for output nodes
  _hit_vector = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if(!_hit_vector) 
  {
    std::cerr << Name() << ": hit vector not found. Abort." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SQTriggerAnalyzer::process_event(PHCompositeNode* topNode) 
{
  if(Verbosity() >= Fun4AllBase::VERBOSITY_EVEN_MORE)
  {
    std::cout << Name() << ": process_event" << std::endl;
  }
  _event_header->Reset();

  //Reset all primitive counters
  for(auto iter = trPrimByName.begin(); iter != trPrimByName.end(); ++iter) iter->second->counter = 0;

  //process the data
  bool allLvlPresent = buildHitPattern();

  //Search for patterns maching Matrix trigger
  if(MatrixAvailable && allLvlPresent)
  {
    for(int i = 0; i < 2; ++i)
    {
      path.clear();
      roads_found[i].clear();
      searchMatrix(matrix[i], 0, i);

      for(auto iter = roads_found[i].begin(); iter != roads_found[i].end(); ++iter)
      {
        ++(trPrimByName[iter->getPrimitiveName().Data()]->counter);
      }
    }
  }

  //Update event header with trigger info
  for(auto iter = trigConds.begin(); iter != trigConds.end(); ++iter)
  {
    if(iter->second.eval()) _event_header->set_trigger(iter->first, true);
  }

  //Update the hit trigger mask flag for future use
  updateHitFlags();
  if(Verbosity() >= Fun4AllBase::VERBOSITY_EVEN_MORE)
  {
    _event_header->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SQTriggerAnalyzer::End(PHCompositeNode* topNode) 
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void SQTriggerAnalyzer::registerTrigger(const SQEvent::TriggerMask bit, const std::string& expr)
{
  if(trigConds.find(bit) != trigConds.end())
  {
    std::cout << Name() << ": requested trigger bit " << bit << " is occupied, will do nothing" << std::endl;
    return;
  }

  TriggerCondition trigcond;
  trigcond.trigBit = bit;
  trigcond.trigExpr.Compile(expr.c_str());

  for(int i = 0; i < trigcond.trigExpr.GetNpar(); ++i)
  {
    std::string primName = trigcond.trigExpr.GetParName(i);
    if(trPrimByName.find(primName) != trPrimByName.end()) 
    {
      trigcond.trigPrims.push_back(trPrimByName[primName]);
      continue;
    }

    TriggerPrimitive* tp = new TriggerPrimitive(primName);
    trPrimByName[primName] = tp;
    trigcond.trigPrims.push_back(tp);
  }

  trigConds[bit] = trigcond;
}

bool SQTriggerAnalyzer::loadMatrixFile()
{
  using namespace std;

  string line;
  ifstream fin(matrixFileName.c_str(), ifstream::in);
  if(fin.is_open()) 
  {
    while(getline(fin, line)) 
    {
      stringstream ss(line);

      int charge, roadID;
      int uIDs[nTrHits];
      double pXmin, sigWeight, bkgRate;

      ss >> charge >> roadID;
      for(int i = 0; i < nTrHits; ++i) ss >> uIDs[i];
      ss >> pXmin >> sigWeight >> bkgRate;

      SQTriggerRoad road(nTrHits);
      road.setCharge(charge);
      road.setRoadID(roadID);
      road.setSigWeight(sigWeight);
      road.setBkgRate(bkgRate);
      road.setPxMin(pXmin);
      road.setTrIDs(uIDs);

      roads[(-charge + 1)/2].insert(std::map<TString, SQTriggerRoad>::value_type(road.getStringID(), road));
    }
    return true;
  }

  return false;
}

void SQTriggerAnalyzer::buildMatrix() 
{
  for(int i = 0; i < 2; ++i) 
  {
    matrix[i] = new MatrixNode(-1);
    for(auto iter = roads[i].begin(); iter != roads[i].end(); ++iter) 
    {
      // build matrix
      MatrixNode* parentNode[nTrHits + 1]; //NOTE: the last entry is useless, just to keep the following code simpler
      parentNode[0] = matrix[i];
      for(int j = 0; j < nTrHits; ++j) 
      {
        int uniqueID = iter->second.getTrID(j);
        if(parentNode[j]->children.find(uniqueID) != parentNode[j]->children.end()) //existing node
        {
          parentNode[j+1] = parentNode[j]->children[uniqueID];
        }
        else //new node
        {
          MatrixNode* newNode = new MatrixNode(uniqueID);
          parentNode[j]->add(newNode);
          parentNode[j+1] = newNode;
        }
      }

      //build lookup tables
      for(int j = 0; j < nTrHits; ++j)
      {
        int detectorID = iter->second.getTrDetectorID(j);
        if(matrixLevelByID.find(detectorID) != matrixLevelByID.end()) continue;
        matrixLevelByID[detectorID] = j;
      }
    }
  }
}

bool SQTriggerAnalyzer::buildHitPattern()
{
  data.clear();
  
  //insert a dummy common root node first
  std::set<int> vertex; vertex.insert(-1);
  data.push_back(vertex);
  
  std::set<int> trHits[nTrHits];
  int nHitsAll = _hit_vector->size();
  for(int i = 0; i < nHitsAll; ++i)
  {
    SQHit* hit = _hit_vector->at(i);

    int detectorID = hit->get_detector_id();
    if(detectorID <= nChamberPlanes) continue;
    
    //Update detector-based trigger primitive counters - for NIM trigger
    int elementID = hit->get_element_id();
    std::string detectorName = p_geomSvc->getDetectorName(detectorID);
    if(trPrimByName.find(detectorName) != trPrimByName.end())
    {
      if(elementID >= trPrimByName[detectorName]->lo && elementID < trPrimByName[detectorName]->hi) ++(trPrimByName[detectorName]->counter);
    }

    //Insert the hit uID for matrix search if applicable - for Matrix trigger
    if(MatrixAvailable && matrixLevelByID.find(detectorID) != matrixLevelByID.end())
    {
      trHits[matrixLevelByID[detectorID]].insert(1000*detectorID + elementID);
    }
  }
  
  for(int i = 0; i < nTrHits; ++i) 
  {
    if(trHits[i].empty()) return false;
    data.push_back(trHits[i]);
  }
  
  return true;
}

int SQTriggerAnalyzer::updateHitFlags()
{
  std::set<int> uIDs_roads;
  for(int i = 0; i < 2; ++i)
  {
    for(auto iter = roads_found[i].begin(); iter != roads_found[i].end(); ++iter)
    {
      for(unsigned int j = 0; j < nTrHits; ++j) uIDs_roads.insert(iter->getTrID(j));
    }
  }

  int nHitsUpdated = 0;
  int nHitsAll = _hit_vector->size();
  for(int i = 0; i < nHitsAll; ++i)
  {
    SQHit* hit = _hit_vector->at(i);
    
    int detectorID = hit->get_detector_id();
    if(detectorID <= nChamberPlanes) continue;

    int uID = 1000*detectorID + hit->get_element_id();
    if(uIDs_roads.find(uID) != uIDs_roads.end()) 
    {
      hit->set_trigger_mask(true);
      ++nHitsUpdated;
    }
  }

  return nHitsUpdated;
}

void SQTriggerAnalyzer::searchMatrix(MatrixNode* node, int level, int index) 
{
  path.push_back(node->uniqueID);
  if(node->children.empty())   //we are at the end of the tree
  {
    SQTriggerRoad road_found(path);
    roads_found[index].push_back(road_found);

    path.pop_back();
    return;
  }

  for(auto iter = node->children.begin(); iter != node->children.end(); ++iter) 
  {
    if(data[level+1].find(iter->first) == data[level+1].end()) continue;
    searchMatrix(iter->second, level+1, index);
  }
  path.pop_back();
}

void SQTriggerAnalyzer::deleteMatrix(MatrixNode* node) 
{
  if(node == nullptr) return;
  if(node->children.empty()) 
  {
    delete node;
    return;
  }

  for(auto iter = node->children.begin(); iter != node->children.end(); ++iter) 
  {
    deleteMatrix(iter->second);
  }
  delete node;
}

void SQTriggerAnalyzer::printHitPattern() 
{
  for(unsigned int i = 1; i < nTrHits; ++i) 
  {
    std::cout << "Lv. " << i << ":  ";
    for(auto iter = data[i].begin(); iter != data[i].end(); ++iter) std::cout << *iter << "  ";
    std::cout << std::endl;
  }
}

void SQTriggerAnalyzer::printPath() 
{
  std::cout << "Found one road: " << std::endl;
  for(auto iter = path.begin(); iter != path.end(); ++iter) std::cout << *iter << " === ";
  std::cout << std::endl;
}

void SQTriggerAnalyzer::initDefault()
{
  initDefaultMatrix();
  initDefaultNIM();
}

void SQTriggerAnalyzer::initDefaultMatrix()
{
  registerTrigger(SQEvent::MATRIX1, "([TP] && [BM]) || ([TM] && [BP])");
  registerTrigger(SQEvent::MATRIX2, "([TP] && [TM]) || ([BP] && [BM])");
  registerTrigger(SQEvent::MATRIX3, "([TP] && [BP]) || ([TM] && [BM])");
  registerTrigger(SQEvent::MATRIX4, "[TP] || [TM] || [BP] || [BM]");
}

void SQTriggerAnalyzer::initDefaultNIM()
{
  registerTrigger(SQEvent::NIM1, "([H1T-4-19] && [H2T] && [H3T] && [H4T]) || ([H1B-4-19] && [H2B] && [H3B] && [H4B])");
  registerTrigger(SQEvent::NIM2, "([H1L] && [H2L] && [H4Y1L] && [H4Y2L]) || ([H1R] && [H2R] && [H4Y1R] && [H4Y1R])");
}

void SQTriggerAnalyzer::print(std::ostream& os) const
{
  using namespace std;

  os << "------------------------------------------------------------------------------" << endl;
  os << "Configuration of trigger analyzer " << Name() << endl;
  os << "------------------------------------------------------------------------------" << endl;

  if(MatrixAvailable)
  {
    os << "Matrix lookup table initialized from " << matrixFileName << endl;
    os << roads[0].size() << " positive roads, " << roads[1].size() << " negative roads." << endl;
    os << "Trigger detector level max = " << nTrHits << endl;
    for(auto iter = matrixLevelByID.begin(); iter != matrixLevelByID.end(); ++iter)
    {
      os << iter->first << " " << p_geomSvc->getDetectorName(iter->first) << " -> " << iter->second << ", ";
    }
    os << endl;
  }

  os << "Following trigger primitives are enabled " << endl;
  for(auto iter = trPrimByName.begin(); iter != trPrimByName.end(); ++iter)
  {
    os << "  " << iter->first << ": " << iter->second->name << " " << iter->second->lo << " " << iter->second->hi << endl;
  }

  os << "Following trigger conditions are enabled " << endl;
  for(auto iter = trigConds.begin(); iter != trigConds.end(); ++iter)
  {
    os << "  " << iter->second.trigBit << ": " << iter->second.trigExpr.GetExpFormula() << endl;
    for(unsigned int i = 0; i < iter->second.trigPrims.size(); ++i)
    {
      os << "    " << i << ": " << iter->second.trigPrims[i]->name << " " << iter->second.trigPrims[i]->lo << " " << iter->second.trigPrims[i]->hi << endl;
    }
  }

  os << "------------------------------------------------------------------------------" << endl;
}

SQTriggerAnalyzer::MatrixNode::MatrixNode(int uID): uniqueID(uID)
{
  children.clear();
}

void SQTriggerAnalyzer::MatrixNode::add(MatrixNode* child)
{
  children[child->uniqueID] = child;
}

SQTriggerAnalyzer::TriggerPrimitive::TriggerPrimitive(TString rawName): lo(-1), hi(99999), counter(0)
{
  TObjArray* arr = rawName.Tokenize("-");
  arr->SetOwner(kTRUE);

  name = ((TObjString*)arr->At(0))->String().Data();
  if(arr->GetEntries() == 3)
  {
    lo = ((TObjString*)arr->At(1))->String().Atoi();
    hi = ((TObjString*)arr->At(2))->String().Atoi();
  }

  delete arr;
}

bool SQTriggerAnalyzer::TriggerCondition::eval()
{
  for(auto iter = trigPrims.begin(); iter != trigPrims.end(); ++iter)
  {
    trigExpr.SetParameter((*iter)->name.c_str(), (*iter)->counter);
  }
  return trigExpr.EvalPar(nullptr) > 0.5;
}

void SQTriggerAnalyzer::TriggerCondition::reset()
{
  for(auto iter = trigPrims.begin(); iter != trigPrims.end(); ++iter)
  {
    (*iter)->counter = 0;
  }
}
