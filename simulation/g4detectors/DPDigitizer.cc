#include "DPDigitizer.h"

#include <interface_main/SQMCHit_v1.h>
#include <interface_main/SQHitVector_v1.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <geom_svc/GeomSvc.h>

#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>

#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/Randomize.hh>

#include <TMath.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include <boost/lexical_cast.hpp>

#define LogDebug(exp)       std::cout<<"DEBUG: "  <<__FUNCTION__<<": "<<__LINE__<<": "<< exp << std::endl

using namespace std;

void DPDigiPlane::preCalculation()
{
    sinth = TMath::Sin(angleFromVert + rZ);
    costh = TMath::Cos(angleFromVert + rZ);
    tanth = TMath::Tan(angleFromVert + rZ);

    wc = getW(xc, yc);

    rotM[0][0] = TMath::Cos(rZ)*TMath::Cos(rY);
    rotM[0][1] = TMath::Cos(rZ)*TMath::Sin(rX)*TMath::Sin(rY) - TMath::Cos(rX)*TMath::Sin(rZ);
    rotM[0][2] = TMath::Cos(rX)*TMath::Cos(rZ)*TMath::Sin(rY) + TMath::Sin(rX)*TMath::Sin(rZ);
    rotM[1][0] = TMath::Sin(rZ)*TMath::Cos(rY);
    rotM[1][1] = TMath::Sin(rZ)*TMath::Sin(rX)*TMath::Sin(rY) + TMath::Cos(rZ)*TMath::Cos(rX);
    rotM[1][2] = TMath::Sin(rZ)*TMath::Sin(rY)*TMath::Cos(rX) - TMath::Cos(rZ)*TMath::Sin(rX);
    rotM[2][0] = -TMath::Sin(rY);
    rotM[2][1] = TMath::Cos(rY)*TMath::Sin(rX);
    rotM[2][2] = TMath::Cos(rY)*TMath::Cos(rX);

    uVec[0] = TMath::Cos(angleFromVert);
    uVec[1] = TMath::Sin(angleFromVert);
    uVec[2] = 0.;

    vVec[0] = -TMath::Sin(angleFromVert);
    vVec[1] = TMath::Cos(angleFromVert);
    vVec[2] = 0.;

    nVec[0] = 0.;
    nVec[1] = 0.;
    nVec[2] = 1.;

    xVec[0] = 1.;
    xVec[1] = 0.;
    xVec[2] = 0.;

    yVec[0] = 0.;
    yVec[1] = 1.;
    yVec[2] = 0.;

    //rotate u/v vector by the rotation matrix
    double temp[3];
    for(int i = 0; i < 3; ++i) temp[i] = uVec[i];
    for(int i = 0; i < 3; ++i)
    {
        uVec[i] = 0.;
        for(int j = 0; j < 3; ++j) uVec[i] += rotM[i][j]*temp[j];
    }

    for(int i = 0; i < 3; ++i) temp[i] = vVec[i];
    for(int i = 0; i < 3; ++i)
    {
        vVec[i] = 0.;
        for(int j = 0; j < 3; ++j) vVec[i] += rotM[i][j]*temp[j];
    }

    for(int i = 0; i < 3; ++i) temp[i] = xVec[i];
    for(int i = 0; i < 3; ++i)
    {
        xVec[i] = 0.;
        for(int j = 0; j < 3; ++j) xVec[i] += rotM[i][j]*temp[j];
    }

    for(int i = 0; i < 3; ++i) temp[i] = yVec[i];
    for(int i = 0; i < 3; ++i)
    {
        yVec[i] = 0.;
        for(int j = 0; j < 3; ++j) yVec[i] += rotM[i][j]*temp[j];
    }

    //n vector is the cross product of u and v
    nVec[0] = uVec[1]*vVec[2] - vVec[1]*uVec[2];
    nVec[1] = uVec[2]*vVec[0] - vVec[2]*uVec[0];
    nVec[2] = uVec[0]*vVec[1] - vVec[0]*uVec[1];
}

bool DPDigiPlane::intercept(double tx, double ty, double x0, double y0, G4ThreeVector& pos, double& w)
{
    //Refer to http://geomalgorithms.com/a05-_intersect-1.html
    //double u[3] = {tx, ty, 1};
    //double p0[3] = {x0, y0, 0};
    //double v0[3] = {xc, yc, zc};
    //double n[3] = nVec[3];
    //double w[3] = p0[3] - v0[3];

    double det = tx*nVec[0] + ty*nVec[1] + nVec[2];
    double dpos[3] = {x0 - xc, y0 - yc, -zc};
    double si = -(nVec[0]*dpos[0] + nVec[1]*dpos[1] + nVec[2]*dpos[2])/det;

    pos[0] = x0 + tx*si;
    pos[1] = y0 + ty*si;
    pos[2] = si;

    //original
    //double vcp[3] = {vVec[1] - vVec[2]*ty, vVec[2]*tx - vVec[0], vVec[0]*ty - vVec[1]*tx};
    //w = (vcp[0]*dpos[0] + vcp[1]*dpos[1] + vcp[2]*dpos[2])/det;

    //yuhw
    w = uVec[0]*(pos[0]-xc) + uVec[1]*(pos[1]-yc) + uVec[2]*(pos[2]-zc);

    return isInPlane(pos[0], pos[1], pos[2]);
}

std::ostream& operator << (std::ostream& os, const DPDigiPlane& plane)
{
    os << "DigiPlane ID = " << plane.detectorID << ", name = " << plane.detectorName << " belongs to group " << plane.detectorGroupName << "\n"
       << "nElements = " << plane.nElements
			 << ", center x = " << plane.xc << ", y = " << plane.yc << ", z = " << plane.zc << "\n"
			 << ", nVec: {" << plane.nVec[0] << ", " << plane.nVec[1] << ", " << plane.nVec[2] << "} "
			 << ", uVec: {" << plane.uVec[0] << ", " << plane.uVec[1] << ", " << plane.uVec[2] << "} "
			 << ", vVec: {" << plane.vVec[0] << ", " << plane.vVec[1] << ", " << plane.vVec[2] << "} "
			 ;
    return os;
}

DPDigitizer::DPDigitizer(const std::string &name) :
		p_geomSvc(nullptr)
{
}

DPDigitizer::~DPDigitizer() {
}

void DPDigitizer::digitize(std::string detectorGroupName, PHG4Hit& g4hit)
{
    if(Verbosity() > 2){
      LogDebug("DPDigitizer::digitize: " << map_groupID[detectorGroupName].size());
    }

    int track_id = g4hit.get_trkid();

    // calculate the central position in each detector group, then linearly extrapolate the hits
    // to each individual plane, this is assuming there is no magnetic field in the detector, or
    // the bending is negligible
    double tx = g4hit.get_px(0)/g4hit.get_pz(0);
    double ty = g4hit.get_py(0)/g4hit.get_pz(0);
    double x0 = (g4hit.get_x(0) - tx*g4hit.get_z(0));///cm;
    double y0 = (g4hit.get_y(0) - ty*g4hit.get_z(0));///cm;

    //temporary variabels
    double w;
    G4ThreeVector pos;
    for(std::vector<int>::iterator dpid = map_groupID[detectorGroupName].begin();
    		dpid != map_groupID[detectorGroupName].end();
    		++dpid)
    {
				if(Verbosity() > 2) {
					cout << "DEBUG: detectorGroupName: " << detectorGroupName << endl;
					cout << "DEBUG: detectorName: " << digiPlanes[*dpid].detectorName << endl;
				}

        //check if the track intercepts the plane
        if(!digiPlanes[*dpid].intercept(tx, ty, x0, y0, pos, w)) continue;

        int DP_elementID = TMath::Nint((digiPlanes[*dpid].nElements + 1.0)/2.0 +
        		//(w - digiPlanes[*dpid].xPrimeOffset - digiPlanes[*dpid].xc*digiPlanes[*dpid].costh + digiPlanes[*dpid].yc*digiPlanes[*dpid].sinth)/digiPlanes[*dpid].spacing) ;
        		(w - digiPlanes[*dpid].xPrimeOffset)/digiPlanes[*dpid].spacing) ;
        double driftDistance = w - digiPlanes[*dpid].spacing*(DP_elementID - digiPlanes[*dpid].nElements/2. - 0.5) - digiPlanes[*dpid].xPrimeOffset;
        if(DP_elementID < 1 || DP_elementID > digiPlanes[*dpid].nElements || fabs(driftDistance) > 0.5*digiPlanes[*dpid].cellWidth) continue;

        SQMCHit_v1 *digiHit = new SQMCHit_v1();

        digiHit->set_track_id(track_id);
        digiHit->set_g4hit_id(g4hit.get_hit_id());

        digiHit->set_truth_x(pos[0]);
        digiHit->set_truth_y(pos[1]);
        digiHit->set_truth_z(pos[2]);

        //digiHit.fPDGCode = vHit.particlePDG;
        //digiHit->set_detector_id(digiPlanes[*dpid].detectorID);
        digiHit->set_element_id(DP_elementID);
        digiHit->set_drift_distance(driftDistance);
        //digiHit.fMomentum.SetXYZ(vHit.get_px(0)/GeV, vHit.get_py(0)/GeV, vHit.get_pz(0)/GeV);
        //digiHit.fPosition.SetXYZ(pos[0], pos[1], pos[2]);
        //digiHit.fDepEnergy = vHit.edep/GeV;

        //if(realize(digiHit)) g4hit.digiHits.push_back(digiHit);

        //see if it also hits the next elements in the overlap region
        if(fabs(driftDistance) > 0.5*digiPlanes[*dpid].cellWidth - digiPlanes[*dpid].overlap)
        {
            if(driftDistance > 0. && DP_elementID != digiPlanes[*dpid].nElements)
            {
                digiHit->set_element_id(DP_elementID + 1);
                digiHit->set_drift_distance(driftDistance - digiPlanes[*dpid].spacing);
                //if(realize(digiHit)) g4hit.digiHits.push_back(digiHit);
            }
            else if(driftDistance < 0. && DP_elementID != 1)
            {
                digiHit->set_element_id(DP_elementID - 1);
                digiHit->set_drift_distance(driftDistance + digiPlanes[*dpid].spacing);
                //if(realize(digiHit)) g4hit.digiHits.push_back(digiHit);
            }
        }

        string detName = digiPlanes[*dpid].detectorName;

        //p_geomSvc->toLocalDetectorName(detName, DP_elementID);
        digiHit->set_detector_id(p_geomSvc->getDetectorID(detName));

        //TODO temp solution
        if(digiHit->get_detector_id()==0) continue;

        digiHit->set_pos(p_geomSvc->getMeasurement(digiHit->get_detector_id(), digiHit->get_element_id()));
        //digiHit->set_pos(w);

        // FIXME figure this out
        digiHit->set_in_time(1);
        digiHit->set_hodo_mask(1);

        digiHit->set_hit_id(digits->size());

        if(Verbosity() > 2) {
        	cout << digiPlanes[*dpid] << endl;
        	cout << "DEBUG: DigiHit: DPSim: ID: " << *dpid << ", Name: " << digiPlanes[*dpid].detectorName << endl;
        	cout << "DEBUG: DigiHit: GeoSvc: ID: " << digiHit->get_detector_id() << ", Name: " << detName << ", "<< endl;
        	cout
					<< "DEBUG: DigiHit: hit_id: " << digiHit->get_hit_id()
					<< ", w: " << w
					<< ", pos: " << digiHit->get_pos()
					<< endl;
        	//digiHit->identify();
        }

        digits->push_back(digiHit);
    }

    //split the energy deposition to all digihits
//    for(std::vector<SQHit>::iterator iter = g4hit.digiHits.begin(); iter != g4hit.digiHits.end(); ++iter)
//    {
//        iter->fDepEnergy = iter->fDepEnergy/g4hit.digiHits.size();
//    }
}

namespace {
	string toGroupName(string in) {

		std::vector<std::regex> regs;
		std::vector<std::string> reps;

		regs.push_back(std::regex("(H)([1-3])([T,B])$"));      //HX
		reps.push_back("$1$2X");

		regs.push_back(std::regex("(H)([1-3])([L,R])$"));      //HY
		reps.push_back("$1$2Y");

		regs.push_back(std::regex("(H4)(T|B)(.*)"));      //H4X
		reps.push_back("$1X");

		regs.push_back(std::regex("(H4Y)([1-2])(.*)"));      //H4Y
		reps.push_back("$1Y$2");

		regs.push_back(std::regex("(P)([0-9])(H|V)(.*)$"));         //photo-tube
		reps.push_back("$1$2$3");

		regs.push_back(std::regex("(D)(.*)(U|X|V|Up|Xp|Vp)$"));         //Drift chamber
		reps.push_back("$1$2");

		regs.push_back(std::regex("(DP)(.*)([L,R])$"));//
		reps.push_back("$1$2");

		for(unsigned int i=0; i<regs.size(); ++i) {
			if(std::regex_match(in, regs[i]))
				return std::regex_replace(in, regs[i], reps[i]);
		}

		return "";
	}
}

int DPDigitizer::InitRun(PHCompositeNode* topNode) {
  if(Verbosity() > 2){
    LogDebug("DPDigitizer::InitRun");
  }

//  const char* mysqlServer = "e906-db1.fnal.gov";
//  const int mysqlPort = 3306;
//  const char* geometrySchema = "user_liuk_geometry_DPTrigger";

  const char* mysqlServer = "seaquestdb01.fnal.gov";
  const int mysqlPort = 3310;
  const char* geometrySchema = "geometry_G17_run3";
  const char* login = "seaguest";
  const char* password = "qqbar2mu+mu-";

  const char* detectorEffResol = "";

  //Load basic setup
  char query[500];

  //user_liuk_geometry_DPTrigger
//  sprintf(query, "SELECT Planes.detectorName,spacing,cellWidth,overlap,numElements,angleFromVert,"
//                 "xPrimeOffset,x0+deltaX,y0+deltaY,z0+deltaZ,planeWidth,planeHeight,theta_x+rotX,"
//                 "theta_y+rotY,theta_z+rotZ,Planes.detectorID,Planes.detectorGroupName,triggerLv "
//                 "FROM %s.Planes,%s.Alignments WHERE Planes.detectorName=Alignments.detectorName",
//								 geometrySchema, geometrySchema);

  //geometry_G17_run3
  sprintf(query, "SELECT Planes.detectorName,spacing,cellWidth,overlap,numElements,angleFromVert,"
                 "xPrimeOffset,x0,y0,z0,planeWidth,planeHeight,theta_x,"
                 "theta_y,theta_z,Planes.detectorID "
                 "FROM %s.Planes",
								 geometrySchema);

  TSQLServer* server = TSQLServer::Connect(Form("mysql://%s:%d", mysqlServer, mysqlPort), login, password);
  TSQLResult* res = server->Query(query);

	if(Verbosity() > 2) {
		LogInfo(query);
	}

  map_groupID.clear();
  map_detectorID.clear();

  unsigned int nRows = res->GetRowCount();
  for(unsigned int i = 0; i < nRows; ++i)
  {
      TSQLRow* row = res->Next();
      int index = boost::lexical_cast<int>(row->GetField(15));

      //assert(index <= NDETPLANES && "detectorID from database exceeds upper limit!");
      if(index > NDETPLANES) {
        continue;
      }

      digiPlanes[index].detectorID = index;
      digiPlanes[index].detectorName = row->GetField(0);
      digiPlanes[index].spacing = boost::lexical_cast<double>(row->GetField(1));
      digiPlanes[index].cellWidth = boost::lexical_cast<double>(row->GetField(2));
      digiPlanes[index].overlap = boost::lexical_cast<double>(row->GetField(3));
      digiPlanes[index].nElements = boost::lexical_cast<int>(row->GetField(4));
      digiPlanes[index].angleFromVert = boost::lexical_cast<double>(row->GetField(5));
      digiPlanes[index].xPrimeOffset = boost::lexical_cast<double>(row->GetField(6));
      digiPlanes[index].xc = boost::lexical_cast<double>(row->GetField(7));
      digiPlanes[index].yc = boost::lexical_cast<double>(row->GetField(8));
      digiPlanes[index].zc = boost::lexical_cast<double>(row->GetField(9));
      digiPlanes[index].planeWidth = boost::lexical_cast<double>(row->GetField(10));
      digiPlanes[index].planeHeight = boost::lexical_cast<double>(row->GetField(11));
      digiPlanes[index].rX = boost::lexical_cast<double>(row->GetField(12));
      digiPlanes[index].rY = boost::lexical_cast<double>(row->GetField(13));
      digiPlanes[index].rZ = boost::lexical_cast<double>(row->GetField(14));

      digiPlanes[index].detectorGroupName = toGroupName(digiPlanes[index].detectorName);

      //user_liuk_geometry_DPTrigger
      //digiPlanes[index].detectorGroupName = row->GetField(16);
      //digiPlanes[index].triggerLv = boost::lexical_cast<int>(row->GetField(17));


      std::regex eP1("(P)([1-2])(H|V)([2-9])(b|f)$");
      if(std::regex_match(digiPlanes[index].detectorName, eP1)) {
      	LogInfo(digiPlanes[index].detectorName);
      	continue;
      }

      std::regex eP2("(P)([1-2])(H|V)(1)(b|f)$");
      if(std::regex_match(digiPlanes[index].detectorName, eP2)) {
      	LogInfo(digiPlanes[index].detectorName);
      	string temp = digiPlanes[index].detectorName;
      	int dummy;
      	p_geomSvc->toLocalDetectorName(temp, dummy);
      	digiPlanes[index].detectorName = temp;
      	LogInfo(digiPlanes[index].detectorName);
        digiPlanes[index].nElements = 72;
        digiPlanes[index].planeWidth = 365.76;
        digiPlanes[index].planeHeight = 365.76;
        // TODO hard coding for now
      	digiPlanes[index].xc = 0;
        digiPlanes[index].yc = 0;
      }

      digiPlanes[index].preCalculation();

      //TODO: implement RT/eff/resolution here
      digiPlanes[index].efficiency.resize(digiPlanes[index].nElements+1, 1.);
      digiPlanes[index].resolution.resize(digiPlanes[index].nElements+1, 0.);

      map_groupID[digiPlanes[index].detectorGroupName].push_back(index);
      map_detectorID[digiPlanes[index].detectorName] = digiPlanes[index].detectorID;

      if(Verbosity() > 2) {
      	cout
				<< "index: " << index << ", "
				<< digiPlanes[index].detectorGroupName << ", "
				<< digiPlanes[index].detectorID << ", "
				<< digiPlanes[index].detectorName << ", "
				<< digiPlanes[index].nElements << " "
				<< endl;
      }

      delete row;
  }

  map_g4name_group["C1X"] = "D1";
//  map_g4name_group["C1V"] = "D1";
//  map_g4name_group["C1U"] = "D1";
//  map_g4name_group["C2U"] = "D2";
  map_g4name_group["C2X"] = "D2";
//  map_g4name_group["C2V"] = "D2";
  map_g4name_group["C3T"] = "D3p";
  map_g4name_group["C3B"] = "D3m";

  map_g4name_group["H1y"] = "H1Y";
  map_g4name_group["H1x"] = "H1X";
  map_g4name_group["H2y"] = "H2Y";
  map_g4name_group["H2x"] = "H2X";
  map_g4name_group["H3x"] = "H3X";

//  map_g4name_group["P1H"] = "P1Y";
//  map_g4name_group["P2V"] = "P1X";
//  map_g4name_group["P2H"] = "P2Y";
//  map_g4name_group["P1V"] = "P2X";

  map_g4name_group["P1H"] = "P1H";
  map_g4name_group["P2V"] = "P2V";
  map_g4name_group["P2H"] = "P2H";
  map_g4name_group["P1V"] = "P1V";

  map_g4name_group["H4y1L"] = "H4Y1";
  map_g4name_group["H4y1R"] = "H4Y1";
  map_g4name_group["H4y2L"] = "H4Y2";
  map_g4name_group["H4y2R"] = "H4Y2";
  map_g4name_group["H4xT"] = "H4X";
  map_g4name_group["H4xB"] = "H4X";


  delete res;
  delete server;

  //Load the detector realization setup
  std::ifstream fin(detectorEffResol);
  if(fin) {
    std::string line;
    while(getline(fin, line))
    {
        std::string detectorName;
        int elementID;
        double eff, res;

        std::stringstream ss(line);
        ss >> detectorName >> elementID >> eff >> res;

        if(eff >= 0. && eff <= 1. && res >= 0.)
        {
            digiPlanes[map_detectorID[detectorName]].efficiency[elementID] = eff;
            digiPlanes[map_detectorID[detectorName]].resolution[elementID] = res;
        }
    }
  };

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      cout << Name() << " DST Node missing, doing nothing." << std::endl;
      exit(1);
    }
  PHNodeIterator dstiter(dstNode);

  string digit_name = "SQHitVector";
  digits = findNode::getClass<SQHitVector>(topNode , digit_name);
	if (!digits){
		digits = new SQHitVector_v1();
		PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(digits, digit_name.c_str() , "PHObject");
		dstNode->addNode(newNode);
	}

	p_geomSvc = GeomSvc::instance();

  if(Verbosity() > 2) {
  	LogInfo(digiPlanes[41].detectorName);
  }

	return Fun4AllReturnCodes::EVENT_OK;
}

int DPDigitizer::process_event(PHCompositeNode* topNode) {
  if(Verbosity() > 2){
    LogDebug("DPDigitizer::process_event");
  }

  if(Verbosity() > 2) {
  	LogInfo(digiPlanes[41].detectorName);
  }

	for(auto detector_iter = map_g4name_group.begin();
			detector_iter != map_g4name_group.end(); ++detector_iter) {
		string g4name = detector_iter->first;
		string hitnodename = "G4HIT_" + g4name;

    if(Verbosity() > 2) {
      LogDebug(g4name);
      LogDebug(hitnodename);
    }

		PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
	  if (!hits)
	  {
	    cout << Name() << " Could not locate g4 hit node " << hitnodename << endl;
	    exit(1);
	  }

    if(Verbosity() > 2) {
      LogDebug(g4name);
    }

	  for(PHG4HitContainer::ConstIterator hit_iter = hits->getHits().first;
	  		hit_iter != hits->getHits().second; ++ hit_iter){
	  	PHG4Hit* hit = hit_iter->second;

	  	int track_id = hit->get_trkid();
	    //FIXME only keep primary hit only for now
	    if(track_id < 0) continue;

	    string group = map_g4name_group[g4name];
      if(Verbosity() > 2) {
        hit->identify();
        LogDebug(g4name);
        LogDebug(group);
      }
	  	try{
	  		digitize(group, *hit);
	  	}catch(...) {
	  		cout << Name() << " Failed digitize " << group << endl;
	  		return Fun4AllReturnCodes::ABORTEVENT;
	  	}
	  }
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

bool DPDigitizer::realize(SQHit& dHit)
{
    if(G4UniformRand() > digiPlanes[dHit.get_detector_id()].efficiency[dHit.get_element_id()]) return false;

    dHit.set_drift_distance(
    		dHit.get_drift_distance() + (G4RandGauss::shoot(0., digiPlanes[dHit.get_detector_id()].resolution[dHit.get_element_id()]))
				);
    return true;
}
