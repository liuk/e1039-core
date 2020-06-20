#include "SQCosmicGen.h"

#include "PHG4Particlev2.h"
#include "PHG4InEvent.h"
#include "PHG4VtxPoint.h"
#include "PHG4TruthInfoContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHRandomSeed.h>

#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleDefinition.hh>

#include <TMath.h>

#include <cstdlib>
#include <cmath>
#include <cassert>

SQCosmicGen::SQCosmicGen(const std::string& name): 
  PHG4ParticleGeneratorBase(name),
  _ineve(nullptr),
  _mom_min(1.),
  _mom_max(100.),
  _theta_min(-0.5*TMath::Pi()),
  _theta_max(0.5*TMath::Pi()),
  _prob_mup(1.3/(1. + 1.3)),
  _prob_mum(1./(1. + 1.3)),
  _p0(55.6),
  _p1(1.04),
  _p2(64.0),
  _altitude(400.),
  _size_x(800.),
  _size_z(4000.),
  _prob_max(1.),
  _rndm(PHRandomSeed())
{}

int SQCosmicGen::InitRun(PHCompositeNode* topNode) 
{
  _ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if(!_ineve) 
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode* dstNode;
    dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
      
    _ineve = new PHG4InEvent();
    PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(_ineve, "PHG4INEVENT", "PHObject");
    dstNode->addNode(newNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void SQCosmicGen::set_mom_range(const double lo, const double hi)
{
  _mom_min = lo;
  _mom_max = hi;

  if(_mom_min < 1.)
  {
    if(Verbosity() > 2) std::cout << Name() << ": mom_min should be larger than 1 GeV" << std::endl;
    _mom_min = 1.;
  }

  if(_theta_max > 1000.)
  {
    if(Verbosity() > 2) std::cout << Name() << ": mom_max should be smaller than 1 TeV" << std::endl;
    _mom_max = 1000.;
  }

  _prob_max = cosmicProb(_mom_min, 0.);
}

void SQCosmicGen::set_theta_range(const double lo, const double hi)
{
  _theta_min = lo;
  _theta_max = hi;

  if(_theta_min < -0.5*TMath::Pi())
  {
    if(Verbosity() > 2) std::cout << Name() << ": theta_min should be larger than -pi/2." << std::endl;
    _theta_min = -0.5*TMath::Pi();
  }

  if(_theta_max > 0.5*TMath::Pi())
  {
    if(Verbosity() > 2) std::cout << Name() << ": theta_max should be smaller than pi/2." << std::endl;
    _theta_max = 0.5*TMath::Pi();
  }
}

void SQCosmicGen::set_charge_ratio(const double p, const double n)
{
  _prob_mup = p/(p + n);
  _prob_mum = n/(p + n);
}

int SQCosmicGen::process_event(PHCompositeNode* topNode) 
{
  //Generate the charge randomly according to the charge ratio of 1.25
  int pdgcode = _rndm.Rndm() < _prob_mup ? -13 : 13;
  std::string pdgname = get_pdgname(pdgcode);

  //Generate a random vertex on the virtual generation plane
  double x_vtx = 0.5*uniformRand(-1., 1.)*_size_x;
  double z_vtx = 0.5*uniformRand(-1., 1.)*_size_z;
  double y_vtx = _altitude;
  int vtxID = _ineve->AddVtx(x_vtx, y_vtx, z_vtx, 0.);

  //Generate a muon and set common features
  PHG4Particle* particle = new PHG4Particlev2();
  particle->set_track_id(0);
  particle->set_vtx_id(vtxID);
  particle->set_parent_id(0);
  particle->set_pid(pdgcode);
  particle->set_name(pdgname);

  //Generate 3-momentum based on the probability function
  bool cosmic = true;
  double p, theta;
  while(cosmic)
  {
    p = uniformRand(_mom_min, _mom_max);
    theta = uniformRand(_theta_min, _theta_max);

    cosmic = _rndm.Rndm() > (cosmicProb(p, theta)/_prob_max);
  }

  double phi = _rndm.Rndm()*2.*TMath::Pi();
  double py = -p*TMath::Cos(theta);
  double px = p*TMath::Sin(theta)*TMath::Cos(phi);
  double pz = p*TMath::Sin(theta)*TMath::Sin(phi);
  double m  = 0.105658;
  double e  = TMath::Sqrt(p*p + m*m);
  particle->set_px(px);
  particle->set_py(py);
  particle->set_pz(pz);
  particle->set_e(e);
  _ineve->AddParticle(vtxID, particle);  //ownership of particle is transferred to _ineve and released in its Reset action

  if(Verbosity() > Fun4AllBase::VERBOSITY_A_LOT)
  {
    _ineve->identify();
    std::cout << "-------------------------------------------------------------------" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

double SQCosmicGen::cosmicProb(double p, double theta)
{
  return TMath::Power(_p0, TMath::Cos(_p1*theta))/p/p/_p2;
}

double SQCosmicGen::uniformRand(const double lo, const double hi)
{
  return lo + (hi - lo)*_rndm.Rndm();
}
