#include "GFTrack.h"

#include <iostream>
#include <algorithm>
#include <stdexcept>

#include <TVectorD.h>
#include <TMatrixDSym.h>

#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/MeasurementOnPlane.h>
#include <GenFit/WireMeasurement.h>
#include <GenFit/PlanarMeasurement.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/AbsFitterInfo.h>
#include <GenFit/KalmanFitterInfo.h>
#include <GenFit/StateOnPlane.h>
#include <GenFit/HMatrixPhi.h>

#include <phool/recoConsts.h>

#include "SRawEvent.h"

namespace
{
  //static flag to indicate the initialized has been done
  static bool inited = false;

  static double Z_TARGET;
  static double Z_DUMP;
  static double Z_UPSTREAM;

  static double X_BEAM;
  static double Y_BEAM;
  static double SIGX_BEAM;
  static double SIGY_BEAM;

  static GeomSvc* p_geomSvc;

  //initialize global variables
  void initGlobalVariables()
  {
    if(!inited) 
    {
      inited = true;
      p_geomSvc = GeomSvc::instance();

      recoConsts* rc = recoConsts::instance();
      Z_TARGET   = rc->get_DoubleFlag("Z_TARGET");
      Z_DUMP     = rc->get_DoubleFlag("Z_DUMP");
      Z_UPSTREAM = rc->get_DoubleFlag("Z_UPSTREAM");

      X_BEAM    = rc->get_DoubleFlag("X_BEAM");
      Y_BEAM    = rc->get_DoubleFlag("Y_BEAM");
      SIGX_BEAM = rc->get_DoubleFlag("SIGX_BEAM");
      SIGY_BEAM = rc->get_DoubleFlag("SIGY_BEAM");
    }
  }
};

namespace SQGenFit
{

GFTrack::GFTrack(): _track(nullptr), _trkrep(nullptr), _propState(nullptr), _virtMeas(nullptr), _trkcand(nullptr), _pdg(0)
{
  initGlobalVariables();
}

GFTrack::GFTrack(SRecTrack& recTrack):  _track(nullptr), _trkrep(nullptr), _propState(nullptr), _virtMeas(nullptr), _trkcand(nullptr), _pdg(0)
{
  _pdg = recTrack.getCharge() > 0 ? -13 : 13;
  _trkrep = new genfit::RKTrackRep(_pdg);
    
  TVector3 seed_mom = recTrack.getMomentumVecSt1();
  TVector3 seed_pos = recTrack.getPositionVecSt1();

  TVectorD seed_state(6);
  seed_state[0] = seed_pos.X();
  seed_state[1] = seed_pos.Y();
  seed_state[2] = seed_pos.Z();
  seed_state[3] = seed_mom.Px();
  seed_state[4] = seed_mom.Py();
  seed_state[5] = seed_mom.Pz();

  TMatrixDSym seed_cov(6);
  double uncertainty[6] = {10., 10., 10., 3., 3., 10.};
  for(int i = 0; i < 6; i++)
  {
    for(int j = 0; j < 6; j++)
    {
      seed_cov[i][j] = uncertainty[i]*uncertainty[j];
    }
  }
  _track = new genfit::Track(_trkrep, seed_state, seed_cov);

  _fitstates.clear();
  int nHits = recTrack.getNHits();
  for(int i = 0; i < nHits; ++i)
  {
    genfit::SharedPlanePtr detPlane(new genfit::DetPlane(recTrack.getGFPlaneO(i), recTrack.getGFPlaneU(i), recTrack.getGFPlaneV(i)));
    _fitstates.push_back(genfit::MeasuredStateOnPlane(recTrack.getGFState(i), recTrack.getGFCov(i), detPlane, _trkrep, recTrack.getGFAuxInfo(i)));
  }
}

GFTrack::GFTrack(Tracklet& tracklet): _track(nullptr), _trkrep(nullptr), _propState(nullptr), _virtMeas(nullptr), _trkcand(nullptr), _pdg(0)
{
  setTracklet(tracklet, 590., false);
}

GFTrack::~GFTrack()
{
  if(_track != nullptr) delete _track;
  if(_trkcand != nullptr) delete _trkcand;
}

GFTrack* GFTrack::clone(bool clean) const
{
  GFTrack* newTrack = new GFTrack();
  if(clean) 
  {
    newTrack->setTracklet(*_trkcand, 1900.);
  }
  else
  {
    newTrack->_track  = new genfit::Track(*_track);
    newTrack->_trkrep = newTrack->_track->getCardinalRep();
    newTrack->_pdg    = newTrack->_trkrep->getPDG();
    if(_trkcand != nullptr) newTrack->_trkcand = _trkcand->Clone();

    newTrack->restoreFromGenFitTrack();
  }

  return newTrack;
}

void GFTrack::restoreFromGenFitTrack(genfit::Track* trkin)
{
  genfit::Track* gftrk = trkin == nullptr ? _track : trkin;

  unsigned int nHits = gftrk->getNumPoints();
  for(unsigned int i = 0; i < nHits; ++i)
  {
    genfit::TrackPoint* tp = gftrk->getPoint(i);
    GFMeasurement* meas = const_cast<GFMeasurement*>(dynamic_cast<GFMeasurement*>(tp->getRawMeasurement()));
    meas->setTrackPtr(this);
    _measurements.insert(meas);
  }
}

void GFTrack::setVerbosity(unsigned int v)
{
  if(_trkrep != nullptr) _trkrep->setDebugLvl(v);
}

void GFTrack::addMeasurements(HitList_t& measurements)
{
  for(auto iter = measurements.begin(); iter != measurements.end(); ++iter)
  {
    if(!(*iter)->isEnabled()) continue;
    addMeasurement(*iter);
  }
}

void GFTrack::addMeasurement(GFMeasurement* measurement, const int id)
{
  measurement->setTrackPtr(this);

  genfit::TrackPoint* tp = new genfit::TrackPoint(measurement, _track);
  tp->setSortingParameter(measurement->getZ());
  _track->insertPoint(tp, id);
}

void GFTrack::addMeasurement(const SignedHit& hit, const int id)
{
  if(_trkcand != nullptr) _trkcand->hits.push_front(hit);

  GFMeasurement* meas = new GFMeasurement(hit);
  _measurements.insert(meas);
  addMeasurement(meas, id);
}

double GFTrack::testOneHitKalman(const SignedHit& hit)
{
  TVector3 ep1, ep2;
  p_geomSvc->getEndPoints(hit.hit.detectorID, hit.hit.elementID, ep1, ep2);

  double width = p_geomSvc->getCellWidth(hit.hit.detectorID);
  TVectorD meas(1); meas[0] = 0.;
  TMatrixDSym V(1); V[0][0] = width*width/12.;

  double chi2 = -1.;
  try
  {
    double len = extrapolateToLine(ep1, ep2);
    if(fabs(len) > 6000.) throw len;

    chi2 = updatePropState(meas, V);
  }
  catch(genfit::Exception& e)
  {
    return -1.;
  }
  catch(double len)
  {
    return -1.;
  }

  return chi2;
}

void GFTrack::getProjection(int detID, double& x, double& y, double& w, double& dw)
{
  x = 1.E6;
  y = 1.E6;
  w = 1.E6;
  dw = -1.E6;

  double z  = p_geomSvc->getPlanePosition(detID);
  double rZ = p_geomSvc->getStereoAngle(detID);

  TVector3 pO(0., 0., z);
  TVector3 pU(1., 0., 0.);
  TVector3 pV(0., 1., 0.);
  genfit::SharedPlanePtr destPlane(new genfit::DetPlane(pO, pU, pV));

  if(!setInitialStateForExtrap()) return;
  double len = -1.;
  try
  {
    genfit::AbsTrackRep* rep = _track->getCardinalRep();
    len = rep->extrapolateToPlane(*_propState, destPlane);
  } 
  catch(...)
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": projection to detectorID = " << detID << " failed." << std::endl;
    return;
  }

  TVector3 pos, mom;
  getExtrapPosMomCov(&pos, &mom, nullptr);

  x = pos.X();
  y = pos.Y();

  genfit::HMatrixPhi H(rZ);
  w = H.Hv(_propState->getState())[0];
  
  TMatrixDSym cov(_propState->getCov());
  H.HMHt(cov);
  dw = sqrt(cov(0, 0));
}

int GFTrack::getProjectedElementID(int detID, int fb)  
{
  //TODO: maybe able to optimize by avoiding the setInitialStateForExtrap by directly getting the fitter info
  // or add some cacheing machnism
  if(!setInitialStateForExtrap(fb)) return -1;

  TVector3 pos, mom;
  getExtrapPosMomCov(&pos, &mom);

  double tx = mom.Px()/mom.Pz();
  double ty = mom.Py()/mom.Pz();
  double w  = p_geomSvc->getInterception(detID, tx, ty, pos.X() - tx*pos.Z(), pos.Y() - ty*pos.Z());

  return p_geomSvc->getExpElementID(detID, w);
}

double GFTrack::getChi2() const
{
  genfit::AbsTrackRep* rep = _track->getCardinalRep();
  if(rep)
  {
    genfit::FitStatus* fs = _track->getFitStatus(rep);
    if(fs) return fs->getChi2();
  }
  
  return -1.;
}

double GFTrack::getNDF() const
{
  genfit::AbsTrackRep* rep = _track->getCardinalRep();
  if(rep)
  {
    genfit::FitStatus* fs = _track->getFitStatus(rep);
    if(fs) return fs->getNdf();
  }
  
  return -1.;
}

double GFTrack::getQuality() const
{
  genfit::AbsTrackRep* rep = _track->getCardinalRep();
  if(rep)
  {
    genfit::FitStatus* fs = _track->getFitStatus(rep);
    if(fs) return fs->getChi2()/fs->getNdf();
  }
  
  return -1.;
}

void GFTrack::getFittedPosMom(TVector3& pos, TVector3& mom) const
{
  _track->getFittedState().getPosMom(pos, mom);
}

double GFTrack::extrapolateToLine(TVector3& endPoint1, TVector3& endPoint2, const int startPtID)
{
  TVector3 linePoint = endPoint1;
  TVector3 lineDir   = endPoint2 - endPoint1;

  if(!setInitialStateForExtrap(startPtID)) return -9999.;

  genfit::AbsTrackRep* rep = _track->getCardinalRep();
  double len = rep->extrapolateToLine(*_propState, linePoint, lineDir);

  TVectorD hitcoord(7);
  hitcoord[0] = endPoint1.X();
  hitcoord[1] = endPoint1.Y();
  hitcoord[2] = endPoint1.Z();
  hitcoord[3] = endPoint2.X();
  hitcoord[4] = endPoint2.Y();
  hitcoord[5] = endPoint2.Z();
  hitcoord[6] = 0.;
  TMatrixDSym tempcov(7);
  tempcov.Zero();
  _virtMeas.reset(new genfit::WireMeasurement(hitcoord, tempcov, 999, 999, nullptr));

  return len;
}

double GFTrack::extrapolateToPlane(TVector3& pO, TVector3& pU, TVector3& pV, const int startPtID)
{
  genfit::SharedPlanePtr destPlane(new genfit::DetPlane(pO, pU, pV));

  if(!setInitialStateForExtrap(startPtID)) return -9999.;

  genfit::AbsTrackRep* rep = _track->getCardinalRep();
  double len = rep->extrapolateToPlane(*_propState, destPlane);

  TVectorD hitcoord(2);
  hitcoord[0] = pO.X();
  hitcoord[1] = pO.Y();
  TMatrixDSym tempcov(2);
  tempcov.UnitMatrix();

  genfit::PlanarMeasurement* pMeas = new genfit::PlanarMeasurement(hitcoord, tempcov, 998, 998, nullptr);
  pMeas->setPlane(destPlane);
  _virtMeas.reset(pMeas);

  return len;
}

double GFTrack::extrapolateToPoint(TVector3& point, const int startPtID)
{
  if(!setInitialStateForExtrap(startPtID)) return -9999.;

  genfit::AbsTrackRep* rep = _track->getCardinalRep();
  double len = rep->extrapolateToPoint(*_propState, point);

  //TODO: implement the virtual spacepoint measurement here

  return len;
}

double GFTrack::updatePropState(const TVectorD& meas, const TMatrixDSym& V)
{
  //get the H matrix from measurement
  std::unique_ptr<const genfit::AbsHMatrix> H(_virtMeas->constructHMatrix(_propState->getRep()));

  //Handcrafted KF
  TVectorD stateVec(_propState->getState());
  TMatrixDSym cov(_propState->getCov());
  double chi2 = 9999.;
  double ndf = meas.GetNrows();

  TVectorD res(meas - H->Hv(stateVec));
  TMatrixDSym covSumInv(cov);
  H->HMHt(covSumInv);
  covSumInv += V;
  genfit::tools::invertMatrix(covSumInv);

  TMatrixD CHt(H->MHt(cov));
  TVectorD update(TMatrixD(CHt, TMatrixD::kMult, covSumInv)*res);

  stateVec += update;
  covSumInv.Similarity(CHt);
  cov -= covSumInv;
  _propState->setStateCov(stateVec, cov);

  //calc chi2
  TVectorD resNew(meas - H->Hv(stateVec));
  TMatrixDSym HCHt(cov);
  H->HMHt(HCHt);
  HCHt -= V;
  HCHt *= -1;
  genfit::tools::invertMatrix(HCHt);
  chi2 = HCHt.Similarity(resNew);

  return chi2/ndf;
}

bool GFTrack::setInitialStateForExtrap(const int startPtID)
{
  genfit::AbsTrackRep* rep = _track->getCardinalRep();
  if(_track->getNumPoints() > 0)
  {
    genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(startPtID, rep);
    if(tp == nullptr) return false;

    genfit::KalmanFitterInfo* info = static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep));
    _propState.reset(new genfit::MeasuredStateOnPlane(info->getFittedState()));
  }
  else
  {
    _propState.reset(new genfit::MeasuredStateOnPlane(_fitstates[startPtID]));
  }

  return true;
}

void GFTrack::getExtrapPosMomCov(TVector3* pos, TVector3* mom, TMatrixDSym* cov)
{
  if(pos != nullptr && mom != nullptr)
  {
    if(cov != nullptr)
      _propState->getPosMomCov(*pos, *mom, *cov);
    else
      _propState->getPosMom(*pos, *mom);
  }
}

double GFTrack::swimToVertex(double z, TVector3* pos, TVector3* mom, TMatrixDSym* cov, bool biased)
{
  //Basic constants
  TVector3 pU(1., 0., 0.);
  TVector3 pV(0., 1., 0.);
  TVector3 pO(0., 0., z );
  
  TVectorD beamCenter(2);
  beamCenter[0] = X_BEAM; beamCenter[1] = Y_BEAM;
  TMatrixDSym beamCov(2);
  beamCov.Zero();
  beamCov(0, 0) = SIGX_BEAM*SIGX_BEAM; beamCov(1, 1) = SIGY_BEAM*SIGY_BEAM;

  double chi2 = -1.;
  try
  {
    double len = extrapolateToPlane(pO, pU, pV);
    if(fabs(len) > 6000.) throw len;
    if(!biased) getExtrapPosMomCov(pos, mom, cov);

    chi2 = updatePropState(beamCenter, beamCov);
    if(biased)  getExtrapPosMomCov(pos, mom, cov);
  }
  catch(genfit::Exception& e)
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": hypo test failed vertex @Z=" << z << ": " << e.what() << std::endl;
    return -1.;
  }
  catch(double len)
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": hypo test failed vertex @Z=" << z << ": " << len << std::endl;
    return -1.;
  }

  return chi2;
}

double GFTrack::getPOCA(SRecTrack* strack)
{
  std::vector<TVector3> pos, mom;
  try
  {
    pos = _trkrep->getPosSteps();
    mom = _trkrep->getMomSteps();
  }
  catch(...)
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": Extrapolation cache is not available." << std::endl;
    return 1.E6;
  }

  unsigned int nSteps = pos.size();
  if(nSteps != mom.size() || nSteps <= 1)
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": Extrapolation cache size illegle or mis-match " << nSteps << " vs. " << mom.size() << std::endl;
    return 1.E6;
  }

  double rx_min = fabs(pos[0].X() - X_BEAM);
  double ry_min = fabs(pos[0].Y() - Y_BEAM);
  double r2_min = rx_min*rx_min + ry_min*ry_min;
  unsigned int idx_r_min = 0;
  unsigned int idx_rx_min = 0;
  unsigned int idx_ry_min = 0;
  for(unsigned int i = 1; i < nSteps; ++i)
  {
    //TODO: add the accidental cross later, e.g. if(FMAGSTR*charge*mom[i].Px() < 0.) continue;
    double rx = fabs(pos[i].X() - X_BEAM);
    double ry = fabs(pos[i].Y() - Y_BEAM);
    double r2 = rx*rx + ry*ry;

    if(rx < rx_min)
    {
      rx_min = rx;
      idx_rx_min = i;
    }

    if(ry < ry_min)
    {
      ry_min = ry;
      idx_ry_min = i;
    }

    if(r2 < r2_min)
    {
      r2_min = r2;
      idx_r_min = i;
    }
  }

  if(strack != nullptr)
  {
    double dz_x = -(pos[idx_rx_min].X() - X_BEAM)/(mom[idx_rx_min].X()/mom[idx_rx_min].Z());
    strack->setXVertexPos(pos[idx_rx_min] + TVector3(mom[idx_rx_min].X()/mom[idx_rx_min].Z()*dz_x, mom[idx_rx_min].Y()/mom[idx_rx_min].Z()*dz_x, dz_x));
    strack->setXVertexMom(mom[idx_rx_min]);

    double dz_y = -(pos[idx_ry_min].Y() - Y_BEAM)/(mom[idx_ry_min].Y()/mom[idx_ry_min].Z());
    strack->setYVertexPos(pos[idx_ry_min] + TVector3(mom[idx_ry_min].X()/mom[idx_ry_min].Z()*dz_y, mom[idx_ry_min].Y()/mom[idx_ry_min].Z()*dz_y, dz_y));
    strack->setYVertexMom(mom[idx_ry_min]);
  }

  // Now we calculate the z position of POCA
  if(idx_r_min == nSteps - 1) return pos[idx_r_min].Z();

  double tx = mom[idx_r_min].X()/mom[idx_r_min].Z();
  double ty = mom[idx_r_min].Y()/mom[idx_r_min].Z();
  double x0 = pos[idx_r_min].X() - tx*pos[idx_r_min].Z();
  double y0 = pos[idx_r_min].Y() - ty*pos[idx_r_min].Z();
  double z_vtx = -(tx*(x0 - X_BEAM) + ty*(y0 - Y_BEAM))/(tx*tx + ty*ty);

  return z_vtx;
}

double GFTrack::getDumpPathLen()
{
  std::vector<genfit::MatStep> mat;
  try
  {
    mat = _trkrep->getSteps();
  }
  catch(...)
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": Extrapolation cache is not available." << std::endl;
    return -1.;
  }

  unsigned int nSteps = mat.size();
  if(nSteps <= 1)
  {
    std::cerr << __FILE__ << " " << __LINE__ << ": Extrapolation cache size illegle " << nSteps << std::endl;
    return -1.;
  }

  double ironLength = 0.;
  for(unsigned int i = 0; i < nSteps; ++i)
  {
    if(fabs(mat[i].material_.Z - 26.) < 1.E-3) ironLength += fabs(mat[i].stepSize_);
  }
  
  return ironLength;
}

void GFTrack::reset(bool updateSeed)
{
  if(updateSeed) _track->udpateSeed(0, _trkrep, true);
  _track->deleteFittedState(_track->getCardinalRep());
}

void GFTrack::checkConsistency()
{
  _track->checkConsistency();
  genfit::TrackPoint* tp = _track->getPointWithFitterInfo(0, _trkrep);
  if(tp != nullptr)
  {
    int fittedCharge = tp->getFitterInfo(_trkrep)->getFittedState().getCharge() > 0 ? 1 : -1;
    if(fittedCharge != _trkcand->getCharge())
    {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fitted sign is different " << fittedCharge << "  " << _trkcand->getCharge() << std::endl;
    }
  }
}

int GFTrack::getCharge() const
{
  genfit::TrackPoint* tp = _track->getPointWithFitterInfo(0, _trkrep);
  if(tp == nullptr) return _trkcand->getCharge();

  return tp->getFitterInfo(_trkrep)->getFittedState().getCharge() > 0 ? 1 : -1;
}

void GFTrack::setTracklet(const Tracklet& tracklet, double z_reference, bool wildseedcov)
{
  _trkcand = tracklet.Clone();
  _pdg = tracklet.getCharge() > 0 ? -13 : 13;
  _trkrep = new genfit::RKTrackRep(_pdg);

  TVectorD seed_state(6);
    
  TVector3 seed_mom = tracklet.getExpMomentum(z_reference);
  TVector3 seed_pos(tracklet.getExpPositionX(z_reference), tracklet.getExpPositionY(z_reference), z_reference);
  seed_state[0] = seed_pos.X();
  seed_state[1] = seed_pos.Y();
  seed_state[2] = seed_pos.Z();
  seed_state[3] = seed_mom.Px();
  seed_state[4] = seed_mom.Py();
  seed_state[5] = seed_mom.Pz();

  TMatrixDSym seed_cov(6);
  double uncertainty[6] = {10., 10., 10., 3., 3., 10.};
  for(int i = 0; i < 6; i++)
  {
    for(int j = 0; j < 6; j++)
    {
      seed_cov[i][j] = uncertainty[i]*uncertainty[j];
    }
  }


  if(!wildseedcov)
  {
    //TODO: Implement a smaller cov based on the chi^2 fit?
  }
  _track = new genfit::Track(_trkrep, seed_state, seed_cov);

  _measurements.clear();
  for(auto iter = tracklet.hits.begin(); iter != tracklet.hits.end(); ++iter)
  {
    if(iter->hit.index < 0) continue;

    GFMeasurement* meas = new GFMeasurement(*iter);
    _measurements.insert(meas);
  }

  addMeasurements(_measurements);
  //checkConsistency();
}

void GFTrack::postFitUpdate(bool updateMeasurements)
{
  if(!updateMeasurements) return;
  for(auto iter = _measurements.begin(); iter != _measurements.end(); ++iter)
  {
    (*iter)->postFitUpdate();
  }
}

SRecTrack GFTrack::getSRecTrack()
{
  //postFitUpdate();
  //The following steps are pretty hacky and should be considered as only a temporary solution
  SRecTrack strack;
  strack.setChisq(getChi2());
  for(auto iter = _measurements.begin(); iter != _measurements.end(); ++iter)
  {
    strack.insertHitIndex((*iter)->getBeforeFitHit().hit.index);

    const genfit::MeasuredStateOnPlane& fitstate = (*iter)->getTrackPoint()->getKalmanFitterInfo()->getFittedState(true);
    strack.insertGFState(fitstate);
    
    TVector3 pos, mom;
    fitstate.getPosMom(pos, mom);
    TMatrixD stateVec(5, 1);
    stateVec[0][0] = getCharge()/mom.Mag();
    stateVec[1][0] = mom.Px()/mom.Pz();
    stateVec[2][0] = mom.Py()/mom.Pz();
    stateVec[3][0] = pos.X();
    stateVec[4][0] = pos.Y();
    strack.insertStateVector(stateVec);
    strack.insertZ(pos.Z());

    TMatrixD cov(fitstate.getCov());
    strack.insertCovariance(cov);
    
    strack.insertChisq((*iter)->getTrackPoint()->getKalmanFitterInfo()->getSmoothedChi2());
  }

  //test Z_UPSTREAM and find POCA as rough vertex point
  strack.setChisqUpstream(swimToVertex(Z_UPSTREAM));
  double z_vtx = getPOCA(&strack);
  //strack.setDumpFacePos(TVector3(0., 0., getDumpPathLen()));

  //We then do the vertex hypothesis tests in order of z to maximize the use of cache
  TVector3 pos, mom;
  if(z_vtx >= Z_UPSTREAM && z_vtx < Z_TARGET)
  {
    //test Z_VERTEX
    strack.setChisqVertex(swimToVertex(z_vtx, &pos, &mom, nullptr, false));
    strack.setVertexPos(pos);
    strack.setVertexMom(mom);
  }

  //test Z_TARGET
  strack.setChisqTarget(swimToVertex(Z_TARGET, &pos, &mom, nullptr, false));
  strack.setTargetPos(pos);
  strack.setTargetMom(mom);

  if(z_vtx >= Z_TARGET && z_vtx < Z_DUMP)
  {
    //test Z_VERTEX
    strack.setChisqVertex(swimToVertex(z_vtx, &pos, &mom, nullptr, false));
    strack.setVertexPos(pos);
    strack.setVertexMom(mom);
  }

  //test Z_DUMP
  strack.setChisqDump(swimToVertex(Z_DUMP, &pos, &mom, nullptr, false));
  strack.setDumpPos(pos);
  strack.setDumpMom(mom);

  if(z_vtx >= Z_DUMP && z_vtx < 9999.)  // protection against crazy numbers 
  {
    //test Z_VERTEX
    strack.setChisqVertex(swimToVertex(z_vtx, &pos, &mom, nullptr, false));
    strack.setVertexPos(pos);
    strack.setVertexMom(mom);
  }

  /*
  //Find POCA to beamline -- it seems to be funky and mostly found some place way upstream or downstream
  // most likely because the cross product of the track direction and beam line direction is way too small 
  // on z axis to provide reasonable calculation of the POCA location. It's disabled for now.
  TVector3 ep1(0., 0., -499.);
  TVector3 ep2(0., 0., 0.);
  try
  {
    extrapolateToLine(ep1, ep2);
    TVectorD beamR(1); beamR(0) = 0.;
    TMatrixDSym beamC(1); beamC(0, 0) = 1000.;
    strack.setChisqVertex(updatePropState(beamR, beamC));
  }
  catch(genfit::Exception& e)
  {
    std::cerr << "Hypothesis test failed at beamline: " << e.what() << std::endl;
    print(0);
  }
  */

  strack.setKalmanStatus(1);
  return strack;
}

void GFTrack::print(unsigned int debugLvl)
{
  std::cout << "================================================ SGTrack ================================================" << std::endl;
  _track->getFitStatus(_trkrep)->Print();

  if(debugLvl > 0)
  {
    std::cout << "------------- Track candidate -----------" << std::endl;
    _trkcand->print();
  }
  else
  {
    std::cout << "Track has " << _track->getNumPoints() << " track points: " << std::endl;
    for(auto it = _measurements.begin(); it != _measurements.end(); ++it)
    {
      SignedHit& hit = (*it)->getBeforeFitHit();
      std::cout << hit.hit.index << "  " << p_geomSvc->getDetectorName(hit.hit.detectorID) << "(" << hit.hit.detectorID << ") " << hit.hit.elementID << " : ";
    }
    std::cout << std::endl;
  }

  try
  {
  if(debugLvl > 0)
  {
    std::cout << "------------- Fit Result ----------------" << std::endl;
    _track->getFittedState().Print();  //default to the first hit
  }
  else
  {
    TVector3 pos, mom;
    _track->getFittedState().getPosMom(pos, mom);
      std::cout << "Fitted pos (X,  Y, Z) = " << std::setprecision(6) << pos.X() << " " << pos.Y() << " " << pos.Z() << " cm\n";
      std::cout << "Fitted mom (Px,Py,Pz) = " << mom.X() << " " << mom.Y() << " " << mom.Z() << " GeV" << std::endl;
  }

  if(debugLvl < 2) return;
  for(auto iter = _measurements.begin(); iter != _measurements.end(); ++iter)
  {
    (*iter)->print(debugLvl-2);
  }
  }
  catch(...)
  {
    return;
  }

  if(debugLvl < 20) return;
  std::cout << "------------ GenFit Track ---------------" << std::endl;
  _track->Print();
}

}
