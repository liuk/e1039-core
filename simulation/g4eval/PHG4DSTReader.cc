// $Id: PHG4DSTReader.cc,v 1.11 2015/01/06 02:52:07 jinhuang Exp $

/*!
 * \file PHG4DSTReader.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.11 $
 * \date $Date: 2015/01/06 02:52:07 $
 */

#include "PHG4DSTReader.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <g4main/PHG4InEvent.h>
//#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitEval.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4VtxPointv1.h>

#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <TClonesArray.h>
#include <TTree.h>

#include <boost/foreach.hpp>
#include <map>
#include <set>
#include <cassert>

#include<sstream>

#include <memory>


using namespace std;

typedef PHG4Particlev2 part_type;
typedef PHG4HitEval hit_type;
typedef PHG4VtxPointv1 vertex_type;

PHG4DSTReader::PHG4DSTReader(const string &filename) :
    SubsysReco("PHG4DSTReader"), nblocks(0), _event(0), //
    _out_file_name(filename), /*_file(nullptr), */_T(nullptr), //
    _save_particle(true), _load_all_particle(false), _load_active_particle(
        true), _save_vertex(true)
{
  // TODO Auto-generated constructor stub

}

PHG4DSTReader::~PHG4DSTReader()
{
  cout << "PHG4DSTReader::destructor - Clean ups" << endl;

  if (_T)
    {
      _T->ResetBranchAddresses();
    }

  _records.clear();
}

int
PHG4DSTReader::Init(PHCompositeNode*)
{

  const int arr_size = 100;

//  BOOST_FOREACH(string nodenam, _node_postfix)
  for (vector<string>::const_iterator it = _node_postfix.begin();
      it != _node_postfix.end(); ++it)
    {
      const char * class_name = hit_type::Class()->GetName();

      const string & nodenam = *it;

      string hname = "G4HIT_" + nodenam;
//      _node_name.push_back(hname);
      cout << "PHG4DSTReader::Init - saving hits from node: " << hname << " - "
          << class_name << endl;

      record rec;
      rec._cnt = 0;
      rec._name = hname;
      rec._arr = std::make_shared<TClonesArray>(class_name, arr_size);
      rec._arr_ptr = rec._arr.get();
      rec._type = record::typ_hit;

      _records.push_back(rec);

      nblocks++;
    }

  if (_save_particle)
    {
      //save particles

      const char * class_name = part_type::Class()->GetName();

      cout << "PHG4DSTReader::Init - saving Particles node:" << class_name
          << endl;

      record rec;
      rec._cnt = 0;
      rec._name = "PHG4Particle";
      rec._arr = std::make_shared<TClonesArray>(class_name, arr_size);
      rec._arr_ptr = rec._arr.get();
      rec._type = record::typ_part;

      _records.push_back(rec);

      nblocks++;
    }

  if (_save_vertex)
    {
      //save particles

      const char * class_name = vertex_type::Class()->GetName();

      cout << "PHG4DSTReader::Init - saving vertex for particles under node:"
          << class_name << endl;

      record rec;
      rec._cnt = 0;
      rec._name = "PHG4VtxPoint";
      rec._arr = std::make_shared<TClonesArray>(class_name, arr_size);
      rec._arr_ptr = rec._arr.get();
      rec._type = record::typ_vertex;

      _records.push_back(rec);

      nblocks++;
    }

  cout << "PHG4DSTReader::Init - requested " << nblocks << " nodes" << endl;

  build_tree();

  return 0;
}

void
PHG4DSTReader::build_tree()
{
  cout << "PHG4DSTReader::build_tree - output to " << _out_file_name << endl;

  static const int BUFFER_SIZE = 32000;

  // open TFile
  PHTFileServer::get().open(_out_file_name, "RECREATE");

  _T = new TTree("T", "PHG4DSTReader");

  nblocks = 0;
  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)
    {
      record & rec = *it;

      cout << "PHG4DSTReader::build_tree - Add " << rec._name << endl;

      const string name_cnt = "n_" + rec._name;
      const string name_cnt_desc = name_cnt + "/I";
      _T->Branch(name_cnt.c_str(), &(rec._cnt), name_cnt_desc.c_str(),
          BUFFER_SIZE);
      _T->Branch(rec._name.c_str(), &(rec._arr_ptr), BUFFER_SIZE, 99);

      nblocks++;
    }

  cout << "PHG4DSTReader::build_tree - added " << nblocks << " nodes" << endl;

  _T->SetAutoSave(16000);
}

int
PHG4DSTReader::process_event(PHCompositeNode* topNode)
{

//  const double significand = _event / TMath::Power(10, (int) (log10(_event)));
//
//  if (fmod(significand, 1.0) == 0 && significand <= 10)
//    cout << "PHG4DSTReader::process_event - " << _event << endl;
  _event++;

  //clean ups
  _particle_set.clear();
  _vertex_set.clear();

  PHG4TruthInfoContainer* truthInfoList = findNode::getClass<
      PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthInfoList)
    {
      if (_event < 2)
        cout
            << "PHG4DSTReader::process_event - Error - can not find node G4TruthInfo. Quit processing!"
            << endl;
      return Fun4AllReturnCodes::DISCARDEVENT;
    }

  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)
    {
      record & rec = *it;

      rec._cnt = 0;
      assert(rec._arr.get() == rec._arr_ptr);
      assert(rec._arr.get());
      rec._arr->Clear();

      if (rec._type == record::typ_hit)
        {
          if (Verbosity() >= 2)
            cout << "PHG4DSTReader::process_event - processing " << rec._name
                << endl;

          PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode,
              rec._name);
          if (!hits)
            {
              if (_event < 2)
                cout
                    << "PHG4DSTReader::process_event - Error - can not find node "
                    << rec._name << endl;

            }
          else
            {
              PHG4HitContainer::ConstRange hit_range = hits->getHits();

              if (Verbosity() >= 2)
                cout << "PHG4DSTReader::process_event - processing "
                    << rec._name << " and received " << hits->size() << " hits"
                    << endl;

              for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first;
                  hit_iter != hit_range.second; hit_iter++)
                {
                  PHG4Hit * hit = hit_iter->second;

//                  hit_type * hit = dynamic_cast<hit_type *>(hit_raw);
//
                  assert(hit);

                  new ((*(rec._arr.get()))[rec._cnt]) hit_type(*hit);

                  hit_type * new_hit =
                      dynamic_cast<hit_type *>(rec._arr.get()->At(rec._cnt));
                  assert(new_hit);

//                  for (int i = 0; i < 2; i++)
//                    {
//                      new_hit->set_x(i, hit->get_x(i));
//                      new_hit->set_y(i, hit->get_y(i));
//                      new_hit->set_z(i, hit->get_z(i));
//                      new_hit->set_t(i, hit->get_t(i));
//                    }
//                  new_hit->set_edep(hit->get_edep());
//                  new_hit->set_layer(hit->get_layer());
//                  new_hit->set_hit_id(hit->get_hit_id());
//                  new_hit->set_trkid(hit->get_trkid());

//                  *new_hit = (*hit);

                  if (_load_active_particle)
                    _particle_set.insert(hit->get_trkid());

                  if (Verbosity() >= 2)
                    cout << "PHG4DSTReader::process_event - processing "
                        << rec._name << " and hit " << hit->get_hit_id()
                        << " with track id " << hit->get_trkid() << endl;

                  rec._cnt++;
                }
            } // if (!hits)
        }
      else if (rec._type == record::typ_part)
        {

          map<int, PHG4Particle *>::const_iterator particle_iter;

          if (_load_all_particle)
            {
              static bool once = true;
              if (once)
                {
                  cout
                      << "PHG4DSTReader::process_event - will load all particle from G4TruthInfo"
                      << endl;

                  once = false;
                }

              for (particle_iter = truthInfoList->GetMap().begin();
                  particle_iter != truthInfoList->GetMap().end();
                  particle_iter++)
                {

                  _particle_set.insert(particle_iter->first);

                }

            } //          if (_load_all_particle)
          else
            {
              static bool once = true;
              if (once)
                {
                  cout
                      << "PHG4DSTReader::process_event - will load primary particle from G4TruthInfo"
                      << endl;

                  once = false;
                }

              PHG4TruthInfoContainer::ConstRange primary_range =
                  truthInfoList->GetPrimaryParticleRange();

              for (PHG4TruthInfoContainer::ConstIterator particle_iter =
                  primary_range.first; particle_iter != primary_range.second;
                  ++particle_iter)
                {

                  _particle_set.insert(particle_iter->first);

                } //          if (_load_all_particle) else
            }
          for (PartSet_t::const_iterator i = _particle_set.begin();
              i != _particle_set.end(); i++)
            {
              particle_iter = truthInfoList->GetMap().find(*i);
              if (particle_iter == truthInfoList->GetMap().end())
                {

                  cout
                      << "PHG4DSTReader::process_event - ERROR - can not find particle/track ID "
                      << *i << " in G4TruthInfo" << endl;

                  continue;
                }

              PHG4Particle * part = particle_iter->second;
              add_particle(rec, part);
            } // for(PartSet_t::const_iterator i = _particle_set.begin();i!=_particle_set.end();i++)

        } //      if (rec._type == record::typ_part)
      else if (rec._type == record::typ_vertex)
        {
          static bool once = true;
          if (once)
            {
              cout
                  << "PHG4DSTReader::process_event - will load vertex from G4TruthInfo"
                  << endl;

              once = false;
            }

          for (PartSet_t::const_iterator i = _vertex_set.begin();
              i != _vertex_set.end(); ++i)
            {

              PHG4VtxPoint * v = truthInfoList->GetVtx(*i);
              if (!v)
                {
                  cout
                      << "PHG4DSTReader::process_event - ERROR - can not find vertex ID "
                      << *i << " in G4TruthInfo" << endl;

                  continue;
                }

              new ((*(rec._arr.get()))[rec._cnt]) vertex_type();

              if (Verbosity() >= 2)
                cout << "PHG4DSTReader::process_event - saving vertex id "
                    << *i  << endl;

              vertex_type * new_v =
                  static_cast<vertex_type *>(rec._arr.get()->At(rec._cnt));
              assert(new_v);

              new_v->set_x(v->get_x());
              new_v->set_y(v->get_y());
              new_v->set_z(v->get_z());
              new_v->set_t(v->get_t());
              new_v->set_id(v->get_id());

              rec._cnt++;
            } // else if (rec._type == record::typ_vertex)

        } //          if (_load_all_particle)

    } //  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)

  if (_T)
    _T->Fill();

  return 0;
} //  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)

void
PHG4DSTReader::add_particle(PHG4DSTReader::record & rec, PHG4Particle * part)
{

  assert(part);

//              if (Verbosity() >= 2)
//                cout << "PHG4DSTReader::process_event - Particle type is "
//                    << p_raw->GetName() << " - " << p_raw->ClassName()
//                    << " get_track_id = " << p_raw->get_track_id() << endl;

//              part_type * part = dynamic_cast<part_type *>(p_raw);

//              assert(part);

  new ((*(rec._arr.get()))[rec._cnt]) part_type();

  part_type * new_part = static_cast<part_type *>(rec._arr.get()->At(rec._cnt));
  assert(new_part);

  new_part->set_track_id(part->get_track_id());
  new_part->set_vtx_id(part->get_vtx_id());
  new_part->set_parent_id(part->get_parent_id());
  new_part->set_primary_id(part->get_primary_id());
  new_part->set_name(part->get_name());
  new_part->set_pid(part->get_pid());
  new_part->set_px(part->get_px());
  new_part->set_py(part->get_py());
  new_part->set_pz(part->get_pz());
  new_part->set_e(part->get_e());

  _vertex_set.insert(part->get_vtx_id());

  rec._cnt++;
}

int
PHG4DSTReader::End(PHCompositeNode * /*topNode*/)
{
  cout << "PHG4DSTReader::End - Clean ups" << endl;

  if (_T)
    {
      PHTFileServer::get().cd(_out_file_name);
      _T->Write();
      _T->ResetBranchAddresses();
    }

  _records.clear();

  return 0;
}
