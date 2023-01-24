/// OnlMonTrigV1495.C
#include <sstream>
#include <iomanip>
#include <cstring>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include <UtilAna/UtilSQHit.h>
#include <UtilAna/UtilHist.h>
#include "OnlMonTrigV1495.h"
using namespace std;

OnlMonTrigV1495::OnlMonTrigV1495()
{
  NumCanvases(7);
  Name("OnlMonTrigV1495" ); 
  Title("Trigger Analysis" );
  
  event_cnt = 1;
}

int OnlMonTrigV1495::InitOnlMon(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int OnlMonTrigV1495::InitRunOnlMon(PHCompositeNode* topNode)
{
  SetDet();

  GeomSvc* geom = GeomSvc::instance();
  ostringstream oss;
 
  const double DT = 5.0; // 1 ns per single count of v1495 TDC
  int NT    = 100;
  double T0 = 100.5*DT;
  double T1 = 200.5*DT;

  int num_tot_ele = 0;
  //Loop through hodoscopes 
  for (int i_det = 0; i_det < N_DET; i_det++) {
    string name = list_det_name[i_det];
    int  det_id = list_det_id  [i_det];
    int n_ele  = geom->getPlaneNElements(det_id);
    num_tot_ele += n_ele;
    if (det_id <= 0 || n_ele <= 0) {
      cout << "OnlMonTrigV1495::InitRunOnlMon():  Invalid det_id or n_ele: " 
           << det_id << " " << n_ele << " at name = " << name << "." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

   /* oss.str("");
    oss << "h2_time_ele_raw_" << i_det;
    h2_time_ele_raw[i_det] = new TH2D(oss.str().c_str(), "", n_ele, 0.5, n_ele+0.5, 80, 720.5, 880.5);
    oss.str("");
    oss << name << ";Element ID;tdcTime;Hit count";
    h2_time_ele_raw[i_det]->SetTitle(oss.str().c_str());
    
    RegisterHist(h2_time_ele_raw[i_det]); */

    oss.str("");
    oss << "h1_ele_" << i_det;
    h1_ele[i_det] = new TH1D(oss.str().c_str(), "", n_ele, 0.5, n_ele+0.5);
    oss.str("");
    oss << name << ";Element ID;Hit count";
    h1_ele[i_det]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_ele_in_" << i_det;
    h1_ele_in[i_det] = new TH1D(oss.str().c_str(), "", n_ele, 0.5, n_ele+0.5);
    oss.str("");
    oss << name << ";Element ID;In-time hit count";
    h1_ele_in[i_det]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_time_" << i_det;
    h1_time[i_det] = new TH1D(oss.str().c_str(), "", NT, T0, T1);

    oss.str("");
    oss << name << ";tdcTime;Hit count";
    h1_time[i_det]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_time_in_" << i_det;
    h1_time_in[i_det] = new TH1D(oss.str().c_str(), "", NT, T0, T1);
    oss.str("");
    oss << name << ";tdcTime;In-time hit count";
    h1_time_in[i_det]->SetTitle(oss.str().c_str());


    oss.str("");
    oss << "h2_time_ele_fpga_" << i_det;
    h2_time_ele_fpga[i_det] = new TH2D(oss.str().c_str(), "", n_ele, 0.5, n_ele+0.5, 80, 720.5, 880.5);
    oss.str("");
    oss << name << ";Element ID;tdcTime;Hit count";
    h2_time_ele_fpga[i_det]->SetTitle(oss.str().c_str());

    RegisterHist(h1_ele    [i_det]);
    RegisterHist(h1_ele_in [i_det]);
    RegisterHist(h1_time   [i_det]);
    RegisterHist(h1_time_in[i_det]);
    RegisterHist(h2_time_ele_fpga[i_det]);
  }
    
    oss.str("");
    oss << "h2_trig_time_" << 1;
    h2_trig_time = new TH2D(oss.str().c_str(), "", 10, 0.5, 10.5, 150, 800.5, 1100.5);
    oss.str("");
    oss << "Trigger Timing" << ";Trigger;tdcTime;Hit count";
    h2_trig_time->SetTitle(oss.str().c_str());

    h2_trig_time->GetXaxis()->SetBinLabel( 1, "FPGA1");
    h2_trig_time->GetXaxis()->SetBinLabel( 2, "FPGA2");
    h2_trig_time->GetXaxis()->SetBinLabel( 3, "FPGA3");
    h2_trig_time->GetXaxis()->SetBinLabel( 4, "FPGA4");
    h2_trig_time->GetXaxis()->SetBinLabel( 5, "FPGA5");
    h2_trig_time->GetXaxis()->SetBinLabel( 6, "NIM1");
    h2_trig_time->GetXaxis()->SetBinLabel( 7, "NIM2");
    h2_trig_time->GetXaxis()->SetBinLabel( 8, "NIM3");
    h2_trig_time->GetXaxis()->SetBinLabel( 9, "NIM4");
    h2_trig_time->GetXaxis()->SetBinLabel(10, "NIM5");

    oss.str("");
    oss << "h1_trig_" << 1;
    h1_trig = new TH1D(oss.str().c_str(), "", 10, 0.5, 10.5);
    oss.str("");
    oss << "Trigger Status" << ";Trigger;N of events";
    h1_trig->SetTitle(oss.str().c_str());

    h1_trig->GetXaxis()->SetBinLabel( 1, "FPGA1");
    h1_trig->GetXaxis()->SetBinLabel( 2, "FPGA2");
    h1_trig->GetXaxis()->SetBinLabel( 3, "FPGA3");
    h1_trig->GetXaxis()->SetBinLabel( 4, "FPGA4");
    h1_trig->GetXaxis()->SetBinLabel( 5, "FPGA5");
    h1_trig->GetXaxis()->SetBinLabel( 6, "NIM1");
    h1_trig->GetXaxis()->SetBinLabel( 7, "NIM2");
    h1_trig->GetXaxis()->SetBinLabel( 8, "NIM3");
    h1_trig->GetXaxis()->SetBinLabel( 9, "NIM4");
    h1_trig->GetXaxis()->SetBinLabel(10, "NIM5");

    oss.str("");
    oss << "h2_RF_" << 1;
    h2_RF = new TH2D(oss.str().c_str(), "",  300, 300.5, 900.5,  9, 0.5, 9.5);
    oss.str("");
    oss << "RF TDC" << ";tdcTime;RF Board;Hit count";
    h2_RF->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h2_RF_raw_" << 1;
    h2_RF_raw = new TH2D(oss.str().c_str(), "",  300, 300.5, 900.5,  9, 0.5, 9.5);
    oss.str("");
    oss << "RF TDC raw" << ";tdcTime;RF Board;Hit count";
    h2_RF_raw->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h2_fpga_nim_time_" << 1;
    h2_fpga_nim_time = new TH2D(oss.str().c_str(), "", 100, 1000.5, 1100.5, 100, 1000.5, 1100.5);
    oss.str("");
    oss << "FPGA & NIM Timing" << ";NIM tdcTime;FPGA tdcTime;Hit count";
    h2_fpga_nim_time->SetTitle(oss.str().c_str());
 
    for(int i = 0; i < 8; i++){
      oss.str("");
      oss << "RF_proj_" << i;
      RF_proj[i] = new TH1D(oss.str().c_str(), "",160 , 720.5, 880.5);
      oss.str("");
      oss << "RF TDC Projection" << ";TDC RF;Hit count";
      RF_proj[i]->SetTitle(oss.str().c_str());
    
      RegisterHist(RF_proj[i]);
    } 
    
    RegisterHist(h2_RF);
    RegisterHist(h2_RF_raw);
    RegisterHist(h1_trig );
    RegisterHist(h2_trig_time);
    RegisterHist(h2_fpga_nim_time);
 // cout << "REGISTERING HISTOGRAMS" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int OnlMonTrigV1495::ProcessEventOnlMon(PHCompositeNode* topNode)
{ 

  SQEvent*      evt     = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  SQHitVector*  hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  SQHitVector*  trig_hit_vec = findNode::getClass<SQHitVector>(topNode, "SQTriggerHitVector");
  if (!evt || !hit_vec  || !trig_hit_vec) return Fun4AllReturnCodes::ABORTEVENT;

  //Determine whether event is FPGA1-4 
  int is_FPGA_event = (evt->get_trigger(SQEvent::MATRIX1) || evt->get_trigger(SQEvent::MATRIX2) ||evt->get_trigger(SQEvent::MATRIX3)||evt->get_trigger(SQEvent::MATRIX4) ) ? 1 : 0; 

//RF *************************************************************************************** 
  //Looping through RF data to determine RF buckets 
  auto vec1 = UtilSQHit::FindTriggerHitsFast(evt, trig_hit_vec, "RF");
  int count = 0;
  for(auto it = vec1->begin(); it != vec1->end(); it++){
    double tdc_time = (*it)->get_tdc_time();
    float element = (*it)->get_element_id();
    
        
    h2_RF_raw->Fill(tdc_time,element);  
    if(is_FPGA_event){
	h2_RF->Fill(tdc_time,element); 
        
        if(count < 8){
           //Saving RF in 1d histograms to find max 
           RF_proj[count]->Fill(tdc_time);
        }
    }else{
    
    }
    
    count ++;
  }

//TDC 4 ************************************************************************************

  auto vec0 = UtilSQHit::FindHitsFast(evt, hit_vec, "AfterInhMatrix");
  for (auto it = vec0->begin(); it != vec0->end(); it++) {
    h2_trig_time->Fill((*it)->get_element_id(),(*it)->get_tdc_time());
  }

  auto vec2 = UtilSQHit::FindHitsFast(evt, hit_vec, "BeforeInhNIM");
  for (auto it = vec2->begin(); it != vec2->end(); it++) {
    h2_trig_time->Fill((*it)->get_element_id()+5,(*it)->get_tdc_time()); // element +5 so nim index start at 6 in histo
  }
 
  if(evt->get_trigger(SQEvent::MATRIX1) && evt->get_trigger(SQEvent::NIM4)){
   // cout << "FPGA1 && NIM4 ***" << endl;
    for(auto it0 = vec0->begin(); it0 != vec0->end(); it0++){ //FPGA
      double ele_FPGA1 = (*it0)->get_element_id();
      double time_FPGA1 = (*it0)->get_tdc_time();
      for(auto it1 = vec2->begin(); it1 != vec2->end(); it1++){//NIM
        double ele_NIM4 = (*it1)->get_element_id();
        double time_NIM4 = (*it1)->get_tdc_time();
        if(ele_FPGA1 == 1 && ele_NIM4 == 4){
          h2_fpga_nim_time->Fill(time_NIM4,time_FPGA1);
        }
      }
    }
    
  }
//Hodoscope Hits *************************************************************************** 
  
  //Loop through all hodoscopes
  for (int i_det = 0; i_det < N_DET; i_det++) {
    int det_id = list_det_id[i_det];

    auto vec = UtilSQHit::FindTriggerHitsFast(evt, trig_hit_vec, det_id);
    //Loop through hodoscope hits from each hodoscope
    for (auto it = vec->begin(); it != vec->end(); it++) {
      if ((*it)->get_level() != 1) continue; //switched m_lvl for 1
      int    ele  = (*it)->get_element_id();
      double time = (*it)->get_tdc_time();
     // h2_time_ele_raw[i_det]->Fill(ele, time); 
       
      //ONLY FPGA 1 THROUGH 4
      if(is_FPGA_event){  

        h1_ele     [i_det]->Fill(ele );
        h1_time    [i_det]->Fill(time);
        if ((*it)->is_in_time()) {
          h1_ele_in [i_det]->Fill(ele );
          h1_time_in[i_det]->Fill(time);
        } 
        
        h2_time_ele_fpga[i_det]->Fill(ele,time);
      }
       

      if (evt->get_trigger(SQEvent::MATRIX1)) h1_trig->Fill( 1);
      if (evt->get_trigger(SQEvent::MATRIX2)) h1_trig->Fill( 2);
      if (evt->get_trigger(SQEvent::MATRIX3)) h1_trig->Fill( 3);
      if (evt->get_trigger(SQEvent::MATRIX4)) h1_trig->Fill( 4);
      if (evt->get_trigger(SQEvent::MATRIX5)) h1_trig->Fill( 5);
      if (evt->get_trigger(SQEvent::NIM1   )) h1_trig->Fill( 6);
      if (evt->get_trigger(SQEvent::NIM2   )) h1_trig->Fill( 7);
      if (evt->get_trigger(SQEvent::NIM3   )) h1_trig->Fill( 8);
      if (evt->get_trigger(SQEvent::NIM4   )) h1_trig->Fill( 9);
      if (evt->get_trigger(SQEvent::NIM5   )) h1_trig->Fill(10); 
    }
  }
  
  if(event_cnt == 10000 || event_cnt == 1){
   // cout << "NEW EVENT " << event_cnt << endl;
  }
  event_cnt++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int OnlMonTrigV1495::EndOnlMon(PHCompositeNode* topNode)
{ 
  return Fun4AllReturnCodes::EVENT_OK;
}

int OnlMonTrigV1495::FindAllMonHist()
{

 // cout << "FIND ALL MON HIST PART" << endl;
  ostringstream oss;
  oss.str("");
  oss << "h1_trig_" << 1;
  h1_trig = FindMonHist(oss.str().c_str());
  if (! h1_trig) return 1;
  
  //Loop through hodoscopes
  for (int i_det = 0; i_det < N_DET; i_det++) {
    /*oss.str("");
    oss << "h2_time_ele_raw_" << i_det;
    h2_time_ele_raw[i_det] = (TH2*)FindMonHist(oss.str().c_str());
    if (! h2_time_ele_raw[i_det]) return 1; 
 */  

    oss.str("");
    oss << "h1_ele_" << i_det;
    h1_ele[i_det] = FindMonHist(oss.str().c_str());
    if (! h1_ele[i_det]) return 1;
    
    oss.str("");
    oss << "h1_ele_in_" << i_det;
    h1_ele_in[i_det] = FindMonHist(oss.str().c_str());
    if (! h1_ele_in[i_det]) return 1;
    
    oss.str("");
    oss << "h1_time_" << i_det;
    h1_time[i_det] = FindMonHist(oss.str().c_str());
    if (! h1_time[i_det]) return 1;
    
    oss.str("");
    oss << "h1_time_in_" << i_det;
    h1_time_in[i_det] = FindMonHist(oss.str().c_str());
    if (! h1_time_in[i_det]) return 1;
    
    oss.str("");
    oss << "h2_time_ele_fpga_" << i_det;
    h2_time_ele_fpga[i_det] = (TH2*)FindMonHist(oss.str().c_str());
    if (! h2_time_ele_fpga[i_det]) return 1;
  }

  
  oss.str("");
  oss << "h2_trig_time_" << 1;
  h2_trig_time = (TH2*)FindMonHist(oss.str().c_str());
  if (! h2_trig_time) return 1;

  oss.str("");
  oss << "h2_RF_" << 1;
  h2_RF = (TH2*)FindMonHist(oss.str().c_str());
  if (! h2_RF) return 1;

  oss.str("");
  oss << "h2_RF_raw_" << 1;
  h2_RF_raw = (TH2*)FindMonHist(oss.str().c_str());
  if (! h2_RF_raw) return 1;
 
  oss.str("");
  oss << "h2_fpga_nim_time_" << 1;
  h2_fpga_nim_time = (TH2*)FindMonHist(oss.str().c_str());
  if (! h2_fpga_nim_time) return 1;
  return 0;
}

int OnlMonTrigV1495::DrawMonitor()
{
  //DRAWING HISTOGRAMS ON .PNG FILES ******************************************
//  cout << "DRAWING HISTOGRAMS ON PNG" << endl;
  int binmax[8];
  double x[8];
  //Determine maximum value from projection histos to determine RF rising edges
  for(int i = 0; i < 8; i++){
   binmax[i] = RF_proj[i]->GetMaximumBin(); 
   x[i] = RF_proj[i]->GetXaxis()->GetBinCenter(binmax[i]); 
  // cout << x[i] << endl;
   //Lines to show RF buckets on 2d histo
   proj_line[i] = new TLine(0.5,x[i],16.5,x[i]);
   proj_line_H1[i]= new TLine(0.5,x[i],23.5,x[i]);
  }

  //can0->SetStatus(OnlMonCanvas::OK);
  OnlMonCanvas* can0 = GetCanvas(0);
  TPad* pad0 = can0->GetMainPad();
  pad0->Divide(1, 2);
  TVirtualPad* pad00 = pad0->cd(1);
  pad00->SetGrid();
 // for (int i_det = 0; i_det < N_DET; i_det++) {
  h2_RF->Draw("colz");
  //ostringstream oss2;
  //oss2 << "pr_" << h2_RF->GetName();
  //TProfile* pr2 = h2_RF->ProfileX(oss2.str().c_str());
  //pr2->SetLineColor(kBlack);
  //pr2->Draw("E1same");
 // }
  //can1->SetStatus(OnlMonCanvas::OK);
 
  TVirtualPad* pad01 = pad0->cd(2);
  pad01->SetGrid();
  h2_RF_raw->Draw("colz");
/*  h2_phase->Draw("colz");
  ostringstream oss;
  oss << "pr_" << h2_phase->GetName();
  TProfile* pr = h2_phase->ProfileX(oss.str().c_str());
  pr->SetLineColor(kBlack);
  pr->Draw("E1same");
*/

  OnlMonCanvas* can1 = GetCanvas(1);
  TPad* pad1 = can1->GetMainPad();
  pad1->Divide(2,2);
  for(int i_det = 0; i_det < 2; i_det++){
    TVirtualPad* pad11 = pad1->cd(2*i_det+1);
    pad11->SetGrid();
    h1_ele[i_det]->SetLineColor(kBlack);
    h1_ele[i_det]->Draw();
    h1_ele_in[i_det]->SetLineColor(kBlue);
    h1_ele_in[i_det]->SetFillColor(kBlue-7);
    h1_ele_in[i_det]->Draw("same");
   
    TVirtualPad* pad10 = pad1->cd(2*i_det+2);
    pad10->SetGrid();
    h2_time_ele_fpga[i_det]->Draw("colz");
    ostringstream oss;
    oss << "pr_" << h2_time_ele_fpga[i_det]->GetName();
    TProfile* pr = h2_time_ele_fpga[i_det]->ProfileX(oss.str().c_str());
    pr->SetLineColor(kBlack);
    pr->Draw("E1same");
    for(int i = 0; i < 8; i++){
      proj_line_H1[i]->SetLineStyle(2);
      proj_line_H1[i]->SetLineWidth(3);
      proj_line_H1[i]->SetLineColor(kRed);
      proj_line_H1[i]->Draw();
    }
  }

  OnlMonCanvas* can2 = GetCanvas(2);
  TPad* pad2 = can2->GetMainPad();
  pad2->Divide(2,2); 
  for(int i_det = 2; i_det < 4; i_det++){
    TVirtualPad* pad22 = pad2->cd(2*(i_det-2)+1);
    pad22->SetGrid();
    h1_ele[i_det]->SetLineColor(kBlack);
    h1_ele[i_det]->Draw();
    h1_ele_in[i_det]->SetLineColor(kBlue);
    h1_ele_in[i_det]->SetFillColor(kBlue-7);
    h1_ele_in[i_det]->Draw("same");

    TVirtualPad* pad20 = pad2->cd(2*(i_det-2)+2);
    pad20->SetGrid();
    h2_time_ele_fpga[i_det]->Draw("colz");
    ostringstream oss;
    oss << "pr_" << h2_time_ele_fpga[i_det]->GetName();
    TProfile* pr = h2_time_ele_fpga[i_det]->ProfileX(oss.str().c_str());
    pr->SetLineColor(kBlack);
    pr->Draw("E1same");
    for(int i = 0; i < 8; i++){
      proj_line[i]->SetLineStyle(2);
      proj_line[i]->SetLineWidth(3);
      proj_line[i]->SetLineColor(kRed);
      proj_line[i]->Draw();
    }
  }
  OnlMonCanvas* can3 = GetCanvas(3);
  TPad* pad3 = can3->GetMainPad();
  pad3->Divide(2,2);
  for(int i_det = 4; i_det < 6; i_det++){
    TVirtualPad* pad33 = pad3->cd(2*(i_det-4)+1);
    pad33->SetGrid();
    h1_ele[i_det]->SetLineColor(kBlack);
    h1_ele[i_det]->Draw();
    h1_ele_in[i_det]->SetLineColor(kBlue);
    h1_ele_in[i_det]->SetFillColor(kBlue-7);
    h1_ele_in[i_det]->Draw("same");


    TVirtualPad* pad30 = pad3->cd(2*(i_det-4)+2);
    pad30->SetGrid();
    h2_time_ele_fpga[i_det]->Draw("colz");
    ostringstream oss;
    oss << "pr_" << h2_time_ele_fpga[i_det]->GetName();
    TProfile* pr = h2_time_ele_fpga[i_det]->ProfileX(oss.str().c_str());
    pr->SetLineColor(kBlack);
    pr->Draw("E1same");
    for(int i = 0; i < 8; i++){
      proj_line[i]->SetLineStyle(2);
      proj_line[i]->SetLineWidth(3);
      proj_line[i]->SetLineColor(kRed);
      proj_line[i]->Draw();
    }
  }

  OnlMonCanvas* can4 = GetCanvas(4);
  TPad* pad4 = can4->GetMainPad();
  pad4->Divide(2,2);
  for(int i_det = 6; i_det < 8; i_det++){
    TVirtualPad* pad44 = pad4->cd(2*(i_det-6)+1);
    pad44->SetGrid();
    h1_ele[i_det]->SetLineColor(kBlack);
    h1_ele[i_det]->Draw();
    h1_ele_in[i_det]->SetLineColor(kBlue);
    h1_ele_in[i_det]->SetFillColor(kBlue-7);
    h1_ele_in[i_det]->Draw("same");

    TVirtualPad* pad40 = pad4->cd(2*(i_det-6)+2);
    pad40->SetGrid();
    h2_time_ele_fpga[i_det]->Draw("colz");
    ostringstream oss;
    oss << "pr_" << h2_time_ele_fpga[i_det]->GetName();
    TProfile* pr = h2_time_ele_fpga[i_det]->ProfileX(oss.str().c_str());
    pr->SetLineColor(kBlack);
    pr->Draw("E1same");
    for(int i = 0; i < 8; i++){
      proj_line[i]->SetLineStyle(2);
      proj_line[i]->SetLineWidth(3);
      proj_line[i]->SetLineColor(kRed);
      proj_line[i]->Draw();
    }
  }

  OnlMonCanvas* can5 = GetCanvas(5);
  TPad* pad5 = can5->GetMainPad();
  pad5->Divide(1,2);
  TVirtualPad* pad50 = pad5->cd(1);
  pad50->SetGrid();
  h1_trig->SetLineColor(kRed);
  h1_trig->Draw();
  
  TVirtualPad* pad51 = pad5->cd(2);
  pad51->SetGrid();
  h2_trig_time->Draw("colz");
  ostringstream oss; 
  oss << "pr_" << h2_trig_time->GetName();
  TProfile* pr = h2_trig_time->ProfileX(oss.str().c_str());
  pr->SetLineColor(kBlack);
  pr->Draw("E1same"); 

  OnlMonCanvas* can6 = GetCanvas(6);
  TPad* pad6 = can6->GetMainPad();
  pad6->Divide(1,2);
  TVirtualPad* pad60 = pad6->cd(1);
  pad60->SetGrid();
  h2_fpga_nim_time->Draw("colz");
 /* ostringstream oss1;
  oss1 << "pr_" << h2_fpga_nim_time->GetName();
  TProfile* pr1 = h2_fpga_nim_time->ProfileX(oss1.str().c_str());
  pr1->SetLineColor(kBlack);
  pr1->Draw("E1same");*/


  return 0;
}

void OnlMonTrigV1495::SetDet()
{
  list_det_name[0] = "H1T";
  list_det_name[1] = "H1B";
  list_det_name[2] = "H2T";
  list_det_name[3] = "H2B";
  list_det_name[4] = "H3T";
  list_det_name[5] = "H3B";
  list_det_name[6] = "H4T";
  list_det_name[7] = "H4B";
   
  GeomSvc* geom = GeomSvc::instance();
  for (int ii = 0; ii < N_DET; ii++) {
    list_det_id[ii] = geom->getDetectorID(list_det_name[ii]);
  }
}

