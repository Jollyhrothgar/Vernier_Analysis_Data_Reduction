#include "VernierDSTReduction.h"

//General PHENIX tools
#include <phool.h>
#include <getClass.h>
#include <PHCompositeNode.h>
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <TriggerHelper.h>

//Data classes I am using in analysis
#include <PHGlobal.h>
#include <PreviousEvent.h>
#include <SpinDataEventOutv2.h>
#include <SpinEvent.h>
#include <BbcOutv1.h>
#include <BbcRaw.h>
#include <ZdcOutv2.h>
#include <ZdcRaw.h>
#include <TrigLvl1.h>
#include <Bbc.hh>

//Root header files
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

using namespace std;
using namespace findNode;

const int VernierDSTReduction::histosize = VernierDSTReduction::asize;
const double VernierDSTReduction::historange = (double)histosize;

// ================================================================
VernierDSTReduction::VernierDSTReduction(const string &outfile) : 
    SubsysReco ("VernierDSTReduction"),
    OutFileName(outfile),
    fout(NULL),
    finpol(NULL),
    VernierTree(NULL),
    cross_id(-1),
    trigraw(-1),
    trigscaled(-1),
    evtnumber(0),
    bbc_z(NAN),
    zdcll1_z(NAN)
{
  return ;
}

//================================================================================
//                         Init
//================================================================================
int VernierDSTReduction::Init(PHCompositeNode *topNode)
{  

  //=============== output file =====================================
  fout = new TFile(OutFileName.c_str(),"RECREATE", "Trig Scalars");

  //================== TTree =========================================

  VernierTree = new TTree("VernierTree","Information needed for BbcVertex efficiency + longtudinal bunch width");

  VernierTree->Branch("cross_id",     &cross_id,       "cross_id/I");
  VernierTree->Branch("evtnumber",    &evtnumber,      "evtnumber/I");
  VernierTree->Branch("trigraw",      &trigraw,        "trigraw/I");
  VernierTree->Branch("triglive",     &triglive,       "triglive/I");
  VernierTree->Branch("trigscaled",   &trigscaled,     "trigscaled/I");
  VernierTree->Branch("bbc_z",        &bbc_z,          "bbc_z/D");
  VernierTree->Branch("zdcll1_z",     &zdcll1_z,       "zdcll1_z/D");
  VernierTree->Branch("BbcNorthTubes",&BbcNorthTubes, "BbcNorthTubes/I");
  VernierTree->Branch("BbcSouthTubes",&BbcSouthTubes, "BbcSouthTubes/I");

  //====================== histograms ==================================
  sprintf (namestring, "BBC_novtxcut_Z");
  BBC_novtxcut_Z    = new TH1F(namestring, namestring, 1200, -300.0, 300.0);  
  sprintf (namestring, "BBC_narrow_Z");
  BBC_narrow_Z      = new TH1F(namestring, namestring, 1200, -300.0, 300.0);  
  sprintf (namestring, "ZDC_wide_Z");
  ZDC_wide_Z        = new TH1F(namestring, namestring, 1200, -300.0, 300.0);
  sprintf (namestring, "ZDC_narrow_Z");
  ZDC_narrow_Z      = new TH1F(namestring, namestring, 1200, -300.0, 300.0);
  sprintf (namestring, "BBC_wide_and_narrow_Z");
  BBC_wide_and_narrow_Z = new TH1F(namestring,namestring,1200,-300.0,300.0);
  sprintf (namestring, "BBC_wide_and_ZDC_wide_Z");
  BBC_wide_and_ZDC_wide_Z = new TH1F(namestring,namestring,1200,-300.0,300.0);

  // These are 'actually' histograms - i.e., with bins filled with "Fill" instead of "Set Bin Content".
  BBC_novtxcut_Z    ->Sumw2();
  BBC_narrow_Z      ->Sumw2();
  ZDC_wide_Z        ->Sumw2();
  ZDC_narrow_Z      ->Sumw2();
  BBC_wide_and_narrow_Z ->Sumw2();
  BBC_wide_and_ZDC_wide_Z ->Sumw2();


  BBC_RATE         = new TH1F("BBC_RATE",       "Crossing combined BBC rate",    asize, 0.0, historange);
  ZDCNarrow_RATE   = new TH1F("ZDCNarrow_RATE", "Crossing combined ZDCNarrow rate", asize, 0.0, historange);
  ZDCWide_RATE     = new TH1F("ZDCWide_RATE",   "Crossing combined ZDCWide rate",  asize, 0.0, historange);

  // These histograms are all filled with "set bin content" - so care has to be taken when
  // integrating or summing. Each bin will contain a number of trigger-firings, and therefore
  // there should not be any statistical weighting until bins are summed (with GetBinContent, 
  // ***NOT*** with "Integrate".).
  for(int k=0;k<120;k++)
  {
    sprintf (namestring, "BBC_RATE_%d",k);
    BBCRATE[k]          = new TH1F(namestring, namestring, histosize, 0.0, historange); 

    sprintf (namestring, "BBC_GL1P_%d",k);
    BBC_GL1P[k]             = new TH1F(namestring, namestring, histosize, 0.0, historange); 

    sprintf (namestring, "ZDCNarrow_RATE_%d",k);
    ZDCNarrowRATE[k]       = new TH1F(namestring, namestring, histosize, 0.0, historange);

    sprintf (namestring, "ZDCNarrow_GL1P_%d",k);
    ZDCNarrow_GL1P[k]          = new TH1F(namestring, namestring, histosize, 0.0, historange); 

    sprintf (namestring, "ZDCWide_RATE_%d",k);
    ZDCWideRATE[k]        = new TH1F(namestring, namestring, histosize, 0.0, historange);

    sprintf (namestring, "ZDCWide_GL1P_%d",k);
    ZDCWide_GL1P[k]           = new TH1F(namestring, namestring, histosize, 0.0, historange); 

    sprintf (namestring, "Clock_GL1P_%d",k);
    Clock_GL1P[k]           = new TH1F(namestring, namestring, histosize, 0.0, historange);    
  }

  BBCGL1P             = new TH1F("BBCGL1P",      "BBCGL1P",       histosize, 0.0, historange);
  ZDCWideGL1P         = new TH1F("ZDCWideGL1P",  "ZDCWideGL1P",   histosize, 0.0, historange);
  ZDCNarrowGL1P       = new TH1F("ZDCNarrowGL1P","ZDCNarrowGL1P", histosize, 0.0, historange);
  ClockGL1P           = new TH1F("ClockGL1P",    "ClockGL1P",     histosize, 0.0, historange);

  //================== intializing variable =============================
  ncalls           = 0;
  bbc_tot          = 0;
  ZDCNarrow_tot    = 0;
  ZDCWide_tot      = 0;

  for(int j = 0; j < asize; j++)
  {
    bbc_a_sum[j]       = 0.0;       
    clock_a_sum[j]     = 0.0;     
    ZDCNarrow_a_sum[j] = 0.0; 
    ZDCWide_a_sum[j]   = 0.0;   
    for(int l = 0; l < 120; l++)
    {
      bbc_a[l][j]       = 0.;
      clock_a[l][j]     = 0.;
      ZDCNarrow_a[l][j] = 0.;
      ZDCWide_a[l][j]   = 0.;  
    }
  }

  cout << " " << endl;
  cout << "Initializing variables." << endl;//debug
  cout << " " << endl;
  return 0;
}

//======================================================================
//                 Process Event
//======================================================================
int VernierDSTReduction::process_event(PHCompositeNode *topNode)
{   
  //============== counting no. of events ===============
  ncalls++;
  if((ncalls % 100000) == 0)
  {
    cout << "No of events processed:  = " << ncalls << endl;
  }

  // ============= reading required data nodes =========
  SpinDataEventOut *d_sde = getClass<SpinDataEventOutv2> (topNode,"SpinDataEventOut");
  if (!d_sde)
  {
    cout << PHWHERE << "Dude, what are you doing?  Your nodes are not in the tree." << endl;
    return 0;
  }
  TrigLvl1 *d_trig = getClass<TrigLvl1>(topNode,"TrigLvl1");
  if(!d_trig)
  {
    cout << PHWHERE << "TrigLevel1 Node not in the tree" << endl;
    return 0;
  }
  BbcOut *bbcout = getClass<BbcOut> (topNode, "BbcOut");
  if (!bbcout) 
  {
    cout<<"BbcOut does not exist"<<endl;
    return 0;
  }
  ZdcOut *zdcout = getClass<ZdcOut> (topNode, "ZdcOut");
  if (!zdcout) 
  {
    cout<<"ZdcOut does not exist"<<endl;
    return 0;
  }

  TriggerHelper* trigHelper = new TriggerHelper(topNode);
  bool trig_bbcnarrow	= trigHelper->didLevel1TriggerFire("BBCLL1(>0 tubes)");
  bool trig_bbcwide 	= trigHelper->didLevel1TriggerFire("BBCLL1(>0 tubes) novertex"); 
  bool trig_zdcnarrow	= trigHelper->didLevel1TriggerFire("ZDCLL1narrow");
  bool trig_zdcwide  	= trigHelper->didLevel1TriggerFire("ZDCLL1wide");
  delete trigHelper;

  // ======= getting event by event info =======================
  /*
     For RUN5 & RUN6:
     GL1P(board number, scaler type)
     type 0 = BBCLL1
     type 1 = Clock
     type 2 = ZDCLL1
     type 3 = ZDCNS
     For RUN8:
     type 0 = BBCLL1
     type 1 = Clock
     type 2 = ZDCNS
     type 3 = ZDCLL1narrow 
     For RUN12:
     type 0 =  BBCLL1(>0 tubes)
     type 1 =  CLOCK
     type 2 =  ZDCLL1wide
     type 3 =  ZDCLL1narrow

*/
  bbc_zvtx       = bbcout->get_VertexPoint();
  zdc_zvtx       = zdcout->get_Zvertex();
  eventsequence  = d_sde->GetEventSequence();
  BbcNorthTubes  = bbcout->get_nPmt(Bbc::North);
  BbcSouthTubes  = bbcout->get_nPmt(Bbc::South);

  gl1crossingID  = d_sde->GetGL1PCrossingID(0);
  bbcll1         = d_sde->GetGL1PScalerCount(0,0);
  gl1clock       = d_sde->GetGL1PScalerCount(0,1);
  ZDCWide        = d_sde->GetGL1PScalerCount(0,2);
  ZDCNarrow      = d_sde->GetGL1PScalerCount(0,3);

  cross_id       = gl1crossingID;
  evtnumber      = eventsequence;
  bbc_z          = bbc_zvtx;
  zdcll1_z       = zdc_zvtx;

  trigscaled     = d_trig->get_lvl1_trigscaled();
  triglive       = d_trig->get_lvl1_triglive();
  trigraw        = d_trig->get_lvl1_trigraw();

  //=======================================================
  if(gl1crossingID >= 120)
  {
    cout << "crossing ID out of bound: " << gl1crossingID << endl;
    return 0;
  }
  ievent = static_cast<int>(eventsequence/1000.0);
  clock_a[gl1crossingID][ievent]     += gl1clock;
  bbc_a[gl1crossingID][ievent]       += bbcll1;
  ZDCWide_a[gl1crossingID][ievent]   += ZDCWide;
  ZDCNarrow_a[gl1crossingID][ievent] += ZDCNarrow;
  clock_a_sum[ievent]     += gl1clock;      
  bbc_a_sum[ievent]       += bbcll1;          
  ZDCWide_a_sum[ievent]   += ZDCWide;     
  ZDCNarrow_a_sum[ievent] += ZDCNarrow; 

  bbc_tot       += bbcll1;
  ZDCNarrow_tot += ZDCNarrow;
  ZDCWide_tot   += ZDCWide;
  clock_tot     += gl1clock;

  //================== filling TTree ======================= 
  VernierTree->Fill(); 

  //======== filling histo's for vtx distribution for diff. triggers ==========
  // Assessing efficiency of LL1 trigger
  if(fabs(bbc_zvtx) <= 300.0)
  {
    if(trig_bbcwide) {
      BBC_novtxcut_Z->Fill(bbc_zvtx);
    }
    if(trig_bbcnarrow) {
      BBC_narrow_Z->Fill(bbc_zvtx);
    }
    if(trig_bbcwide && trig_bbcnarrow) {
      BBC_wide_and_narrow_Z->Fill(bbc_zvtx);
    }
  }
  if(fabs(zdc_zvtx) <= 300.0)
  {
    if(trig_zdcwide ) {
      ZDC_wide_Z->Fill(zdc_zvtx);
      if( trig_bbcwide ) {
        BBC_wide_and_ZDC_wide_Z->Fill(zdc_zvtx);
      }
    }
    if(trig_zdcnarrow) {
      ZDC_narrow_Z->Fill(zdc_zvtx);
    }
  }

  return 0;
}

//========================================================================================================
//                                    End
//========================================================================================================
int VernierDSTReduction::End(PHCompositeNode *topNode){
  // ============ Filling Histos =====================================
  for(int event_i = 0; event_i < asize; event_i++){ 
    //collecting the total rate summing over all bunch crossings 
    BBCGL1P->Fill(event_i, bbc_a_sum[event_i]);
    ZDCNarrowGL1P->Fill(event_i, ZDCNarrow_a_sum[event_i]);
    ZDCWideGL1P->Fill(event_i, ZDCWide_a_sum[event_i]);
    ClockGL1P->Fill(event_i,clock_a_sum[event_i]);
    if(clock_a_sum[event_i] > 0) {
      BBC_RATE->Fill(event_i, bbc_a_sum[event_i]/clock_a_sum[event_i]);
      ZDCNarrow_RATE->Fill(event_i, ZDCNarrow_a_sum[event_i]/clock_a_sum[event_i]);
      ZDCWide_RATE->Fill(event_i, ZDCWide_a_sum[event_i]/clock_a_sum[event_i]);
    }
    for(int bunch_i=0; bunch_i<120; bunch_i++){
      //=== filling histogram with scalers for each bunch crossing ID ===
      Clock_GL1P[bunch_i]->Fill(event_i, clock_a[bunch_i][event_i]);
      BBC_GL1P[bunch_i]->Fill(event_i, bbc_a[bunch_i][event_i]);
      ZDCNarrow_GL1P[bunch_i]->Fill(event_i, ZDCNarrow_a[bunch_i][event_i]);
      ZDCWide_GL1P[bunch_i]->Fill(event_i, ZDCWide_a[bunch_i][event_i]);

      if(clock_a[bunch_i][event_i] > 0) {
        //=== filling histogram with rates for each bunch crossing ID ===
        BBCRATE[bunch_i]->Fill(event_i, bbc_a[bunch_i][event_i]/clock_a[bunch_i][event_i]);
        ZDCNarrowRATE[bunch_i]->Fill(event_i, ZDCNarrow_a[bunch_i][event_i]/clock_a[bunch_i][event_i]);
        ZDCWideRATE[bunch_i]->Fill(event_i, ZDCWide_a[bunch_i][event_i]/clock_a[bunch_i][event_i]);
      }
    }
  } 


  fout->Write();  
  fout->Close();

  cout << " " << endl;
  cout << "Writing Output File." << endl;//debug
  cout << " " << endl;  
  return 0;
}
