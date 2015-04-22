#ifndef __VernierDSTReduction_H__
#define __VernierDSTReduction_H__

#include <SubsysReco.h>
#include <string>

class Fun4AllServer;
class Fun4AllHistoManager;
class PHCompositeNode;
class TH1;
class TFile;
class TTree;

class VernierDSTReduction: public SubsysReco
{
 public:
  VernierDSTReduction(const std::string &outfilename);
  virtual ~VernierDSTReduction() {}
 
  int  Init         (PHCompositeNode *topNode); 
  int  process_event(PHCompositeNode *topNode);
  int  End          (PHCompositeNode *topNode);

 protected:
  std::string OutFileName;
  TFile *fout;
  TFile *finpol;
  TTree *VernierTree;

  int cross_id;
  int trigraw;
  int triglive;
  int trigscaled;
  int evtnumber;
  double bbc_z;
  double zdcll1_z; 

  // asize needs to be events_in_run / 1000 at a minimum
  // daq=> select runnumber,eventsinrun,triggerconfig,runtype from run 
  // where triggerconfig like '%Run12_VS%' and runtype like 'VERNIERSCAN' 
  // order by eventsinrun desc limit 100;
  // 
  // runnumber | eventsinrun | triggerconfig |   runtype   
  //-----------+-------------+---------------+-------------
  //    364636 |    35164069 | PP510Run12_VS | VERNIERSCAN
  //    365866 |    31203288 | PP510Run12_VS | VERNIERSCAN
  //    366605 |    27547634 | PP510Run12_VS | VERNIERSCAN
  //    367138 |    24595822 | PP510Run12_VS | VERNIERSCAN
  //    362492 |    11951269 | PP200Run12_VS | VERNIERSCAN
  //    360877 |     9273002 | PP200Run12_VS | VERNIERSCAN
  //    360879 |     8724848 | PP200Run12_VS | VERNIERSCAN
  //    359711 |     8075075 | PP200Run12_VS | VERNIERSCAN
  //    367136 |     4992073 | PP510Run12_VS | VERNIERSCAN
  //    365682 |      465229 | PP510Run12_VS | VERNIERSCAN
  //    371329 |        3909 | PP200Run12_VS | VERNIERSCAN
  //    371331 |        3262 | PP200Run12_VS | VERNIERSCAN
  //    371328 |        1825 | PP200Run12_VS | VERNIERSCAN
  //    360878 |          67 | PP200Run12_VS | VERNIERSCAN
  // asize should be 40,000 just to be safe 
  static const int asize = 40000;
  static const int histosize;
  static const double historange;
  double bbc_a[120][asize];       // these arrays will sum GL1P counts 
  double clock_a[120][asize];     // subdivided by event index (summed over every
  double ZDCNarrow_a[120][asize]; // 1000 event_sequence) and GL1 crossing id.
  double ZDCWide_a[120][asize];   //
  double bbc_a_sum[asize];        // these arrays sum over the _a arrays' Gl1P index
  double clock_a_sum[asize];      // ..
  double ZDCNarrow_a_sum[asize];  // ..
  double ZDCWide_a_sum[asize];    // ..

  TH1 *BBC_GL1P[120];
  TH1 *ZDCNarrow_GL1P[120];
  TH1 *ZDCWide_GL1P[120];
  TH1 *Clock_GL1P[120];

  TH1 *BBCGL1P;
  TH1 *ClockGL1P;
  TH1 *ZDCWideGL1P;
  TH1 *ZDCNarrowGL1P;
  

  TH1 *BBCRATE[120];
  TH1 *ZDCNarrowRATE[120];
  TH1 *ZDCWideRATE[120];
  TH1 *BBC_RATE;
  TH1 *ZDCNarrow_RATE;
  TH1 *ZDCWide_RATE;

  TH1 *BBC_novtxcut_Z;
  TH1 *BBC_narrow_Z;
  TH1 *BBC_wide_and_narrow_Z;
  TH1 *BBC_wide_and_ZDC_wide_Z;
  TH1 *ZDC_wide_Z;
  TH1 *ZDC_narrow_Z;
  
  unsigned long int bbc_tot, ZDCNarrow_tot, ZDCWide_tot, clock_tot, BbcNorthTubes, BbcSouthTubes;
  double bbc_zvtx, zdc_zvtx;
  int gl1clock, bbcll1, ZDCNarrow, ZDCWide, eventsequence;
  int dohh, double_dohh;
  char namestring[64];
  long int ncalls;
  int ievent, gl1crossingID; 
  unsigned int scaled_trig;
};

#endif
