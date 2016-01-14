#ifndef __SIMPLETREESVXCNT__
#define __SIMPLETREESVXCNT__

#include <SubsysReco.h>
class PHCompositeNode;
class TFile;
class TH1;
class TTree;

class SimpleTreeSVXCNT : public SubsysReco
{

 public:

  SimpleTreeSVXCNT(const char *outfilename = "test.root") : SubsysReco("SimpleTreeSVXCNT"), d_outfilename(outfilename) {}
  virtual ~SimpleTreeSVXCNT(){}


  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_zvertex_cut(float tmp){zvertex_cut = tmp;};
  void set_emcdz_cut(float tmp){emcdz_cut = tmp;};
  void set_emcdphi_cut(float tmp){emcdphi_cut = tmp;};
  void set_min_mom_cut(float tmp){min_mom_cut = tmp;};
  void set_max_mom_cut(float tmp){max_mom_cut = tmp;};
  void set_n0_cut(int tmp){n0_cut = tmp;};
  void set_chi2npe0_cut(float tmp){chi2npe0_cut = tmp;};
  void set_disp_cut(float tmp){disp_cut = tmp;};
  void set_min_dep_cut(float tmp){min_dep_cut = tmp;};
  void set_max_dep_cut(float tmp){max_dep_cut = tmp;};

 private:

  TFile *d_OutputFile;

  std::string d_outfilename;
  float zvertex_cut;
  float emcdz_cut;
  float emcdphi_cut;
  float min_mom_cut;
  float max_mom_cut;
  int n0_cut;
  float chi2npe0_cut;
  float disp_cut;
  float min_dep_cut;
  float max_dep_cut;

  int d_nevent;

  // --- conversion from base class to derived class allowed
  // --- but should consider being explicit here...
  TTree *tree_cnt;
  TTree *tree_svxcnt;

  TH1 *th1f_cent;
  TH1 *th1f_vtxz;
  TH1 *th1f_bbcz;
  TH1 *th1f_nsvx;
  TH1 *th1f_ndch;
  TH1 *th1f_bbcsq;
  TH1 *th1f_bbcnq;

  TH1 *th1f_dcat_B3_ladder[24];
  TH1 *th1f_dcat_B2_ladder[16];
  TH1 *th1f_dcat_B1_ladder[20];
  TH1 *th1f_dcat_B0_ladder[10];
  TH1 *th1f_dcal_B3_ladder[24];
  TH1 *th1f_dcal_B2_ladder[16];
  TH1 *th1f_dcal_B1_ladder[20];
  TH1 *th1f_dcal_B0_ladder[10];

  TH1 *th1f_posdcat_B3_ladder[24];
  TH1 *th1f_posdcat_B2_ladder[16];
  TH1 *th1f_posdcat_B1_ladder[20];
  TH1 *th1f_posdcat_B0_ladder[10];
  TH1 *th1f_posdcal_B3_ladder[24];
  TH1 *th1f_posdcal_B2_ladder[16];
  TH1 *th1f_posdcal_B1_ladder[20];
  TH1 *th1f_posdcal_B0_ladder[10];

  TH1 *th1f_negdcat_B3_ladder[24];
  TH1 *th1f_negdcat_B2_ladder[16];
  TH1 *th1f_negdcat_B1_ladder[20];
  TH1 *th1f_negdcat_B0_ladder[10];
  TH1 *th1f_negdcal_B3_ladder[24];
  TH1 *th1f_negdcal_B2_ladder[16];
  TH1 *th1f_negdcal_B1_ladder[20];
  TH1 *th1f_negdcal_B0_ladder[10];


  int d_runn;
  int d_cent;
  int d_ntrk;
  float d_bbcz;

  static const int maxn = 5000;

  //global position
  float d_phi0[maxn];
  float d_the0[maxn];

  //dc
  int d_quality[maxn];
  int d_charge[maxn];
  int d_dcarm[maxn];
  float d_alpha[maxn];
  float d_phi[maxn];
  float d_zed[maxn];
  float d_mom[maxn];
  float d_pT[maxn];

  //emc
  int d_deadmap[maxn];
  int d_warnmap[maxn];
  float d_emcphi[maxn];
  float d_emcz[maxn];
  float d_emcdphi[maxn];
  float d_emcdz[maxn];
  float d_emcsdphi[maxn];
  float d_emcsdz[maxn];
  float d_ecent[maxn];
  float d_ecore[maxn];
  float d_emce[maxn];
  float d_prob[maxn];
  float d_emcchi2[maxn];

  //rich
  float d_disp[maxn];
  int d_n0[maxn];
  float d_npe0[maxn];
  float d_crkchi2[maxn];

  // --- svxcnt
  int d_idch[maxn];
  float d_chi2[maxn];
  int d_ndf[maxn];
  float d_dcat[maxn];
  float d_dcal[maxn];




};

#endif
