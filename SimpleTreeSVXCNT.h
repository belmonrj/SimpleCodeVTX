#ifndef __SIMPLETREESVXCNT__
#define __SIMPLETREESVXCNT__

#include <SubsysReco.h>
class PHCompositeNode;
class TFile;
class TH1;
class TNtuple;

class SimpleTreeSVXCNT : public SubsysReco
{

 public:

  SimpleTreeSVXCNT(const char *outfilename = "test.root") : SubsysReco("SimpleTreeSVXCNT"), d_outfilename(outfilename) {}
  virtual ~SimpleTreeSVXCNT(){}


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
  TH1 *th1f_cent;
  TH1 *th1f_vtxz;
  TH1 *th1f_bbcz;
  TH1 *th1f_nsvx;
  TH1 *th1f_ndch;
  TH1 *th1f_bbcsq;
  TH1 *th1f_bbcnq;


  TNtuple *ntpedch;
  TNtuple *ntpesvx;

};

#endif
