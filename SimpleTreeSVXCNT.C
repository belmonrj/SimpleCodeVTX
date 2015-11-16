#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include <TFile.h>

#include "SimpleTreeSVXCNT.h"
#include "Fun4AllServer.h"
#include "Fun4AllHistoManager.h"
#include "PHCentralTrack.h"
#include "PHGlobal.h"
#include "getClass.h"
#include "TMath.h"
#include "RunHeader.h"
#include "Fun4AllReturnCodes.h"

#include "SvxCentralTrack.h"
#include "SvxCentralTrackReco.h"
#include "SvxCentralTrackList.h"
#include "VtxOut.h"
#include "PHPoint.h"
#include "TrigLvl1.h"


#define eID_nPbins 50
#define minPhi -0.7
#define maxPhi 3.8
#define minAlpha -1
#define maxAlpha 1 
#define eID_minP 0.2
#define eID_maxP 15


using namespace std;

int d_nevent;

TH1 *centrality = NULL;
TH1 *zvertex = NULL;

TH2 *elcut = NULL;

TH3 *eldep = NULL;
TH3 *elsdep = NULL;

TH3 *eln0 = NULL;
TH3 *elsn0 = NULL;

// TH3 *elchi2npe0 = NULL;
// TH3 *elschi2npe0 = NULL;

TH3 *elemcdphimatch = NULL;
TH3 *elemcdphimatch_e = NULL;
TH3 *elemcdphimatch_se = NULL;

TH3 *elemcdzmatch = NULL;
TH3 *elemcdzmatch_e = NULL;
TH3 *elemcdzmatch_se = NULL;

TH3 *elemcmatch = NULL;
TH3 *elemcmatch_e = NULL;
TH3 *elemcmatch_se = NULL;

TH2 *elacceptance = NULL;
TH2 *elacceptance_e = NULL;
TH2 *elacceptance_se = NULL;
TH3 *elep[8];
// TH3 *elsep[8];


TH2 *elmom2 = NULL;
//TH2 *elep2 =NULL;//for ghost sharing leaving out for now
TH3 *eldisp = NULL;
TH3 *elsdisp = NULL;

TH3 *elricht0 = NULL;
TH3 *elemct0 = NULL;
TH3 *elsricht0 = NULL;
TH3 *elsemct0 = NULL;

//added for svxcentral tracks
TH2 *svxchisqndfe = NULL;
TH2 *svxdcate = NULL;
TH2 *svxdcale = NULL;
TH2 *svxdepe = NULL;
TH2 *svxeope = NULL;
TH2 *svxn0e = NULL;
TH2 *svxdispe = NULL;
TH2 *svxchisqndfw = NULL;
TH2 *svxdcatw = NULL;
TH2 *svxdcalw = NULL;
TH2 *svxdepw = NULL;
TH2 *svxeopw = NULL;
TH2 *svxn0w = NULL;
TH2 *svxdispw = NULL;

TH1 *nsvx = NULL;
TH1 *ndch = NULL;

TNtuple* ntpedch = NULL;  
TNtuple* ntpesvx = NULL;  
TFile*   d_OutputFile = NULL;

int SimpleTreeSVXCNT::InitRun(PHCompositeNode *topNode)
{
  d_nevent=0;

  zvertex_cut = 10;
  //  emcmatch_cut = 3;
  min_mom_cut = 0.5;
  max_mom_cut = 15;
  emcdz_cut = 5;
  emcdphi_cut = 5;

  n0_cut = 1;
  // chi2npe0_cut = 50;
  disp_cut = 10;
  min_dep_cut = -10;
  max_dep_cut = 10;
  //  pfoa_cut = 2.5; //in degrees
  //  deltaZ_cut = 1.0;
  //  deltaPhi_cut = 0.1;


  // Fun4AllServer *se = Fun4AllServer::instance();
  // Fun4AllHistoManager *hm = new Fun4AllHistoManager("hists");
  // se->registerHistoManager(hm);
  // hm->setOutfileName(d_outfilename.c_str());
  //Fun4AllHistoManager *hm = se->getHistoManager("SimpleTreeSVXCNT");

  //-----------------------------
  //tim testing somet nutpple
  //-----------------------------
  d_OutputFile = new TFile(d_outfilename.c_str(),"RECREATE");




  ///create and register histograms here
  centrality = new TH1I("centrality","Event Centrality", 100,0,100);
  //hm->registerHisto(centrality);
  zvertex = new TH1F("zvertex", "Event Z-Vertex Position",100,-50,50);
  //hm->registerHisto(zvertex);

  elcut = new TH2F("elcut","Electron Cut", 15, -0.5, 14.5, 50, 0, 15);
  //hm->registerHisto(elcut);

  ////x-axis is variable, y-axis is momentum, z-axis is centrality
  eldep = new TH3F("eldep", "DEP", 200, -10, 10, eID_nPbins, min_mom_cut, max_mom_cut, 100, 0, 100);
  //hm->registerHisto(eldep);
  elsdep = new TH3F("elsdep", "DEP(swapped)", 200, -10, 10, eID_nPbins, min_mom_cut, max_mom_cut, 100, 0, 100);
  //hm->registerHisto(elsdep);
  eln0 = new TH3F("eln0", "RICH n0 ", 12, -2, 10, eID_nPbins, min_mom_cut, max_mom_cut, 100, 0, 100);
  //hm->registerHisto(eln0);
  elsn0 = new TH3F("elsn0", "RICH n0 (swapped)", 12, -2, 10, eID_nPbins, min_mom_cut, max_mom_cut, 100, 0, 100);
  //hm->registerHisto(elsn0);
  // elchi2npe0 = new TH3F("elchi2npe0", "chi2/npe0", 200, 0, 30, eID_nPbins, min_mom_cut, max_mom_cut, 100, 0, 100);
  // hm->registerHisto(elchi2npe0);
  // elschi2npe0 = new TH3F("elschi2npe0", "schi2/snpe0 (swapped)", 200, 0, 30, eID_nPbins, min_mom_cut, max_mom_cut, 100, 0, 100);
  // hm->registerHisto(elschi2npe0);

  elemcdphimatch = new TH3F("elemcdphimatch", "EMC dphi matching", 200,-10,10, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elemcdphimatch);
  elemcdphimatch_e = new TH3F("elemcdphimatch_e", "EMC dphi matching (electron candidate", 200,-10,10, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elemcdphimatch_e);
  elemcdphimatch_se = new TH3F("elemcdphimatch_se", "EMC dphi matching (swapped electron candidate)", 200,-10,10, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elemcdphimatch_se);

  elemcdzmatch = new TH3F("elemcdzmatch", "EMC dz matching", 200,-10,10, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elemcdzmatch);
  elemcdzmatch_e = new TH3F("elemcdzmatch_e", "EMC dz matching (electron candidate", 200,-10,10, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elemcdzmatch_e);
  elemcdzmatch_se = new TH3F("elemcdzmatch_se", "EMC dz matching (swapped electron candidate)", 200,-10,10, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elemcdzmatch_se);

  ///elemcmatch: x-axis = dphi, y-axis = dz, z-axis = momentum
  elemcmatch = new TH3F("elemcmatch", "EMC matching", 200,-10,10,500,-10,10, eID_nPbins, min_mom_cut, max_mom_cut);
  //hm->registerHisto(elemcmatch);
  elemcmatch_e = new TH3F("elemcmatch_e", "EMC matching (electron candidate", 200,-10,10,500,-10,10, eID_nPbins, min_mom_cut, max_mom_cut);
  //hm->registerHisto(elemcmatch_e);
  elemcmatch_se = new TH3F("elemcmatch_se", "EMC matching (swapped electron candidate)", 200,-10,10,500,-10,10, eID_nPbins, min_mom_cut, max_mom_cut);
  //hm->registerHisto(elemcmatch_se);


  elacceptance = new TH2F("elacceptance", "#phi x #alpha",100, minPhi, maxPhi, 100, minAlpha, maxAlpha);
  //hm->registerHisto(elacceptance);
  elacceptance_e = new TH2F("elacceptance_e", "#phi x #alpha",100, minPhi, maxPhi, 100, minAlpha, maxAlpha);
  //hm->registerHisto(elacceptance_e);
    elacceptance_se = new TH2F("elacceptance_se", "#phi x #alpha",100, minPhi, maxPhi, 100, minAlpha, maxAlpha);
  //hm->registerHisto(elacceptance_se);

  char name[10], title[40];
  for (int iarm = 0; iarm < 2; iarm++)
    for (int isector = 0; isector < 4; isector++)
      {
        int iarmsect = iarm * 4 + isector;
        sprintf(name, "elep_%d", iarmsect);
        if (iarm == 0)
          sprintf(title, "DEP x p x cent E%d", isector);
        if (iarm == 1)
          sprintf(title, "DEP x p x cent W%d", isector);
        elep[iarmsect] = new TH3F(name, title, 200, -10, 10, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
	//hm->registerHisto(elep[iarmsect]);

      }


  eldisp = new TH3F("eldisp", "RICH ring displacement", 40, 0, 10, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(eldisp);
  elsdisp = new TH3F("elsdisp", "RICH shared ring displacement", 40, 0, 10, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elsdisp);
  //elep2 = new TH2F("elep2", "#frac{E}{P}", 100, -5, 5,10,0,100);
  //hm->registerHisto(elep2);
  elmom2 = new TH2F("elmom2", "Electron momentum", 100, 0, 5,100,0,100);
  //hm->registerHisto(elmom2);
  elricht0 = new TH3F("elricht0", "RICH T0 (eID)", 100, -5, 5, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elricht0);
  elsricht0 = new TH3F("elsricht0", "RICH T0 (swapped eID)", 100, -5, 5, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elsricht0);
  elemct0 = new TH3F("elemct0", "EmCal T0 (eID)", 100, -40, 40, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elemct0);
  elsemct0 = new TH3F("elsemct0", "EmCal T0 (swapped eID)", 100, -40, 40, eID_nPbins, min_mom_cut, max_mom_cut,100,0,100);
  //hm->registerHisto(elsemct0);


  //added for loop over svxcentraltracks
  int nptbins = 20;
  int ptlowedge = 0;
  int pthighedge = 10;
  svxchisqndfe = new TH2F("svxchisqndfe","chisq/ndf",100,0,10,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxchisqndfe);
  svxchisqndfw = new TH2F("svxchisqndfw","chisq/ndf",100,0,10,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxchisqndfw);

  svxdcate = new TH2F("svxdcate","dcat",400,-0.4,0.4,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxdcate);
  svxdcatw = new TH2F("svxdcatw","dcat",400,-0.4,0.4,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxdcatw);

  svxdcale = new TH2F("svxdcale","dcal",400,-0.4,0.4,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxdcale);
  svxdcalw = new TH2F("svxdcalw","dcal",400,-0.4,0.4,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxdcalw);

  svxdepe = new TH2F("svxdepe","dep",200,-5,5,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxdepe);
  svxdepw = new TH2F("svxdepw","dep",200,-5,5,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxdepw);

  svxeope = new TH2F("svxeope","E/mom",200,0,5,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxeope);
  svxeopw = new TH2F("svxeopw","E/mom",200,0,5,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxeopw);

  svxn0e = new TH2F("svxn0e","n0",200,-5,5,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxn0e);
  svxn0w = new TH2F("svxn0w","n0",200,-5,5,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxn0w);

  svxdispe = new TH2F("svxdispe","n0",200,-10,10,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxdispe);
  svxdispw = new TH2F("svxdispw","n0",200,-10,10,nptbins,ptlowedge,pthighedge);
  //hm->registerHisto(svxdispw);


  nsvx = new TH1F("nsvx",";number of svx trks per evnt;",30,0,10);
  //hm->registerHisto(nsvx);
  ndch = new TH1F("ndch",";number of dch trks per evnt;",60,0,20);


  ntpedch = new TNtuple("ntpedch","","mom:emce:quality:phi:alpha:n0:sn0:disp:sdisp:dep:emcdphi:emcdz:pt:zvtx:scaledtrigbit:rawtrigbit");
  ntpesvx = new TNtuple("ntpesvx","","chisq:ndf:dcat:dcal:pt:disp:sdisp:dep:mom:emce:n0:sn0:quality:phi0:zvtx:scaledtrigbit:rawtrigbit:emcdphi:emcdz");

   
  

  return 0;
}

int SimpleTreeSVXCNT::process_event(PHCompositeNode *topNode)
{

  float ntpdch[99]; for(int i=0; i<99; i++) {ntpdch[i]=-9999.;} 
  float ntpsvx[99]; for(int i=0; i<99; i++) {ntpsvx[i]=-9999.;} 

  d_nevent++;
  if(d_nevent%10000==0)cout<<"Events processed: "<<d_nevent<<endl;
  
  PHCentralTrack *trk = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  PHGlobal *global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  SvxCentralTrackList *svxcnttrklist = findNode::getClass<SvxCentralTrackList>(topNode,"SvxCentralTrackList");
  if(!svxcnttrklist) {cerr <<  PHWHERE  << " ERROR: Can't find SVXcentralTrack "<< endl; return DISCARDEVENT;}
  VtxOut *vtxout = findNode::getClass<VtxOut>(topNode,"VtxOut");
  if(vtxout == NULL) { cerr << "VtxOut node not found." << endl; return DISCARDEVENT;}
  //float bbcQ = global->getBbcChargeN()+global->getBbcChargeS();
  float cent=global->getCentrality();
  
  TrigLvl1 *lvl1 = findNode::getClass<TrigLvl1>(topNode,"TrigLvl1");
  if(lvl1 == NULL) { cerr << "TrigLvl1 node not found" << endl; return DISCARDEVENT;}

  int trigbitscaled = lvl1->get_lvl1_trigscaled();
  int trigbitraw = lvl1->get_lvl1_trigraw();

  // cout << "raw bit " << trigbitraw << " scaled bit " << trigbitscaled << endl;

  string vertexsourcecode = vtxout->which_Vtx();
  if (vertexsourcecode != "SVX_PRECISE") {return DISCARDEVENT;}
  centrality->Fill(cent);

  // float zvert = global->getBbcZVertex();

  float zvert   =  (vtxout->get_Vertex()).getZ();

  elcut->Fill(0.,min_mom_cut);
  if(abs(zvert) > 10) return 0;

  // cout << zvert << endl;
  zvertex->Fill(zvert);
  float numberdchtrks = 0;
  bool isElectronCandidate;
  for (unsigned int i=0; i<trk->get_npart();i++)
    {

      float mom = trk->get_mom(i);
      isElectronCandidate = 1;
      if(isElectronCandidate) elcut->Fill(2,mom);
      float dchpx = trk->get_px(i);
      float dchpy = trk->get_py(i);
      float pt = sqrt(dchpx*dchpx+dchpy*dchpy);
      float emce = trk->get_emce(i);
      if (pt > 1 && pt < 4)
	{
	  numberdchtrks = numberdchtrks+1;
	}
      int quality = trk->get_quality(i);
      if(quality != 31 && quality != 63) continue;

      if(isElectronCandidate) elcut->Fill(3,mom);

      float phi = trk->get_phi(i);
      float alpha = trk->get_alpha(i);
      elacceptance->Fill(phi, alpha);

      int n0 = trk->get_n0(i);
      int sn0 = trk->get_sn0(i);

      // float chi2 = trk->get_chi2(i);
      // float npe0 = trk->get_npe0(i);
      // float chi2npe0 = chi2/npe0;
      // float schi2 = trk->get_schi2(i);
      // float snpe0 = trk->get_snpe0(i);
      // float schi2npe0 = schi2/snpe0;

      float disp = trk->get_disp(i);
      float sdisp = trk->get_sdisp(i);
      int dcarm = trk->get_dcarm(i);
      int sect = trk->get_sect(i);
      int iarmsect = dcarm * 4 + sect;

      //float E = trk->get_ecore(i);
      //float dep = E/mom-1;
      float dep = trk->get_dep(i);

      bool epcut = (dep > min_dep_cut && dep < max_dep_cut);

      float emcdphi = trk->get_emcsdphi_e(i);
      float emcdz = trk->get_emcsdz_e(i);
      bool emcmatchcut = (fabs(emcdphi)<emcdphi_cut && fabs(emcdz)<emcdz_cut);
      float zvtx   =  (vtxout->get_Vertex()).getZ();

      ntpdch[0] = mom;
      ntpdch[1] = emce;
      ntpdch[2] = quality;
      ntpdch[3] = phi;
      ntpdch[4] = alpha;
      ntpdch[5] = n0;
      ntpdch[6] = sn0;
      ntpdch[7] = disp;
      ntpdch[8] = sdisp;
      ntpdch[9] = dep;
      ntpdch[10] = emcdphi;
      ntpdch[11] = emcdz;
      ntpdch[12] = pt;
      ntpdch[13] = zvtx;
      ntpdch[14] = trigbitscaled;
      ntpdch[15] = trigbitraw;

      ntpedch->Fill(ntpdch);

      if(isElectronCandidate && abs(emcdz) < emcdz_cut) elcut->Fill(4,mom);
      elemcdzmatch->Fill(emcdz,mom,cent);
      if(n0>n0_cut && disp<disp_cut && epcut) elemcdzmatch_e->Fill(emcdz,mom,cent);
      if(sn0>n0_cut && sdisp<disp_cut && epcut) elemcdzmatch_se->Fill(emcdz,mom,cent);

      if(isElectronCandidate && abs(emcdphi) < emcdphi_cut) elcut->Fill(5,mom);
      elemcdphimatch->Fill(emcdphi,mom,cent);
      if(n0>n0_cut && disp<disp_cut && epcut) elemcdphimatch_e->Fill(emcdphi,mom,cent);
      if(sn0>n0_cut && sdisp<disp_cut && epcut) elemcdphimatch_se->Fill(emcdphi,mom,cent);

      elemcmatch->Fill(emcdphi,emcdz,mom);
      if(n0>n0_cut && disp<disp_cut && epcut) elemcmatch_e->Fill(emcdphi,emcdz,mom);
      if(sn0>n0_cut && sdisp<disp_cut && epcut) elemcmatch_se->Fill(emcdphi,emcdz,mom);

      if(!emcmatchcut) continue;
      if(isElectronCandidate) elcut->Fill(6,mom);

      if(mom < min_mom_cut || mom > max_mom_cut) continue;
      if(isElectronCandidate) elcut->Fill(7,mom);

      bool isSwappedElectronCandidate = isElectronCandidate;
      if(isElectronCandidate) eln0->Fill(n0,mom,cent);
      if(isSwappedElectronCandidate) elsn0->Fill(sn0,mom,cent);

      isElectronCandidate &= (n0>n0_cut);
      isSwappedElectronCandidate &= (sn0>n0_cut);

      if(isElectronCandidate) elcut->Fill(8,mom);

      // if(isElectronCandidate && epcut) elchi2npe0->Fill(chi2npe0,mom,cent);
      // if(isSwappedElectronCandidate && epcut) elschi2npe0->Fill(schi2npe0,mom,cent);

      // isElectronCandidate &= (chi2npe0 < chi2npe0_cut);
      // isSwappedElectronCandidate &= (schi2npe0 < chi2npe0_cut);

      if(isElectronCandidate) elcut->Fill(9,mom);


      if(isElectronCandidate && epcut) eldisp->Fill(disp,mom,cent);
      if(isSwappedElectronCandidate && epcut) elsdisp->Fill(sdisp,mom,cent);

      isElectronCandidate &= (disp < disp_cut);
      isSwappedElectronCandidate &= (sdisp < disp_cut);

      if(isElectronCandidate) elcut->Fill(10,mom);
      if(isElectronCandidate) eldep->Fill(dep,mom,cent);
      if(isSwappedElectronCandidate) elsdep->Fill(dep,mom,cent);
      if(isElectronCandidate) {elep[iarmsect]->Fill(dep,mom,cent);}
      // if(isSwappedElectronCandidate){ cout << dep << " ," << cent << endl; elsep[iarmsect]->Fill(dep,mom,cent);}


      isElectronCandidate &= (epcut);
      isSwappedElectronCandidate &= (epcut);

      if(isElectronCandidate) elcut->Fill(11,mom);
      
      if(isElectronCandidate) elmom2->Fill(mom,cent);

      if(isElectronCandidate) elacceptance_e->Fill(phi,alpha);
      if(isSwappedElectronCandidate) elacceptance_se->Fill(phi,alpha);
      
      float emc_t0 = trk->get_temc(i);
      float crk_t0 = trk->get_tcrk(i);
      float scrk_t0 = trk->get_stcrk(i);
      if(isElectronCandidate) elricht0->Fill(crk_t0,mom,cent);
      if(isElectronCandidate) elemct0->Fill(emc_t0,mom,cent);
      if(isSwappedElectronCandidate) elsricht0->Fill(scrk_t0,mom,cent);
      if(isSwappedElectronCandidate) elsemct0->Fill(emc_t0,mom,cent);      
     
    }
   if ((trigbitscaled&0x00000010) == 16)
    {
      ndch->Fill(numberdchtrks);
    }
  //count the number of SVXCentralTracks with nhit>2, and are good quality
  int nsvxtrks = 0;
  // int maxnumbertracks = svxcnttrklist->get_nCentralTracks();
  for (int i=0; i<svxcnttrklist->get_nCentralTracks();i++)
    {

      SvxCentralTrack* svxcntrltrk = svxcnttrklist->getCentralTrack(i);

      if (!svxcntrltrk) {cout << " no central track" << endl; continue;}

      int dchindex = svxcntrltrk->getDchIndex();
      float nhit = svxcntrltrk->getNhits();
      if (nhit<=2) {continue;}

      float chisq = svxcntrltrk->getChiSquare();
      float ndf = svxcntrltrk->getNDF();
      float dcat = svxcntrltrk->getDCA2D();
      float dcal = svxcntrltrk->getDCAZ();
      float dchpx = trk->get_px(dchindex);
      float dchpy = trk->get_py(dchindex);
      float pt = sqrt(dchpx*dchpx+dchpy*dchpy);
      float disp = trk->get_disp(dchindex);
      float sdisp = trk->get_sdisp(dchindex);
      float dep = trk->get_dep(dchindex);
      float emce = trk->get_emce(dchindex);
      float mom = trk->get_mom(dchindex);
      float n0 = trk->get_n0(dchindex);
      float emcdphi = trk->get_emcsdphi_e(dchindex);
      float emcdz = trk->get_emcsdz_e(dchindex);
      float sn0 = trk->get_sn0(dchindex);
      float dchquality = trk->get_quality(dchindex);
      float phi0 = trk->get_phi0(dchindex);
      float zvtx   =  (vtxout->get_Vertex()).getZ();
      if (dchquality != 31 && dchquality != 63) {continue;}

      ntpsvx[0] = chisq;
      ntpsvx[1] = ndf;
      ntpsvx[2] = dcat;
      ntpsvx[3] = dcal;
      ntpsvx[4] = pt;
      ntpsvx[5] = disp;
      ntpsvx[6] = sdisp;
      ntpsvx[7] = dep;
      ntpsvx[8] = mom;
      ntpsvx[9] = emce;
      ntpsvx[10] = n0;
      ntpsvx[11] = sn0;
      ntpsvx[12] = dchquality;
      ntpsvx[13] = phi0;
      ntpsvx[14] = zvtx;
      ntpsvx[15] = trigbitscaled;
      ntpsvx[16] = trigbitraw;
      ntpsvx[17] = emcdphi;
      ntpsvx[18] = emcdz;
      ntpesvx->Fill(ntpsvx);



      if (pt > 1 && pt < 4)
	{
      nsvxtrks = nsvxtrks+1;   //iterate the SVXCentralTrack counter
	}
      // float emcdphi = trk->get_emcsdphi_e(dchindex);
      // float emcdz = trk->get_emcsdz_e(dchindex);
      // bool emcmatchcut = (fabs(emcdphi)<emcdphi_cut && fabs(emcdz)<emcdz_cut);
      if (n0<=0)	{continue;}
      if (fabs(phi0) < 1.5)
	{
	  svxchisqndfe->Fill(chisq/ndf,pt);
	  svxdcate->Fill(dcat,pt);
	  svxdcale->Fill(dcal,pt);
	  svxdepe->Fill(dep,pt);
	  svxeope->Fill(emce/mom,pt);
	  svxn0e->Fill(n0,pt);
	  svxdispe->Fill(disp,pt);
	}
      else 
	{
	  svxchisqndfw->Fill(chisq/ndf,pt);
	  svxdcatw->Fill(dcat,pt);
	  svxdcalw->Fill(dcal,pt);
	  svxdepw->Fill(dep,pt);
	  svxeopw->Fill(emce/mom,pt);
	  svxn0w->Fill(n0,pt);
	  svxdispw->Fill(disp,pt);
	}
    }
  if ((trigbitscaled&0x00000010) == 16)
    {
  nsvx->Fill(nsvxtrks);
    }
  return 0;
}

int SimpleTreeSVXCNT::End(PHCompositeNode *topNode)
{
    // Fun4AllServer *se = Fun4AllServer::instance();
    // Fun4AllHistoManager *hm = se->getHistoManager("hists");
    // hm->dumpHistos();
    cout << "Writing out..." << endl;
    d_OutputFile->Write();
    cout << "Closing output file..." << endl;
    d_OutputFile->Close();
    delete d_OutputFile;
    
    return 0;
}

