#include "TH1.h"
#include "TH2.h"
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












//added for svxcentral tracks

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
  // hm->setOutfileName(d_outfilename.c_str());
  //Fun4AllHistoManager *hm = se->getHistoManager("SimpleTreeSVXCNT");

  //-----------------------------
  //tim testing somet nutpple
  //-----------------------------
  d_OutputFile = new TFile(d_outfilename.c_str(),"RECREATE");




  ///create and register histograms here
  centrality = new TH1I("centrality","Event Centrality", 100,0,100);
  zvertex = new TH1F("zvertex", "Event Z-Vertex Position",100,-50,50);


  ////x-axis is variable, y-axis is momentum, z-axis is centrality



  ///elemcmatch: x-axis = dphi, y-axis = dz, z-axis = momentum



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

      }




  //added for loop over svxcentraltracks
  int nptbins = 20;
  int ptlowedge = 0;
  int pthighedge = 10;








  nsvx = new TH1F("nsvx",";number of svx trks per evnt;",30,0,10);
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

  if(abs(zvert) > 10) return 0;

  // cout << zvert << endl;
  zvertex->Fill(zvert);
  float numberdchtrks = 0;
  bool isElectronCandidate;
  for (unsigned int i=0; i<trk->get_npart();i++)
    {

      float mom = trk->get_mom(i);
      isElectronCandidate = 1;
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


      float phi = trk->get_phi(i);
      float alpha = trk->get_alpha(i);

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


      if(mom < min_mom_cut || mom > max_mom_cut) continue;

      
      float emc_t0 = trk->get_temc(i);
      float crk_t0 = trk->get_tcrk(i);
      float scrk_t0 = trk->get_stcrk(i);
     
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
	}
      else 
	{
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

