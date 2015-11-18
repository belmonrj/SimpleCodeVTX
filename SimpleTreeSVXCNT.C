#include <iostream>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "SimpleTreeSVXCNT.h"

#include "getClass.h"
#include "Fun4AllReturnCodes.h"
#include "RunHeader.h"
#include "PHGlobal.h"
#include "VtxOut.h"
#include "PHPoint.h" // needed for VtxOut.h ...
#include "TrigLvl1.h"
#include "PHCentralTrack.h"
#include "SvxCentralTrack.h"
#include "SvxCentralTrackList.h"


#define eID_nPbins 50
#define minPhi -0.7
#define maxPhi 3.8
#define minAlpha -1
#define maxAlpha 1
#define eID_minP 0.2
#define eID_maxP 15


using namespace std;



int SimpleTreeSVXCNT::Init(PHCompositeNode *topNode)
{

  d_nevent = 0;

  d_OutputFile = new TFile(d_outfilename.c_str(),"RECREATE");


  // --- tree
  tree_cnt = new TTree("tree_cnt","tree_cnt");
  // --- basic stuff
  tree_cnt->Branch("runn",&d_runn,"runn/I");
  tree_cnt->Branch("cent",&d_cent,"cent/S");
  tree_cnt->Branch("bbcz",&d_bbcz,"bbcz/F");
  tree_cnt->Branch("ntrk",&d_ntrk,"ntrk/I");
  tree_cnt->Branch("phi0",d_phi0,"phi0[ntrk]/F");
  tree_cnt->Branch("the0",d_the0,"the0[ntrk]/F");
  tree_cnt->Branch("quality",d_quality,"quality[ntrk]/S");
  tree_cnt->Branch("charge",d_charge,"charge[ntrk]/S");
  tree_cnt->Branch("dcarm",d_dcarm,"dcarm[ntrk]/S");
  tree_cnt->Branch("alpha",d_alpha,"alpha[ntrk]/F");
  tree_cnt->Branch("phi",d_phi,"phi[ntrk]/F");
  tree_cnt->Branch("zed",d_zed,"zed[ntrk]/F");
  tree_cnt->Branch("mom",d_mom,"mom[ntrk]/F");
  tree_cnt->Branch("pT",d_pT,"pT[ntrk]/F");
  // --- pc
  tree_cnt->Branch("ppc1x",d_ppc1x,"ppc1x[ntrk]/F");
  tree_cnt->Branch("ppc1y",d_ppc1y,"ppc1y[ntrk]/F");
  tree_cnt->Branch("ppc1z",d_ppc1z,"ppc1z[ntrk]/F");
  tree_cnt->Branch("ppc1phi",d_ppc1phi,"ppc1phi[ntrk]/F");
  tree_cnt->Branch("ppc2x",d_ppc2x,"ppc2x[ntrk]/F");
  tree_cnt->Branch("ppc2y",d_ppc2y,"ppc2y[ntrk]/F");
  tree_cnt->Branch("ppc2z",d_ppc2z,"ppc2z[ntrk]/F");
  tree_cnt->Branch("pc2z",d_pc2z,"pc2z[ntrk]/F");
  tree_cnt->Branch("pc2phi",d_pc2phi,"pc2phi[ntrk]/F");
  tree_cnt->Branch("pc2dz",d_pc2dz,"pc2dz[ntrk]/F");
  tree_cnt->Branch("pc2dphi",d_pc2dphi,"pc2dphi[ntrk]/F");
  tree_cnt->Branch("pc2sdz",d_pc2sdz,"pc2sdz[ntrk]/F");
  tree_cnt->Branch("pc2sdphi",d_pc2sdphi,"pc2sdphi[ntrk]/F");
  tree_cnt->Branch("ppc3x",d_ppc3x,"ppc3x[ntrk]/F");
  tree_cnt->Branch("ppc3y",d_ppc3y,"ppc3y[ntrk]/F");
  tree_cnt->Branch("ppc3z",d_ppc3z,"ppc3z[ntrk]/F");
  tree_cnt->Branch("pc3z",d_pc3z,"pc3z[ntrk]/F");
  tree_cnt->Branch("pc3phi",d_pc3phi,"pc3phi[ntrk]/F");
  tree_cnt->Branch("pc3dz",d_pc3dz,"pc3dz[ntrk]/F");
  tree_cnt->Branch("pc3dphi",d_pc3dphi,"pc3dphi[ntrk]/F");
  tree_cnt->Branch("pc3sdz",d_pc3sdz,"pc3sdz[ntrk]/F");
  tree_cnt->Branch("pc3sdphi",d_pc3sdphi,"pc3sdphi[ntrk]/F");
  // --- emc
  tree_cnt->Branch("deadmap",d_deadmap,"deadmap[ntrk]/I");
  tree_cnt->Branch("warnmap",d_warnmap,"warnmap[ntrk]/I");
  tree_cnt->Branch("pemcx",d_pemcx,"pemcx[ntrk]/F");
  tree_cnt->Branch("pemcy",d_pemcy,"pemcy[ntrk]/F");
  tree_cnt->Branch("pemcz",d_pemcz,"pemcz[ntrk]/F");
  tree_cnt->Branch("emcz",d_emcz,"emcz[ntrk]/F");
  tree_cnt->Branch("emcphi",d_emcdphi,"emcphi[ntrk]/F");
  tree_cnt->Branch("emcdz",d_emcdz,"emcdz[ntrk]/F");
  tree_cnt->Branch("emcdphi",d_emcdphi,"emcdphi[ntrk]/F");
  tree_cnt->Branch("emcsdz",d_emcsdz,"emcsdz[ntrk]/F");
  tree_cnt->Branch("emcsdphi",d_emcsdphi,"emcsdphi[ntrk]/F");
  tree_cnt->Branch("ecent",d_ecent,"ecent[ntrk]/F");
  tree_cnt->Branch("ecore",d_ecore,"ecore[ntrk]/F");
  tree_cnt->Branch("emce",d_emce,"emce[ntrk]/F");
  tree_cnt->Branch("prob",d_prob,"prob[ntrk]/F");
  tree_cnt->Branch("emcchi2",d_emcchi2,"emcchi2[ntrk]/F");
  tree_cnt->Branch("sect",d_sect,"sect[ntrk]/F");
  tree_cnt->Branch("ysect",d_ysect,"ysect[ntrk]/F");
  tree_cnt->Branch("zsect",d_zsect,"zsect[ntrk]/F");
  // --- rich
  tree_cnt->Branch("disp",d_disp,"disp[ntrk]/S");
  tree_cnt->Branch("n0",d_n0,"n0[ntrk]/S");
  tree_cnt->Branch("n1",d_n1,"n1[ntrk]/S");
  tree_cnt->Branch("n2",d_n2,"n2[ntrk]/S");
  tree_cnt->Branch("npe0",d_npe0,"npe0[ntrk]/F");
  tree_cnt->Branch("npe1",d_npe1,"npe1[ntrk]/F");
  tree_cnt->Branch("npe2",d_npe2,"npe2[ntrk]/F");


  // --- tree
  tree_svxcnt = new TTree("tree_svxcnt","tree_svxcnt");
  // --- basic stuff
  tree_svxcnt->Branch("runn",&d_runn,"runn/I");
  tree_svxcnt->Branch("cent",&d_cent,"cent/S");
  tree_svxcnt->Branch("bbcz",&d_bbcz,"bbcz/F");
  tree_svxcnt->Branch("ntrk",&d_ntrk,"ntrk/I");
  tree_svxcnt->Branch("phi0",d_phi0,"phi0[ntrk]/F");
  tree_svxcnt->Branch("the0",d_the0,"the0[ntrk]/F");
  tree_svxcnt->Branch("quality",d_quality,"quality[ntrk]/S");
  tree_svxcnt->Branch("charge",d_charge,"charge[ntrk]/S");
  tree_svxcnt->Branch("dcarm",d_dcarm,"dcarm[ntrk]/S");
  tree_svxcnt->Branch("alpha",d_alpha,"alpha[ntrk]/F");
  tree_svxcnt->Branch("phi",d_phi,"phi[ntrk]/F");
  tree_svxcnt->Branch("zed",d_zed,"zed[ntrk]/F");
  tree_svxcnt->Branch("mom",d_mom,"mom[ntrk]/F");
  tree_svxcnt->Branch("pT",d_pT,"pT[ntrk]/F");
  // ---
  tree_svxcnt->Branch("idch",d_idch,"idch[ntrk]/F");
  tree_svxcnt->Branch("chi2",d_chi2,"chi2[ntrk]/F");
  tree_svxcnt->Branch("ndf",d_ndf,"ndf[ntrk]/F");
  tree_svxcnt->Branch("dcat",d_dcat,"dcat[ntrk]/F");
  tree_svxcnt->Branch("dcal",d_dcal,"dcal[ntrk]/F");
  // --- pc
  tree_svxcnt->Branch("ppc1x",d_ppc1x,"ppc1x[ntrk]/F");
  tree_svxcnt->Branch("ppc1y",d_ppc1y,"ppc1y[ntrk]/F");
  tree_svxcnt->Branch("ppc1z",d_ppc1z,"ppc1z[ntrk]/F");
  tree_svxcnt->Branch("ppc1phi",d_ppc1phi,"ppc1phi[ntrk]/F");
  tree_svxcnt->Branch("ppc2x",d_ppc2x,"ppc2x[ntrk]/F");
  tree_svxcnt->Branch("ppc2y",d_ppc2y,"ppc2y[ntrk]/F");
  tree_svxcnt->Branch("ppc2z",d_ppc2z,"ppc2z[ntrk]/F");
  tree_svxcnt->Branch("pc2z",d_pc2z,"pc2z[ntrk]/F");
  tree_svxcnt->Branch("pc2phi",d_pc2phi,"pc2phi[ntrk]/F");
  tree_svxcnt->Branch("pc2dz",d_pc2dz,"pc2dz[ntrk]/F");
  tree_svxcnt->Branch("pc2dphi",d_pc2dphi,"pc2dphi[ntrk]/F");
  tree_svxcnt->Branch("pc2sdz",d_pc2sdz,"pc2sdz[ntrk]/F");
  tree_svxcnt->Branch("pc2sdphi",d_pc2sdphi,"pc2sdphi[ntrk]/F");
  tree_svxcnt->Branch("ppc3x",d_ppc3x,"ppc3x[ntrk]/F");
  tree_svxcnt->Branch("ppc3y",d_ppc3y,"ppc3y[ntrk]/F");
  tree_svxcnt->Branch("ppc3z",d_ppc3z,"ppc3z[ntrk]/F");
  tree_svxcnt->Branch("pc3z",d_pc3z,"pc3z[ntrk]/F");
  tree_svxcnt->Branch("pc3phi",d_pc3phi,"pc3phi[ntrk]/F");
  tree_svxcnt->Branch("pc3dz",d_pc3dz,"pc3dz[ntrk]/F");
  tree_svxcnt->Branch("pc3dphi",d_pc3dphi,"pc3dphi[ntrk]/F");
  tree_svxcnt->Branch("pc3sdz",d_pc3sdz,"pc3sdz[ntrk]/F");
  tree_svxcnt->Branch("pc3sdphi",d_pc3sdphi,"pc3sdphi[ntrk]/F");
  // --- emc
  tree_svxcnt->Branch("deadmap",d_deadmap,"deadmap[ntrk]/I");
  tree_svxcnt->Branch("warnmap",d_warnmap,"warnmap[ntrk]/I");
  tree_svxcnt->Branch("pemcx",d_pemcx,"pemcx[ntrk]/F");
  tree_svxcnt->Branch("pemcy",d_pemcy,"pemcy[ntrk]/F");
  tree_svxcnt->Branch("pemcz",d_pemcz,"pemcz[ntrk]/F");
  tree_svxcnt->Branch("emcz",d_emcz,"emcz[ntrk]/F");
  tree_svxcnt->Branch("emcphi",d_emcdphi,"emcphi[ntrk]/F");
  tree_svxcnt->Branch("emcdz",d_emcdz,"emcdz[ntrk]/F");
  tree_svxcnt->Branch("emcdphi",d_emcdphi,"emcdphi[ntrk]/F");
  tree_svxcnt->Branch("emcsdz",d_emcsdz,"emcsdz[ntrk]/F");
  tree_svxcnt->Branch("emcsdphi",d_emcsdphi,"emcsdphi[ntrk]/F");
  tree_svxcnt->Branch("ecent",d_ecent,"ecent[ntrk]/F");
  tree_svxcnt->Branch("ecore",d_ecore,"ecore[ntrk]/F");
  tree_svxcnt->Branch("emce",d_emce,"emce[ntrk]/F");
  tree_svxcnt->Branch("prob",d_prob,"prob[ntrk]/F");
  tree_svxcnt->Branch("emcchi2",d_emcchi2,"emcchi2[ntrk]/F");
  tree_svxcnt->Branch("sect",d_sect,"sect[ntrk]/F");
  tree_svxcnt->Branch("ysect",d_ysect,"ysect[ntrk]/F");
  tree_svxcnt->Branch("zsect",d_zsect,"zsect[ntrk]/F");
  // --- rich
  tree_svxcnt->Branch("disp",d_disp,"disp[ntrk]/S");
  tree_svxcnt->Branch("n0",d_n0,"n0[ntrk]/S");
  tree_svxcnt->Branch("n1",d_n1,"n1[ntrk]/S");
  tree_svxcnt->Branch("n2",d_n2,"n2[ntrk]/S");
  tree_svxcnt->Branch("npe0",d_npe0,"npe0[ntrk]/F");
  tree_svxcnt->Branch("npe1",d_npe1,"npe1[ntrk]/F");
  tree_svxcnt->Branch("npe2",d_npe2,"npe2[ntrk]/F");


  // --- create and register histograms here
  th1f_cent = new TH1I("th1f_cent","", 100,0,100);
  th1f_bbcz = new TH1F("th1f_bbcz","",100,-50,50);
  th1f_vtxz = new TH1F("th1f_vtxz","",100,-50,50);

  th1f_bbcsq = new TH1F("th1f_bbcsq","",3000,0,3000);
  th1f_bbcnq = new TH1F("th1f_bbcnq","",3000,0,3000);

  th1f_nsvx = new TH1F("th1f_nsvx","",20,0,20);
  th1f_ndch = new TH1F("th1f_ndch","",20,0,20);

  return 0;

}



int SimpleTreeSVXCNT::InitRun(PHCompositeNode *topNode)
{
  RunHeader *runhdr = findNode::getClass<RunHeader>(topNode,"RunHeader");
  if(!runhdr)
    {
      cout<<"RunHeader class not in Node Tree"<<endl;
      return 0;
    }
  d_runn = runhdr->get_RunNumber();
  return 0;
}



int SimpleTreeSVXCNT::process_event(PHCompositeNode *topNode)
{

  d_nevent++;
  if ( d_nevent % 1000 == 0 ) cout << "Events processed: " << d_nevent << endl;


  // --- get the PHGlobal node
  PHGlobal *global = findNode::getClass<PHGlobal>(topNode,"PHGlobal");
  if ( global == NULL )
    {
      cerr << PHWHERE << " ERROR: cannot find PHGlobal, discarding event" << endl;
      return DISCARDEVENT;
    }
  float bbcsq = global->getBbcChargeN();
  float bbcnq = global->getBbcChargeS();
  float cent = global->getCentrality(); // should be -9999 for p+p
  float bbcz = global->getBbcZVertex(); // vertex from bbc
  th1f_cent->Fill(cent);
  th1f_bbcz->Fill(bbcz);
  th1f_bbcsq->Fill(bbcsq);
  th1f_bbcnq->Fill(bbcnq);


  // --- get the TrigLvl1 node
  TrigLvl1 *lvl1 = findNode::getClass<TrigLvl1>(topNode,"TrigLvl1");
  if ( lvl1 == NULL )
    {
      cerr << "ERROR: TrigLvl1 node not found, discarding event" << endl;
      return DISCARDEVENT;
    }
  int trigbitscaled = lvl1->get_lvl1_trigscaled();
  //int trigbitraw = lvl1->get_lvl1_trigraw();


  // --- get the VtxOut node
  VtxOut *vtxout = findNode::getClass<VtxOut>(topNode,"VtxOut");
  if ( vtxout == NULL )
    {
      cerr << "VtxOut node not found." << endl;
      return DISCARDEVENT;
    }
  if ( (string)vtxout->which_Vtx() != "SVX_PRECISE" ) return DISCARDEVENT;
  float vtxz = (vtxout->get_Vertex()).getZ();
  if ( abs(vtxz) > 10 ) return DISCARDEVENT; // abs is overloaded in C++ (must use fabs in C)
  th1f_vtxz->Fill(vtxz);


  // --- get the PHCentralTrack node
  PHCentralTrack *node_cnt = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if ( node_cnt == NULL )
    {
      cerr << PHWHERE << "ERROR: can't find PHCentralTrack node, discarding event" << endl;
      return DISCARDEVENT;
    }


  // --- get the SvxCentralTrackList node
  SvxCentralTrackList *node_svxcnt = findNode::getClass<SvxCentralTrackList>(topNode,"SvxCentralTrackList");
  if(!node_svxcnt)
    {
      cerr << PHWHERE << " ERROR: can't find SVXcentralTrack, discarding event " << endl;
      return DISCARDEVENT;
    }


  // --- loop over Central Tracks
  int ntrk = 0;
  int numberdchtrks = 0;
  for ( unsigned int itrk = 0; itrk < node_cnt->get_npart(); itrk++)
    {

      // global position
      d_phi0[ntrk] = node_cnt->get_phi0(itrk);
      d_the0[ntrk] = node_cnt->get_the0(itrk);

      // dc
      d_quality[ntrk] = node_cnt->get_quality(itrk);
      d_charge[ntrk] = node_cnt->get_charge(itrk);
      d_dcarm[ntrk] = node_cnt->get_dcarm(itrk);
      d_alpha[ntrk] = node_cnt->get_alpha(itrk);
      d_phi[ntrk] = node_cnt->get_phi(itrk);
      d_zed[ntrk] = node_cnt->get_zed(itrk);
      d_mom[ntrk] = node_cnt->get_mom(itrk);
      d_pT[ntrk] = node_cnt->get_mom(itrk)*sin(node_cnt->get_the0(itrk));

      // pc
      d_ppc1x[ntrk] = node_cnt->get_ppc1x(itrk);
      d_ppc1y[ntrk] = node_cnt->get_ppc1y(itrk);
      d_ppc1z[ntrk] = node_cnt->get_ppc1z(itrk);
      d_ppc1phi[ntrk] = atan2(node_cnt->get_ppc1y(itrk),node_cnt->get_ppc1x(itrk));
      d_ppc2x[ntrk] = node_cnt->get_ppc2x(itrk);
      d_ppc2y[ntrk] = node_cnt->get_ppc2y(itrk);
      d_ppc2z[ntrk] = node_cnt->get_ppc2z(itrk);
      d_pc2z[ntrk] = node_cnt->get_ppc2z(itrk) - node_cnt->get_pc2dz(itrk);
      d_pc2phi[ntrk] = atan2(node_cnt->get_ppc2y(itrk),node_cnt->get_ppc2x(itrk)) - node_cnt->get_pc2dphi(itrk);
      d_pc2dz[ntrk] = node_cnt->get_pc2dz(itrk);
      d_pc2dphi[ntrk] = node_cnt->get_pc2dphi(itrk);
      d_pc2sdz[ntrk] = node_cnt->get_pc2sdz(itrk);
      d_pc2sdphi[ntrk] = node_cnt->get_pc2sdphi(itrk);
      d_ppc3x[ntrk] = node_cnt->get_ppc3x(itrk);
      d_ppc3y[ntrk] = node_cnt->get_ppc3y(itrk);
      d_ppc3z[ntrk] = node_cnt->get_ppc3z(itrk);
      d_pc3z[ntrk] = node_cnt->get_ppc3z(itrk) - node_cnt->get_pc3dz(itrk);
      d_pc3phi[ntrk] = atan2(node_cnt->get_ppc3y(itrk),node_cnt->get_ppc3x(itrk)) - node_cnt->get_pc3dphi(itrk);
      d_pc3dz[ntrk] = node_cnt->get_pc3dz(itrk);
      d_pc3dphi[ntrk] = node_cnt->get_pc3dphi(itrk);
      d_pc3sdz[ntrk] = node_cnt->get_pc3sdz(itrk);
      d_pc3sdphi[ntrk] = node_cnt->get_pc3sdphi(itrk);

      // emc
      d_deadmap[ntrk] = node_cnt->get_deadmap(itrk);
      d_warnmap[ntrk] = node_cnt->get_warnmap(itrk);
      d_pemcx[ntrk] = node_cnt->get_pemcx(itrk);
      d_pemcy[ntrk] = node_cnt->get_pemcy(itrk);
      d_pemcz[ntrk] = node_cnt->get_pemcz(itrk);
      d_emcz[ntrk] = node_cnt->get_pemcz(itrk) - node_cnt->get_emcdz(itrk);
      d_emcphi[ntrk] = atan2(node_cnt->get_pemcy(itrk),node_cnt->get_pemcx(itrk)) - node_cnt->get_emcdphi(itrk);
      d_emcdz[ntrk] = node_cnt->get_emcdz(itrk);
      d_emcdphi[ntrk] = node_cnt->get_emcdphi(itrk);
      d_emcsdz[ntrk] = node_cnt->get_emcsdz(itrk);
      d_emcsdphi[ntrk] = node_cnt->get_emcsdphi(itrk);
      d_ecore[ntrk] = node_cnt->get_ecore(itrk);
      d_ecent[ntrk] = node_cnt->get_ecent(itrk);
      d_emce[ntrk] = node_cnt->get_emce(itrk);
      d_prob[ntrk] = node_cnt->get_prob(itrk);
      d_emcchi2[ntrk] = node_cnt->get_emcchi2(itrk);
      d_sect[ntrk] = node_cnt->get_sect(itrk);
      d_ysect[ntrk] = node_cnt->get_ysect(itrk);
      d_zsect[ntrk] = node_cnt->get_zsect(itrk);

      // rich
      d_disp[ntrk] = node_cnt->get_disp(itrk);
      d_n0[ntrk] = node_cnt->get_n0(itrk);
      d_n1[ntrk] = node_cnt->get_n1(itrk);
      d_n2[ntrk] = node_cnt->get_n2(itrk);
      d_npe0[ntrk] = node_cnt->get_npe0(itrk);
      d_npe1[ntrk] = node_cnt->get_npe1(itrk);
      d_npe2[ntrk] = node_cnt->get_npe2(itrk);

      // --- if all track cuts passed, increment counter
      ntrk++;

    } // loop over Central Tracks

  d_ntrk = ntrk; // assign tree variable of number of tracks to counted variable
  if ( ntrk == 0 ) return DISCARDEVENT;
  tree_cnt->Fill();



  // --- trigger selection
  if ((trigbitscaled&0x00000010) == 16)
    {
      th1f_ndch->Fill(numberdchtrks);
    }


  ntrk = 0;
  int numbersvxtrks = 0;
  for ( int i = 0; i < node_svxcnt->get_nCentralTracks(); i++ )
    {

      SvxCentralTrack *sngl_svxcnt = node_svxcnt->getCentralTrack(i);

      if (!sngl_svxcnt) {cout << " no central track" << endl; continue;}
      float nhit = sngl_svxcnt->getNhits();
      if ( nhit <= 2 ) continue; // require 3 hit tracks...?

      int idch = sngl_svxcnt->getDchIndex();

      d_idch[ntrk] = idch;
      d_chi2[ntrk] = sngl_svxcnt->getChiSquare();
      d_ndf[ntrk] = sngl_svxcnt->getNDF();
      d_dcat[ntrk] = sngl_svxcnt->getDCA2D();
      d_dcal[ntrk] = sngl_svxcnt->getDCAZ();



      // global position
      d_phi0[ntrk] = node_cnt->get_phi0(idch);
      d_the0[ntrk] = node_cnt->get_the0(idch);

      // dc
      d_quality[ntrk] = node_cnt->get_quality(idch);
      d_charge[ntrk] = node_cnt->get_charge(idch);
      d_dcarm[ntrk] = node_cnt->get_dcarm(idch);
      d_alpha[ntrk] = node_cnt->get_alpha(idch);
      d_phi[ntrk] = node_cnt->get_phi(idch);
      d_zed[ntrk] = node_cnt->get_zed(idch);
      d_mom[ntrk] = node_cnt->get_mom(idch);
      d_pT[ntrk] = node_cnt->get_mom(idch)*sin(node_cnt->get_the0(idch));

      // pc
      d_ppc1x[ntrk] = node_cnt->get_ppc1x(idch);
      d_ppc1y[ntrk] = node_cnt->get_ppc1y(idch);
      d_ppc1z[ntrk] = node_cnt->get_ppc1z(idch);
      d_ppc1phi[ntrk] = atan2(node_cnt->get_ppc1y(idch),node_cnt->get_ppc1x(idch));
      d_ppc2x[ntrk] = node_cnt->get_ppc2x(idch);
      d_ppc2y[ntrk] = node_cnt->get_ppc2y(idch);
      d_ppc2z[ntrk] = node_cnt->get_ppc2z(idch);
      d_pc2z[ntrk] = node_cnt->get_ppc2z(idch) - node_cnt->get_pc2dz(idch);
      d_pc2phi[ntrk] = atan2(node_cnt->get_ppc2y(idch),node_cnt->get_ppc2x(idch)) - node_cnt->get_pc2dphi(idch);
      d_pc2dz[ntrk] = node_cnt->get_pc2dz(idch);
      d_pc2dphi[ntrk] = node_cnt->get_pc2dphi(idch);
      d_pc2sdz[ntrk] = node_cnt->get_pc2sdz(idch);
      d_pc2sdphi[ntrk] = node_cnt->get_pc2sdphi(idch);
      d_ppc3x[ntrk] = node_cnt->get_ppc3x(idch);
      d_ppc3y[ntrk] = node_cnt->get_ppc3y(idch);
      d_ppc3z[ntrk] = node_cnt->get_ppc3z(idch);
      d_pc3z[ntrk] = node_cnt->get_ppc3z(idch) - node_cnt->get_pc3dz(idch);
      d_pc3phi[ntrk] = atan2(node_cnt->get_ppc3y(idch),node_cnt->get_ppc3x(idch)) - node_cnt->get_pc3dphi(idch);
      d_pc3dz[ntrk] = node_cnt->get_pc3dz(idch);
      d_pc3dphi[ntrk] = node_cnt->get_pc3dphi(idch);
      d_pc3sdz[ntrk] = node_cnt->get_pc3sdz(idch);
      d_pc3sdphi[ntrk] = node_cnt->get_pc3sdphi(idch);

      // emc
      d_deadmap[ntrk] = node_cnt->get_deadmap(idch);
      d_warnmap[ntrk] = node_cnt->get_warnmap(idch);
      d_pemcx[ntrk] = node_cnt->get_pemcx(idch);
      d_pemcy[ntrk] = node_cnt->get_pemcy(idch);
      d_pemcz[ntrk] = node_cnt->get_pemcz(idch);
      d_emcz[ntrk] = node_cnt->get_pemcz(idch) - node_cnt->get_emcdz(idch);
      d_emcphi[ntrk] = atan2(node_cnt->get_pemcy(idch),node_cnt->get_pemcx(idch)) - node_cnt->get_emcdphi(idch);
      d_emcdz[ntrk] = node_cnt->get_emcdz(idch);
      d_emcdphi[ntrk] = node_cnt->get_emcdphi(idch);
      d_emcsdz[ntrk] = node_cnt->get_emcsdz(idch);
      d_emcsdphi[ntrk] = node_cnt->get_emcsdphi(idch);
      d_ecore[ntrk] = node_cnt->get_ecore(idch);
      d_ecent[ntrk] = node_cnt->get_ecent(idch);
      d_emce[ntrk] = node_cnt->get_emce(idch);
      d_prob[ntrk] = node_cnt->get_prob(idch);
      d_emcchi2[ntrk] = node_cnt->get_emcchi2(idch);
      d_sect[ntrk] = node_cnt->get_sect(idch);
      d_ysect[ntrk] = node_cnt->get_ysect(idch);
      d_zsect[ntrk] = node_cnt->get_zsect(idch);

      // rich
      d_disp[ntrk] = node_cnt->get_disp(idch);
      d_n0[ntrk] = node_cnt->get_n0(idch);
      d_n1[ntrk] = node_cnt->get_n1(idch);
      d_n2[ntrk] = node_cnt->get_n2(idch);
      d_npe0[ntrk] = node_cnt->get_npe0(idch);
      d_npe1[ntrk] = node_cnt->get_npe1(idch);
      d_npe2[ntrk] = node_cnt->get_npe2(idch);

      ntrk++;
    } // loop over SVX Central Tracks
  if ( ntrk == 0 ) return DISCARDEVENT;
  tree_svxcnt->Fill();



  // --- trigger selection
  if ((trigbitscaled&0x00000010) == 16)
    {
      th1f_nsvx->Fill(numbersvxtrks);
    }

  return EVENT_OK;

}



int SimpleTreeSVXCNT::End(PHCompositeNode *topNode)
{
    cout << "Writing out..." << endl;
    d_OutputFile->Write();
    cout << "Closing output file..." << endl;
    d_OutputFile->Close();
    delete d_OutputFile;

    return 0;
}

