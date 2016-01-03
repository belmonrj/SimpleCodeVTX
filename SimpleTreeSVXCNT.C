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
#include "SvxClusterInfo.h"


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
  // --- emc
  tree_cnt->Branch("deadmap",d_deadmap,"deadmap[ntrk]/I");
  tree_cnt->Branch("warnmap",d_warnmap,"warnmap[ntrk]/I");
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
  // --- rich
  tree_cnt->Branch("disp",d_disp,"disp[ntrk]/S");
  tree_cnt->Branch("n0",d_n0,"n0[ntrk]/S");
  tree_cnt->Branch("npe0",d_npe0,"npe0[ntrk]/F");
  tree_cnt->Branch("crkchi2",d_crkchi2,"crkchi2[ntrk]/F");


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
  // --- emc
  tree_svxcnt->Branch("deadmap",d_deadmap,"deadmap[ntrk]/I");
  tree_svxcnt->Branch("warnmap",d_warnmap,"warnmap[ntrk]/I");
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
  // --- rich
  tree_svxcnt->Branch("disp",d_disp,"disp[ntrk]/S");
  tree_svxcnt->Branch("n0",d_n0,"n0[ntrk]/S");
  tree_svxcnt->Branch("npe0",d_npe0,"npe0[ntrk]/F");
  tree_svxcnt->Branch("crkchi2",d_crkchi2,"crkchi2[ntrk]/F");


  // --- create and register histograms here
  th1f_cent = new TH1F("th1f_cent","", 100,0,100);
  th1f_bbcz = new TH1F("th1f_bbcz","",100,-50,50);
  th1f_vtxz = new TH1F("th1f_vtxz","",100,-50,50);

  th1f_bbcsq = new TH1F("th1f_bbcsq","",3000,0,3000);
  th1f_bbcnq = new TH1F("th1f_bbcnq","",3000,0,3000);

  th1f_nsvx = new TH1F("th1f_nsvx","",20,0,20);
  th1f_ndch = new TH1F("th1f_ndch","",20,0,20);

  // --- histogram arrays for ladder study...

  for(int i=0; i<24; i++)
    {
      th1f_dcat_B3_ladder[i] = new TH1F(Form("th1f_dcat_B3_ladder%d",i),"",200,-1,1);
      if(i<16) th1f_dcat_B2_ladder[i] = new TH1F(Form("th1f_dcat_B2_ladder%d",i),"",200,-1,1);
      if(i<20) th1f_dcat_B1_ladder[i] = new TH1F(Form("th1f_dcat_B1_ladder%d",i),"",200,-1,1);
      if(i<10) th1f_dcat_B0_ladder[i] = new TH1F(Form("th1f_dcat_B0_ladder%d",i),"",200,-1,1);
      // ---
      th1f_dcal_B3_ladder[i] = new TH1F(Form("th1f_dcal_B3_ladder%d",i),"",200,-1,1);
      if(i<16) th1f_dcal_B2_ladder[i] = new TH1F(Form("th1f_dcal_B2_ladder%d",i),"",200,-1,1);
      if(i<20) th1f_dcal_B1_ladder[i] = new TH1F(Form("th1f_dcal_B1_ladder%d",i),"",200,-1,1);
      if(i<10) th1f_dcal_B0_ladder[i] = new TH1F(Form("th1f_dcal_B0_ladder%d",i),"",200,-1,1);
    }

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

      // emc
      d_deadmap[ntrk] = node_cnt->get_deadmap(itrk);
      d_warnmap[ntrk] = node_cnt->get_warnmap(itrk);
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

      // rich
      d_disp[ntrk] = node_cnt->get_disp(itrk);
      d_n0[ntrk] = node_cnt->get_n0(itrk);
      d_npe0[ntrk] = node_cnt->get_npe0(itrk);
      d_crkchi2[ntrk] = node_cnt->get_chi2(itrk);

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
      if ( !sngl_svxcnt ) {cout << " no central track" << endl; continue;}
      float nhits = sngl_svxcnt->getNhits();
      d_chi2[ntrk] = sngl_svxcnt->getChiSquare();
      d_ndf[ntrk] = sngl_svxcnt->getNDF();
      d_dcat[ntrk] = sngl_svxcnt->getDCA2D();
      d_dcal[ntrk] = sngl_svxcnt->getDCAZ();


      float cx[4] = {-999,-999,-999,-999};
      float cy[4] = {-999,-999,-999,-999};
      float cz[4] = {-999,-999,-999,-999};
      float cr[4] = {-999,-999,-999,-999};
      float cphi[4] = {-999,-999,-999,-999};
      //cout << i << " " << nhits << endl;
      for ( int ihit = 0; ihit < nhits; ihit++ )
	{
	  SvxClusterInfo *cluster = (SvxClusterInfo *)sngl_svxcnt->getClusterInfo(ihit);
	  float x = cluster->getPosition(0);
	  float y = cluster->getPosition(1);
	  float z = cluster->getPosition(2);
	  //cout << ihit << " " << x << " " << y << " " << z << " " << endl;
	  cx[ihit] = x;
	  cy[ihit] = y;
	  cz[ihit] = z;
	  cr[ihit] = sqrt(x*x + y*y);
	  cphi[ihit] = atan2(x,y);
	  int layer = cluster->getLayer();
	  int ladder = cluster->getLadder();
	  int sensor = cluster->getSensor();
	  if(layer==0) th1f_dcat_B0_ladder[ladder]->Fill(d_dcat[ntrk]);
	  if(layer==1) th1f_dcat_B1_ladder[ladder]->Fill(d_dcat[ntrk]);
	  if(layer==2) th1f_dcat_B2_ladder[ladder]->Fill(d_dcat[ntrk]);
	  if(layer==3) th1f_dcat_B3_ladder[ladder]->Fill(d_dcat[ntrk]);
	  if(false) cout << ihit << " out of " << nhits << ", " << layer << " " << ladder << " " << sensor << " " << endl;
	  if(layer==0) th1f_dcal_B0_ladder[ladder]->Fill(d_dcal[ntrk]);
	  if(layer==1) th1f_dcal_B1_ladder[ladder]->Fill(d_dcal[ntrk]);
	  if(layer==2) th1f_dcal_B2_ladder[ladder]->Fill(d_dcal[ntrk]);
	  if(layer==3) th1f_dcal_B3_ladder[ladder]->Fill(d_dcal[ntrk]);

	}


      int idch = sngl_svxcnt->getDchIndex();
      d_idch[ntrk] = idch;

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

      // emc
      d_deadmap[ntrk] = node_cnt->get_deadmap(idch);
      d_warnmap[ntrk] = node_cnt->get_warnmap(idch);
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

      // rich
      d_disp[ntrk] = node_cnt->get_disp(idch);
      d_n0[ntrk] = node_cnt->get_n0(idch);
      d_npe0[ntrk] = node_cnt->get_npe0(idch);
      d_crkchi2[ntrk] = node_cnt->get_chi2(idch);

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

