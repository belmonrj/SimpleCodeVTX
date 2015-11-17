#include <iostream>
#include <cmath>
#include <string>

#include "TH1.h"
#include "TFile.h"

#include "SimpleTreeSVXCNT.h"

#include "getClass.h"
#include "Fun4AllReturnCodes.h"
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



int SimpleTreeSVXCNT::InitRun(PHCompositeNode *topNode)
{

  d_nevent = 0;

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
  // pfoa_cut = 2.5; //in degrees
  // deltaZ_cut = 1.0;
  // deltaPhi_cut = 0.1;


  d_OutputFile = new TFile(d_outfilename.c_str(),"RECREATE");

  // --- create and register histograms here
  th1f_cent = new TH1I("th1f_cent","", 100,0,100);
  th1f_bbcz = new TH1F("th1f_bbcz","",100,-50,50);
  th1f_vtxz = new TH1F("th1f_vtxz","",100,-50,50);

  th1f_bbcsq = new TH1F("th1f_bbcsq","",3000,0,3000);
  th1f_bbcnq = new TH1F("th1f_bbcnq","",3000,0,3000);

  // x-axis is variable, y-axis is momentum, z-axis is centrality
  // elemcmatch: x-axis = dphi, y-axis = dz, z-axis = momentum

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

  th1f_nsvx = new TH1F("th1f_nsvx","",20,0,20);
  th1f_ndch = new TH1F("th1f_ndch","",20,0,20);

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
  int trigbitraw = lvl1->get_lvl1_trigraw();


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
  PHCentralTrack *trk = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if ( trk == NULL )
    {
      cerr << PHWHERE << "ERROR: can't find PHCentralTrack node, discarding event" << endl;
      return DISCARDEVENT;
    }


  // --- get the SvxCentralTrackList node
  SvxCentralTrackList *svxcnttrklist = findNode::getClass<SvxCentralTrackList>(topNode,"SvxCentralTrackList");
  if(!svxcnttrklist)
    {
      cerr << PHWHERE << " ERROR: can't find SVXcentralTrack, discarding event " << endl;
      return DISCARDEVENT;
    }


  // --- loop over Central Tracks
  int numberdchtrks = 0;
  for (unsigned int i=0; i<trk->get_npart();i++)
    {

      float mom = trk->get_mom(i);
      float dchpx = trk->get_px(i);
      float dchpy = trk->get_py(i);
      float pt = sqrt(dchpx*dchpx+dchpy*dchpy);
      float emce = trk->get_emce(i);
      if (pt > 0.1 && pt < 4) numberdchtrks++;
      int quality = trk->get_quality(i);
      if(quality != 31 && quality != 63) continue;

      float phi = trk->get_phi(i);
      float alpha = trk->get_alpha(i);

      int n0 = trk->get_n0(i);
      int sn0 = trk->get_sn0(i);

      float disp = trk->get_disp(i);
      float sdisp = trk->get_sdisp(i);
      int dcarm = trk->get_dcarm(i);
      int sect = trk->get_sect(i);
      int iarmsect = dcarm * 4 + sect;
      int verbosity = 0;
      if(iarmsect<0&&verbosity>1) cout << "Houston, we have a problem" << endl;

      //float E = trk->get_ecore(i);
      //float dep = E/mom-1;
      float dep = trk->get_dep(i);

      //bool epcut = (dep > min_dep_cut && dep < max_dep_cut);

      float emcdphi = trk->get_emcsdphi_e(i);
      float emcdz = trk->get_emcsdz_e(i);
      //bool emcmatchcut = (fabs(emcdphi)<emcdphi_cut && fabs(emcdz)<emcdz_cut);
      float zvtx   =  (vtxout->get_Vertex()).getZ();


      if(mom < min_mom_cut || mom > max_mom_cut) continue;

      // float emc_t0 = trk->get_temc(i);
      // float crk_t0 = trk->get_tcrk(i);
      // float scrk_t0 = trk->get_stcrk(i);

    } // loop over Central Tracks

  // --- trigger selection
  if ((trigbitscaled&0x00000010) == 16)
    {
      th1f_ndch->Fill(numberdchtrks);
    }



  int nsvxtrks = 0;
  for (int i=0; i<svxcnttrklist->get_nCentralTracks();i++)
    {

      SvxCentralTrack *svxcntrltrk = svxcnttrklist->getCentralTrack(i);

      if (!svxcntrltrk) {cout << " no central track" << endl; continue;}

      int dchindex = svxcntrltrk->getDchIndex();
      float nhit = svxcntrltrk->getNhits();
      if ( nhit <= 2 ) continue; // require 3 hit tracks...?

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

      if ( pt > 0.1 && pt < 4 ) nsvxtrks++;

    } // loop over SVX Central Tracks

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

