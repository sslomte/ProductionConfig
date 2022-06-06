#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TEfficiency.h"
#include <TLorentzVector.h>


//sig only
TH1F* h_Ngenjets=new TH1F("s_N_genjets","N gen jets;N_{gen jets}",10,0,10);
TH1F* h_Ngoodgenjets=new TH1F("s_N_goodgenjets","N good gen jets;N_{gen jets}",10,0,10);
TH1F* h_Ntightgenjets=new TH1F("s_N_tightgenjets","N tight gen jets;N_{gen jets}",10,0,10);
TH1F* h_Nrecojets=new TH1F("s_N_recojets","N reco jets;N_{reco jets}",40,0,40);
TH1F* h_Ngoodrecojets=new TH1F("s_N_goodrecojets","N good reco jets;N_{good reco jets}",40,0,40);
TH1F* h_Nmatchedrecojets=new TH1F("s_N_matchrecojets","N matched reco jets to gen jets in min dR;N_{matched jets}",10,0,10);
TH1F* h_Nmatchedgoodrecojets=new TH1F("s_N_matchgoodrecojets","N matched good reco jets to gen jets in min dR;N_{matched good jets}",10,0,10);

TH1F* h_Ndau_gj=new TH1F("s_Ndau_gj","N daughters for gen jets;N_{dau}",30,0,30);
TH1F* h_Ndau_ggj=new TH1F("s_Ndau_ggj","N daughters for good gen jets;N_{dau}",30,0,30);
TH1F* h_Ndau_rj=new TH1F("s_Ndau_rj","N daughters for reco jets;N_{dau}",30,0,30);
TH1F* h_Ndau_grj=new TH1F("s_Ndau_grj","N daughters for good reco jets;N_{dau}",30,0,30);
TH1F* h_Ndau_mrj=new TH1F("s_Ndau_mrj","N daughters for matched reco jets;N_{dau}",30,0,30);
TH1F* h_Ndau_mgrj=new TH1F("s_Ndau_mgrj","N daughters for matched good reco jets;N_{dau}",30,0,30);

TH1F* h_gj_Pt=new TH1F("s_gj_Pt","gen jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_gj_Eta=new TH1F("s_gj_Eta","gen jet Eta;#eta",100,-5,5);
TH1F* h_gj_Phi=new TH1F("s_gj_Phi","gen jet Phi;#phi",100,-5,5);
TH1F* h_gj_E=new TH1F("s_gj_E","gen jet E;E [GeV]",250,0,500);

TH1F* h_ggj_Pt=new TH1F("s_ggj_Pt","good gen jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_ggj_Eta=new TH1F("s_ggj_Eta","good gen jet Eta;#eta",100,-5,5);
TH1F* h_ggj_Phi=new TH1F("s_ggj_Phi","good gen jet Phi;#phi",100,-5,5);
TH1F* h_ggj_E=new TH1F("s_ggj_E","good gen jet E;E [GeV]",250,0,500);

TH1F* h_rj_Pt=new TH1F("s_rj_Pt","reco jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_rj_Eta=new TH1F("s_rj_Eta","reco jet Eta;#eta",100,-5,5);
TH1F* h_rj_Phi=new TH1F("s_rj_Phi","reco jet Phi;#phi",100,-5,5);
TH1F* h_rj_E=new TH1F("s_rj_E","reco jet E;E [GeV]",250,0,500);

TH1F* h_grj_Pt=new TH1F("s_grj_Pt","good reco jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_grj_Eta=new TH1F("s_grj_Eta","good reco jet Eta;#eta",100,-5,5);
TH1F* h_grj_Phi=new TH1F("s_grj_Phi","good reco jet Phi;#phi",100,-5,5);
TH1F* h_grj_E=new TH1F("s_grj_E","good reco jet E;E [GeV]",250,0,500);

TH1F* h_mrj_Pt=new TH1F("s_mrj_Pt","matched reco jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_mrj_Eta=new TH1F("s_mrj_Eta","matched reco jet Eta;#eta",100,-5,5);
TH1F* h_mrj_Phi=new TH1F("s_mrj_Phi","matched reco jet Phi;#phi",100,-5,5);
TH1F* h_mrj_E=new TH1F("s_mrj_E","matched reco jet E;E [GeV]",250,0,500);

TH1F* h_M_j0j1_rj=new TH1F("s_M_j0j1_rj","reco leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_M_j0j1_grj=new TH1F("s_M_j0j1_grj","good reco leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_M_j0j1_mrj=new TH1F("s_M_j0j1_mrj","matched reco leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_M_j0j1_mgrj=new TH1F("s_M_j0j1_mgrj","matched good reco leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_M_j0j1_gj=new TH1F("s_M_j0j1_gj","gen leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_M_j0j1_ggj=new TH1F("s_M_j0j1_ggj","good gen leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_M_j0j1_tgj=new TH1F("s_M_j0j1_tgj","tight gen leading di-jets inv mass;M [GeV]",125,0,250);

TH1F* h_M_jj_rj=new TH1F("s_M_jj_rj","reco di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_M_jj_grj=new TH1F("s_M_jj_grj","good reco di-jets inv mass;M [GeV]",125,0,250);

TH1F* h_cmrj_Pt=new TH1F("s_cmrj_Pt","corr matched reco jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_cmrj_Eta=new TH1F("s_cmrj_Eta","corr matched reco jet Eta;#eta",100,-5,5);
TH1F* h_cmrj_Phi=new TH1F("s_cmrj_Phi","corr matched reco jet Phi;#phi",100,-5,5);
TH1F* h_cmrj_E=new TH1F("s_cmrj_E","corr matched reco jet E;E [GeV]",250,0,500);
TH1F* h_cmgrj_Pt=new TH1F("s_cmgrj_Pt","corr matched reco jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_cmgrj_Eta=new TH1F("s_cmgrj_Eta","corr matched reco jet Eta;#eta",100,-5,5);
TH1F* h_cmgrj_Phi=new TH1F("s_cmgrj_Phi","corr matched reco jet Phi;#phi",100,-5,5);
TH1F* h_cmgrj_E=new TH1F("s_cmgrj_E","corr matched reco jet E;E [GeV]",250,0,500);


TH1F* h_M_cmrj=new TH1F("s_M_j0j1_cmrj","corr match reco leading di-jets inv mass;M[GeV]",200,0,400);
TH1F* h_M_cmgrj=new TH1F("s_M_j0j1_cmgrj","corr match good reco leading di-jets inv mass;M[GeV]",200,0,400);

TH1F* h_maxPTratio_mrj=new TH1F("s_maxPT_mrj","max PT ratio",62,-0.1,3.0);
TH1F* h_maxEratio_mrj=new TH1F("s_maxE_mrj","max PT ratio",62,-0.1,3.0);
TH1F* h_maxPTratio_mgrj=new TH1F("s_maxPT_mgrj","max PT ratio",62,-0.1,3.0);
TH1F* h_maxEratio_mgrj=new TH1F("s_maxE_mgrj","max PT ratio",62,-0.1,3.0);

TH1F* h_maxPTratio_grj=new TH1F("s_maxPT_grj","max PT ratio of daughter and jet",62,-0.1,3.0);
TH1F* h_maxEratio_grj=new TH1F("s_maxE_grj","max E ratio of daughter and jet",62,-0.1,3.0);


//sig and BIB
TH1F* h_sb_Ngenjets=new TH1F("sb_N_genjets","N gen jets;N_{gen jets}",10,0,10);
TH1F* h_sb_Ngoodgenjets=new TH1F("sb_N_goodgenjets","N good gen jets;N_{gen jets}",10,0,10);
TH1F* h_sb_Ntightgenjets=new TH1F("sb_N_tightgenjets","N tight gen jets;N_{gen jets}",10,0,10);
TH1F* h_sb_Nrecojets=new TH1F("sb_N_recojets","N reco jets;N_{reco jets}",40,0,40);
TH1F* h_sb_Ngoodrecojets=new TH1F("sb_N_goodrecojets","N good reco jets;N_{good reco jets}",40,0,40);
TH1F* h_sb_Nmatchedrecojets=new TH1F("sb_N_matchrecojets","N matched reco jets to gen jets in min dR;N_{matched jets}",10,0,10);
TH1F* h_sb_Nmatchedgoodrecojets=new TH1F("sb_N_matchgoodrecojets","N matched good reco jets to gen jets in min dR;N_{matched good jets}",10,0,10);
TH1F* h_sb_Ngen_matchreco_jets=new TH1F("sb_N_gen_matchreco_jets","N {gen, matchreco} jets",10,0,10);

TH1F* h_sb_Ndau_gj=new TH1F("sb_Ndau_gj","N daughters for gen jets;N_{dau}",30,0,30);
TH1F* h_sb_Ndau_ggj=new TH1F("sb_Ndau_ggj","N daughters for good gen jets;N_{dau}",30,0,30);
TH1F* h_sb_Ndau_rj=new TH1F("sb_Ndau_rj","N daughters for reco jets;N_{dau}",30,0,30);
TH1F* h_sb_Ndau_grj=new TH1F("sb_Ndau_grj","N daughters for good reco jets;N_{dau}",30,0,30);
TH1F* h_sb_Ndau_mrj=new TH1F("sb_Ndau_mrj","N daughters for matched reco jets;N_{dau}",30,0,30);
TH1F* h_sb_Ndau_mgrj=new TH1F("sb_Ndau_mgrj","N daughters for matched good reco jets;N_{dau}",30,0,30);

TH1F* h_sb_Ncmrj=new TH1F("sb_N_corrmatchrecojets","N corr matched reco jets",10,0,10);
TH1F* h_sb_Ncmgrj=new TH1F("sb_N_corrmatchgoodrecojets","N corr matched good reco jets",10,0,10);

TH1F* h_sb_gj_Pt=new TH1F("sb_gj_Pt","gen jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_sb_gj_Eta=new TH1F("sb_gj_Eta","gen jet Eta;#eta",100,-5,5);
TH1F* h_sb_gj_Phi=new TH1F("sb_gj_Phi","gen jet Phi;#phi",100,-5,5);
TH1F* h_sb_gj_E=new TH1F("sb_gj_E","gen jet E;E [GeV]",250,0,500);

TH1F* h_sb_ggj_Pt=new TH1F("sb_ggj_Pt","good gen jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_sb_ggj_Eta=new TH1F("sb_ggj_Eta","good gen jet Eta;#eta",100,-5,5);
TH1F* h_sb_ggj_Phi=new TH1F("sb_ggj_Phi","good gen jet Phi;#phi",100,-5,5);
TH1F* h_sb_ggj_E=new TH1F("sb_ggj_E","good gen jet E;E [GeV]",250,0,500);

TH1F* h_sb_rj_Pt=new TH1F("sb_rj_Pt","reco jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_sb_rj_Eta=new TH1F("sb_rj_Eta","reco jet Eta;#eta",100,-5,5);
TH1F* h_sb_rj_Phi=new TH1F("sb_rj_Phi","reco jet Phi;#phi",100,-5,5);
TH1F* h_sb_rj_E=new TH1F("sb_rj_E","reco jet E;E [GeV]",250,0,500);

TH1F* h_sb_grj_Pt=new TH1F("sb_grj_Pt","good reco jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_sb_grj_Eta=new TH1F("sb_grj_Eta","good reco jet Eta;#eta",100,-5,5);
TH1F* h_sb_grj_Phi=new TH1F("sb_grj_Phi","good reco jet Phi;#phi",100,-5,5);
TH1F* h_sb_grj_E=new TH1F("sb_grj_E","good reco jet E;E [GeV]",250,0,500);

TH1F* h_sb_mrj_Pt=new TH1F("sb_mrj_Pt","matched reco jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_sb_mrj_Eta=new TH1F("sb_mrj_Eta","matched reco jet Eta;#eta",100,-5,5);
TH1F* h_sb_mrj_Phi=new TH1F("sb_mrj_Phi","matched reco jet Phi;#phi",100,-5,5);
TH1F* h_sb_mrj_E=new TH1F("sb_mrj_E","matched reco jet E;E [GeV]",250,0,500);

TH1F* h_sb_cmrj_Pt=new TH1F("sb_cmrj_Pt","corr matched reco jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_sb_cmrj_Eta=new TH1F("sb_cmrj_Eta","corr matched reco jet Eta;#eta",100,-5,5);
TH1F* h_sb_cmrj_Phi=new TH1F("sb_cmrj_Phi","corr matched reco jet Phi;#phi",100,-5,5);
TH1F* h_sb_cmrj_E=new TH1F("sb_cmrj_E","corr matched reco jet E;E [GeV]",250,0,500);
TH1F* h_sb_cmgrj_Pt=new TH1F("sb_cmgrj_Pt","corr matched reco jet PT;p_{T} [GeV]",250,0,500);
TH1F* h_sb_cmgrj_Eta=new TH1F("sb_cmgrj_Eta","corr matched reco jet Eta;#eta",100,-5,5);
TH1F* h_sb_cmgrj_Phi=new TH1F("sb_cmgrj_Phi","corr matched reco jet Phi;#phi",100,-5,5);
TH1F* h_sb_cmgrj_E=new TH1F("sb_cmgrj_E","corr matched reco jet E;E [GeV]",250,0,500);


TH1F* h_sb_M_j0j1_rj=new TH1F("sb_M_j0j1_rj","reco leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_sb_M_j0j1_grj=new TH1F("sb_M_j0j1_grj","good reco leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_sb_M_j0j1_mrj=new TH1F("sb_M_j0j1_mrj","matched reco leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_sb_M_j0j1_mgrj=new TH1F("sb_M_j0j1_mgrj","matched good reco leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_sb_M_j0j1_gj=new TH1F("sb_M_j0j1_gj","gen leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_sb_M_j0j1_ggj=new TH1F("sb_M_j0j1_ggj","good gen leading di-jets inv mass;M [GeV]",125,0,250);
TH1F* h_sb_M_j0j1_tgj=new TH1F("sb_M_j0j1_tgj","tight gen leading di-jets inv mass;M [GeV]",125,0,250);

TH1F* h_sb_M_cmrj=new TH1F("sb_M_j0j1_cmrj","corr match reco leading di-jets inv mass;M[GeV]",150,0,300);
TH1F* h_sb_M_cmgrj=new TH1F("sb_M_j0j1_cmgrj","corr match good reco leading di-jets inv mass;M[GeV]",150,0,300);

TH1F* h_sb_maxPTratio_mrj=new TH1F("sb_maxPT_mrj","max PT ratio",62,-0.1,3.0);
TH1F* h_sb_maxEratio_mrj=new TH1F("sb_maxE_mrj","max PT ratio",62,-0.1,3.0);
TH1F* h_sb_maxPTratio_mgrj=new TH1F("sb_maxPT_mgrj","max PT ratio",62,-0.1,3.0);
TH1F* h_sb_maxEratio_mgrj=new TH1F("sb_maxE_mgrj","max PT ratio",62,-0.1,3.0);


TH1F* h_sb_maxPTratio_grj=new TH1F("sb_maxPT_grj","max PT ratio of daughter and jet",62,-0.1,3.0);
TH1F* h_sb_maxEratio_grj=new TH1F("sb_maxE_grj","max E ratio of daughter and jet",62,-0.1,3.0);



//Jet
/*TH3F* h3_theta_phi_ene_j=new TH3F("J_Theta_Phi_E",";#theta;#phi;E [GeV]",45,-0.5,4,40,-4,4,100,0,200);
TH2F* h2_theta_phi_j=new TH2F("J_Theta_Phi",";#theta;#phi",45,-0.5,4,40,-4,4);
TH1F* h_t_jet=new TH1F("J_T","Jet T;T[ns]",1002,-0.2,100);
TH1F* h_z_jet=new TH1F("J_Z","Jet Z;Z[mm]",100,-100,100);
TH1F* h_R_jet=new TH1F("J_R","Jet R;R[mm]",50,0,100);
TH2F* h2_z_R_jet=new TH2F("J_Z_R","Jet Z vs R;Z[mm];R[mm]",100,-100,100,50,0,100);
TH3F* h3_z_R_E_jet=new TH3F("J_Z_R_E",";Z[mm];R[mm];E[GeV]",100,-100,100,50,0,100,100,0,200);
TH3F* h3_z_R_t_jet=new TH3F("J_Z_R_T",";Z[mm];R[mm];T[ns]",100,-100,100,50,0,100,102,-0.2,10);
TH1F* h_Ndau_jet=new TH1F("J_Ndau",";N_{dau}",20,0,20);

//jpar
TH3F* h3_theta_phi_ene_jpar=new TH3F("j_Theta_Phi_E",";#theta;#phi;E [GeV]",45,-0.5,4,40,-4,4,100,0,200);
TH2F* h2_theta_phi_jpar=new TH2F("j_Theta_Phi",";#theta;#phi",45,-0.5,4,40,-4,4);
TH2F* h2_z_R_jpar=new TH2F("j_Z_R",";Z[mm];R[mm]",100,-100,100,50,0,100);
TH3F* h3_z_R_E_jpar=new TH3F("j_Z_R_E",";Z[mm];R[mm];E[GeV]",100,-100,100,50,0,100,100,0,200);
TH1F* h_pt_ratio_max=new TH1F("j_maxPT_ratio",";max p_{T} ratio",44,-0.1,1); 
TH1F* h_ene_ratio_max=new TH1F("j_maxE_ratio",";max E ratio",44,-0.1,1); 
TH1F* h_min_z_jpar=new TH1F("j_minZ",";Z_{min} [mm]",100,-5,5);
TH1F* h_max_z_jpar=new TH1F("j_maxZ",";Z_{max} [mm]",100,-5,5);
TH1F* h_min_R_jpar=new TH1F("j_minR",";R_{min} [mm]",400,0,20);
TH1F* h_max_R_jpar=new TH1F("j_maxR",";R_{max} [mm]",400,0,20);
TH1F* h_delZ_jpar=new TH1F("j_delZ",";#Delta Z [mm]",150,0,30);
TH1F* h_delR_jpar=new TH1F("j_delR",";#Delta R [mm]",150,0,30);
TH2F* h2_delZ_delR_jpar=new TH2F("j_delZ_delR",";#Delta Z [mm];#Delta R [mm]",100,0,5,100,0,5);
TH1F* h_delT_z_jpar=new TH1F("j_delT_z",";#Delta T [ns]",102,-0.2,10);
TH1F* h_delT_R_jpar=new TH1F("j_delT_R",";#Delta T [ns]",102,-0.2,10);
TH1F* h_delR_delZg05=new TH1F("j_delR_delZg05","#Delta R when #Delta Z > 0.5;#Delta R [mm]",100,0,20);
TH1F* h_Ndau_delZg05=new TH1F("J_Ndau_delZg05","N daughters when #Delta Z > 0.5;N_{dau}",20,0,20);
TH1F* h_delR_delZl05=new TH1F("j_delR_delZl05","#Delta R when #Delta Z <= 0.5;#Delta R [mm]",100,0,20);
TH1F* h_Ndau_delZl05=new TH1F("J_Ndau_delZl05","N daughters when #Delta Z <= 0.5;N_{dau}",20,0,20);


TH1F* h_Nj=new TH1F("J_Nj","N jets;N_{jets}",40,0,40);

//after dZ and dR
TH1F* h_Ngj=new TH1F("J_Ngj","N jets with 2<dZ<3 and dR>2.5;N_{jets}",40,0,40);
TH1F* h_jj_mass=new TH1F("J_jj_mass","di-jet mass;M [GeV]",100,0,200);
TH1F* h_jj_p=new TH1F("J_jj_p","di-jet momentum;P [GeV]",100,0,200);

TH1F* h_Ngjs=new TH1F("J_Ngjs","N jets with 2<dZ<3 and dR>2.5;N_{jets}",40,0,40);
TH1F* h_dijet_M=new TH1F("J_dijet_mass","di-jet mass for good jets;M [GeV]",100,0,200);
TH1F* h_dijet_Pt=new TH1F("J_dijet_pt","di-jet pt for good jets;M [GeV]",100,0,200);
*/
 

//find min from a list
int FindMin(double input_arr[], int size)
  {
    float min = 100000000.0; int ind=0;
    for(int k=0;k<size;k++){
      if(input_arr[k]<min){
	min = input_arr[k];
	ind = k;
      }
    }   
    return ind;
  }
//find max from a list
int FindMax(double input_arr[], int size)
  {
    float max = 0.0000000001; int ind=0;
    for(int k=0;k<size;k++){
      if(input_arr[k]>max){
	max = input_arr[k];
	ind = k;
      }
    }   
    return ind;
  }
//find max from vector
int findmax(vector <double> vec)
{
  double max = 0.0; int ind=0;
  for(int i=0; i<vec.size(); i++){
    if(vec.at(i)>max){
      max = vec.at(i);
      ind=i;
    }
  }
  return ind;
}


//function to find dR given P4 of two candidates                                                                                                                                                            
float find_delR(TLorentzVector P4_1, TLorentzVector P4_2)
{
  float dEta = abs(P4_1.Eta() - P4_2.Eta());
  float dPhi = abs(P4_1.Phi() - P4_2.Phi());
  if (dPhi>TMath::Pi()){dPhi = 2*TMath::Pi() - dPhi;}
  float dR = sqrt( pow(dEta, 2) + pow(dPhi, 2) );
  return dR;
}


//function to find best match mcpar for given reco candidate                                                                                                                                                
int FindRecoMatch(vector <TLorentzVector> recos, TLorentzVector genCandP4)
{
  float dR_min = 999.0;
  int ind = -1;
  for(int i=0; i<recos.size(); i++){
    float dR_tmp = find_delR(recos.at(i), genCandP4);
    if(dR_tmp<dR_min){
      dR_min=dR_tmp;
      ind = i;
    }
  }
  return ind;
}

//find average
double findAvg(vector <double> arr)
{
  double sum=0.0; int n=arr.size();
  for(int a=0; a<arr.size(); a++){
    sum = sum + arr.at(a);
  }
  double avg = sum/n;
  return avg;
}

//find std dev
double findSD(vector <double> arr)
{
  int n=arr.size(); double avgg = findAvg(arr); double sumsq=0.0;
  for(int b=0; b<arr.size(); b++){
    sumsq = sumsq + pow((arr.at(b) - avgg),2); 
  }
  double stddev = sqrt(sumsq/(n-1));
  return stddev;
}


void analysis()
{
  

  double mh=125.;
 // Declaration of leaf types

  Int_t           ndaughters_s_gen[10];
  Float_t         daughters_PX_s_gen[10][60];   
  Float_t         daughters_PY_s_gen[10][60];   
  Float_t         daughters_PZ_s_gen[10][60];   
  Float_t         daughters_M_s_gen[10][60];   
  Float_t         daughters_E_s_gen[10][60];   

  Int_t           ndaughters_s_reco[10];
  Float_t         daughters_PX_s_reco[10][60];   
  Float_t         daughters_PY_s_reco[10][60];   
  Float_t         daughters_PZ_s_reco[10][60];   
  Float_t         daughters_M_s_reco[10][60];   
  Float_t         daughters_E_s_reco[10][60];   

  Int_t           nj_s_gen;
  Float_t         jmox_s_gen[10];   //[njet]                                                                                                                                                              
  Float_t         jmoy_s_gen[10];   //[njet]                                                                                                                                                                
  Float_t         jmoz_s_gen[10];   //[njet]                                                                                                                                                                
  Float_t         jmas_s_gen[10];   //[njet]                                                                                                                                                                
  Float_t         jene_s_gen[10];   //[njet]                                                                                                                                                                
  Float_t         jcost_s_gen[10];   //[njet]                                                                                                                                                               
  Float_t         jmom_s_gen[10];   //[njet]                                                                                                                                                                

  Int_t           nj_s_reco;
  Float_t         jmox_s_reco[10];   //[njet]                                                                                                                                                               
  Float_t         jmoy_s_reco[10];   //[njet]                                                                                                                                                               
  Float_t         jmoz_s_reco[10];   //[njet]                                                                                                                                                               
  Float_t         jmas_s_reco[10];   //[njet]                                                                                                                                                               
  Float_t         jene_s_reco[10];   //[njet]                                                                                                                                                               
  Float_t         jcost_s_reco[10];   //[njet]                                                                                                                                                            
  Float_t         jmom_s_reco[10];   //[njet]                                                                                                                                                               
  


  // List of branches
  TBranch        *b_njet;   
  TBranch        *b_jmox;   
  TBranch        *b_jmoy;   
  TBranch        *b_jmoz;   
  TBranch        *b_jmas;   
  TBranch        *b_jene;   
  TBranch        *b_jcost;   
  TBranch        *b_jmom;   

  TBranch        *b_ndaughters;   
  TBranch        *b_daughters_PX;
  TBranch        *b_daughters_PY;
  TBranch        *b_daughters_PZ;
  TBranch        *b_daughters_M;
  TBranch        *b_daughters_E;


  TChain* fChain_s = new TChain("fChain_s");
  fChain_s->Add("JetHistogramsHbb_sig_only_bib100_2mev_250ps_ak5_10kEvents.root/JetHistogramRecoJetTuple");   

  TChain* fChain_s_gen = new TChain("fChain_s_gen");
  TChain* fChain_s_reco = new TChain("fChain_s_reco");
  fChain_s_gen->Add("JetHistogramsHbb_sig_only_bib100_2mev_250ps_ak5_10kEvents.root/JetHistogramGenJetTuple");
  fChain_s_reco->Add("JetHistogramsHbb_sig_only_bib100_2mev_250ps_ak5_10kEvents.root/JetHistogramRecoJetTuple");
  

  fChain_s_gen->SetBranchAddress("ndaughters", &ndaughters_s_gen, &b_ndaughters);
  fChain_s_gen->SetBranchAddress("daughters_PX", daughters_PX_s_gen, &b_daughters_PX);
  fChain_s_gen->SetBranchAddress("daughters_PY", daughters_PY_s_gen, &b_daughters_PY);
  fChain_s_gen->SetBranchAddress("daughters_PZ", daughters_PZ_s_gen, &b_daughters_PZ);
  fChain_s_gen->SetBranchAddress("daughters_M", daughters_M_s_gen, &b_daughters_M);
  fChain_s_gen->SetBranchAddress("daughters_E", daughters_E_s_gen, &b_daughters_E);

  fChain_s_reco->SetBranchAddress("ndaughters", &ndaughters_s_reco, &b_ndaughters);
  fChain_s_reco->SetBranchAddress("daughters_PX", daughters_PX_s_reco, &b_daughters_PX);
  fChain_s_reco->SetBranchAddress("daughters_PY", daughters_PY_s_reco, &b_daughters_PY);
  fChain_s_reco->SetBranchAddress("daughters_PZ", daughters_PZ_s_reco, &b_daughters_PZ);
  fChain_s_reco->SetBranchAddress("daughters_M", daughters_M_s_reco, &b_daughters_M);
  fChain_s_reco->SetBranchAddress("daughters_E", daughters_E_s_reco, &b_daughters_E);

  
  fChain_s_gen->SetBranchAddress("nj", &nj_s_gen, &b_njet);
  fChain_s_gen->SetBranchAddress("jmox", jmox_s_gen, &b_jmox);
  fChain_s_gen->SetBranchAddress("jmoy", jmoy_s_gen, &b_jmoy);
  fChain_s_gen->SetBranchAddress("jmoz", jmoz_s_gen, &b_jmoz);
  fChain_s_gen->SetBranchAddress("jmas", jmas_s_gen, &b_jmas);
  fChain_s_gen->SetBranchAddress("jene", jene_s_gen, &b_jene);
  fChain_s_gen->SetBranchAddress("jcost", jcost_s_gen, &b_jcost);
  fChain_s_gen->SetBranchAddress("jmom", jmom_s_gen, &b_jmom);
  

  fChain_s_reco->SetBranchAddress("nj", &nj_s_reco, &b_njet);
  fChain_s_reco->SetBranchAddress("jmox", jmox_s_reco, &b_jmox);
  fChain_s_reco->SetBranchAddress("jmoy", jmoy_s_reco, &b_jmoy);
  fChain_s_reco->SetBranchAddress("jmoz", jmoz_s_reco, &b_jmoz);
  fChain_s_reco->SetBranchAddress("jmas", jmas_s_reco, &b_jmas);
  fChain_s_reco->SetBranchAddress("jene", jene_s_reco, &b_jene);
  fChain_s_reco->SetBranchAddress("jcost", jcost_s_reco, &b_jcost);
  fChain_s_reco->SetBranchAddress("jmom", jmom_s_reco, &b_jmom);
  

  //----------------------------------------------------------------------------------------------------------------------------

  vector <TLorentzVector> jets;
  vector <TLorentzVector> goodjets;
  vector <TLorentzVector> genjets;
  vector <TLorentzVector> goodgenjets;
  vector <TLorentzVector> tightgenjets;
  vector <TLorentzVector> recojets;
  vector <TLorentzVector> goodrecojets;
  vector <TLorentzVector> matchedrecojets;
  vector <TLorentzVector> matchedgoodrecojets;
  vector <TLorentzVector> goodcorrmatchrecojets;
  vector <vector<TLorentzVector>> gen_matchreco_jets;

  vector <double> pt_ratio_s_reco;
  vector <double> e_ratio_s_reco;
  vector <double> pt_ratio_mrj;
  vector <double> e_ratio_mrj;


  vector <double> Tf_eta0To0p5_pt10To50; vector <double> Tf_eta0p5To1p0_pt10To50; vector <double> Tf_eta1p0To2p5_pt10To50; 
  vector <double> Tf_eta0To0p5_pt50To100; vector <double> Tf_eta0p5To1p0_pt50To100; vector <double> Tf_eta1p0To2p5_pt50To100; 
  vector <double> Tf_eta0To0p5_pt100; vector <double> Tf_eta0p5To1p0_pt100; vector <double> Tf_eta1p0To2p5_pt100; 

  Tf_eta0To0p5_pt10To50.clear(); Tf_eta0p5To1p0_pt10To50.clear(); Tf_eta1p0To2p5_pt10To50.clear();
  Tf_eta0To0p5_pt50To100.clear(); Tf_eta0p5To1p0_pt50To100.clear(); Tf_eta1p0To2p5_pt50To100.clear();
  Tf_eta0To0p5_pt100.clear(); Tf_eta0p5To1p0_pt100.clear(); Tf_eta1p0To2p5_pt100.clear();

  //s_jetpar

  /*  for(unsigned int ientry=0; ientry<fChain_s->GetEntries(); ++ientry){
    fChain_s->GetEntry(ientry);
   
    jets.clear();
    goodjets.clear();


    //cout<<"-----------Event no:  "<<ientry<<endl;
    //for each jet, access its daughter particles.
    //Access their Px, Py, Pz, E. Construct theta and phi
    //then find dTheta and dPhi from jet axis 
    //construct cone variables: width and length of cone, E ratio of daughters with parent jet 
    //for jetpar no x,y,z info since pf object, not hit. But can construct x,y,z using TLorentzvector

    
    h_Nj->Fill(nj_s);
    for(int i=0; i<nj_s; i++){

      double jPhi=TMath::ATan(jmoy_s[i]/jmox_s[i]);
      double jTheta=TMath::ACos(jcost_s[i]);
      double jPt=sqrt((jmox_s[i]*jmox_s[i]) + (jmoy_s[i]*jmoy_s[i]));
      //double jEta=-TMath::LogE(TMath::Tan(jTheta/2));

      h2_theta_phi_j->Fill(jTheta, jPhi);
      h3_theta_phi_ene_j->Fill(jTheta, jPhi ,jene_s[i]);

      //Jet TLorentzVector
      TLorentzVector Jet; Jet.SetPxPyPzE(jmox_s[i], jmoy_s[i], jmoz_s[i], jene_s[i]);
      double Jet_x=Jet.X(); double Jet_y=Jet.Y(); double Jet_z=Jet.Z(); double Jet_t=Jet.T();
      double Jet_R=sqrt(Jet_x*Jet_x + Jet_y*Jet_y);
      h_t_jet->Fill(Jet_t);
      h3_z_R_t_jet->Fill(Jet_z, Jet_R, Jet_t);
      h_z_jet->Fill(Jet_z);
      h_R_jet->Fill(Jet_R);
      h2_z_R_jet->Fill(Jet_z, Jet_R);
      h3_z_R_E_jet->Fill(Jet_z, Jet_R, jene_s[i]);

      jets.push_back(Jet);


      h_Ndau_jet->Fill(ndaughters[i]);


      double pt_ratio_arr[ndaughters[i]]; double ene_ratio_arr[ndaughters[i]]; double z_jpar_arr[ndaughters[i]]; double R_jpar_arr[ndaughters[i]]; double t_jpar_arr[ndaughters[i]];
      for(int j=0; j<ndaughters[i]; j++){     
	double jpar_mom=sqrt(daughters_PX[i][j]*daughters_PX[i][j] + daughters_PY[i][j]*daughters_PY[i][j] + daughters_PZ[i][j]*daughters_PZ[i][j] );
	double jpar_pt=sqrt(daughters_PX[i][j]*daughters_PX[i][j] + daughters_PY[i][j]*daughters_PY[i][j] );
	double jpar_theta=TMath::ACos(daughters_PZ[i][j]/jpar_mom);
	double jpar_phi=TMath::ATan(daughters_PY[i][j]/daughters_PX[i][j]);

	h2_theta_phi_jpar->Fill(jpar_theta, jpar_phi);
	h3_theta_phi_ene_jpar->Fill(jpar_theta, jpar_phi, daughters_E[i][j]);

	//cout<<"dau_Pt=  "<<jpar_pt<<"    dau_P= "<<jpar_mom<<"   and dau_Pz= "<<daughters_PZ[i][j]<<endl;

	//TLorentzVector
	TLorentzVector jpar; jpar.SetPxPyPzE(daughters_PX[i][j], daughters_PY[i][j], daughters_PZ[i][j], daughters_E[i][j]);
	double jpar_x=jpar.X(); double jpar_y=jpar.Y(); double jpar_z=jpar.Z(); double jpar_t=jpar.T();
	double jpar_R=sqrt(jpar_x*jpar_x + jpar_y*jpar_y);
	z_jpar_arr[j]=jpar_z; R_jpar_arr[j]=jpar_R; t_jpar_arr[j]=jpar_t;

	h2_z_R_jpar->Fill(jpar_z, jpar_R);
	h3_z_R_E_jpar->Fill(jpar_z, jpar_R, daughters_E[i][j]);
       

	//dTheta and dPhi for each daughter from parent jet
	double dTheta=abs(jpar_theta - jTheta);
	double dPhi=abs(jpar_phi - jPhi);
	double dPt=abs(jpar_pt - jPt);
	double dP=abs(jpar_mom - jmom_s[i]);
	double dEne=abs(daughters_E[i][j] - jene_s[i]);

	//PT and E ratio 
	double PT_ratio=jpar_pt/jPt;
	double E_ratio=daughters_E[i][j]/jene_s[i];
	pt_ratio_arr[j]=PT_ratio;
	ene_ratio_arr[j]=E_ratio;

      }//daughter loop

      //find max of pt_ratio_arr and ene_ratio_array
      double max_pt_ratio=pt_ratio_arr[FindMax(pt_ratio_arr,sizeof(pt_ratio_arr)/sizeof(pt_ratio_arr[0]))];
      double max_ene_ratio=ene_ratio_arr[FindMax(ene_ratio_arr,sizeof(ene_ratio_arr)/sizeof(ene_ratio_arr[0]))];

      if(jPt>5){
	h_pt_ratio_max->Fill(max_pt_ratio);
	h_ene_ratio_max->Fill(max_ene_ratio);
      }
      //find min of dTheta and dPhi


      //find min and max of z_jpar
      double min_z_jpar=z_jpar_arr[FindMin(z_jpar_arr,sizeof(z_jpar_arr)/sizeof(z_jpar_arr[0]))]; double min_z_t_jpar=t_jpar_arr[FindMin(z_jpar_arr,sizeof(z_jpar_arr)/sizeof(z_jpar_arr[0]))];
      double max_z_jpar=z_jpar_arr[FindMax(z_jpar_arr,sizeof(z_jpar_arr)/sizeof(z_jpar_arr[0]))]; double max_z_t_jpar=t_jpar_arr[FindMax(z_jpar_arr,sizeof(z_jpar_arr)/sizeof(z_jpar_arr[0]))];
      //find min and max of R_jpar
      double min_R_jpar=R_jpar_arr[FindMin(R_jpar_arr,sizeof(R_jpar_arr)/sizeof(R_jpar_arr[0]))]; double min_R_t_jpar=t_jpar_arr[FindMin(R_jpar_arr,sizeof(R_jpar_arr)/sizeof(R_jpar_arr[0]))];
      double max_R_jpar=R_jpar_arr[FindMax(R_jpar_arr,sizeof(R_jpar_arr)/sizeof(R_jpar_arr[0]))]; double max_R_t_jpar=t_jpar_arr[FindMax(R_jpar_arr,sizeof(R_jpar_arr)/sizeof(R_jpar_arr[0]))];
      h_min_z_jpar->Fill(min_z_jpar);
      h_max_z_jpar->Fill(max_z_jpar);
      h_min_R_jpar->Fill(min_R_jpar);
      h_max_R_jpar->Fill(max_R_jpar);
      double delZ=abs(max_z_jpar - min_z_jpar);
      double delR=abs(max_R_jpar - min_R_jpar);
      h_delZ_jpar->Fill(delZ);
      h_delR_jpar->Fill(delR);
      double delT_z=abs(max_z_t_jpar - min_z_t_jpar);
      double delT_R=abs(max_R_t_jpar - min_R_t_jpar);
      h_delT_z_jpar->Fill(delT_z);
      h_delT_R_jpar->Fill(delT_R);
      h2_delZ_delR_jpar->Fill(delZ, delR);


      if(delZ>0.5){
	h_delR_delZg05->Fill(delR);
	h_Ndau_delZg05->Fill(ndaughters[i]);
      }
      else {
	h_delR_delZl05->Fill(delR);
	h_Ndau_delZl05->Fill(ndaughters[i]);
      }

      //cout<<"dZ=="<<delZ<<"    dR=="<<delR<<endl;
      //apply delZ and delR cuts and then plot dijet mass
      if(delZ>2 && delZ<3){
	if(delR>2.5){
	  Ngj++;
	  goodjets.push_back(Jet);
	}
      }





    }//jet loop


    //dijet mass
    if(jets.size()>1){
      double jj_mass = (jets.at(0) + jets.at(1)).M();
      double jj_p = (jets.at(0) + jets.at(1)).P();
      h_jj_mass->Fill(jj_mass);
      h_jj_p->Fill(jj_p);
      
    }


    //number of jets satisfying dZ and dR cuts:
    h_Ngj->Fill(Ngj);
    h_Ngjs->Fill(goodjets.size());
    if(goodjets.size()>1){
      for(int h=0;h<goodjets.size();h++){
	for(int k=h+1;k<goodjets.size();k++){
	  double dijet_mass = (goodjets.at(h) + goodjets.at(k)).M();
	  double dijet_pt = (goodjets.at(h) + goodjets.at(k)).Pt();
	  h_dijet_M->Fill(dijet_mass);
	  h_dijet_Pt->Fill(dijet_pt);
	}
      }
    }

  }//event loop
  */

  
  for(unsigned int ientry=0; ientry<fChain_s_gen->GetEntries(); ++ientry){
    fChain_s_gen->GetEntry(ientry);
    fChain_s_reco->GetEntry(ientry);

    genjets.clear();
    goodgenjets.clear();
    tightgenjets.clear();
    recojets.clear();
    goodrecojets.clear();
    matchedrecojets.clear();
    matchedgoodrecojets.clear();
    gen_matchreco_jets.clear();

    pt_ratio_s_reco.clear();
    e_ratio_s_reco.clear();
    pt_ratio_mrj.clear();
    e_ratio_mrj.clear();

    for(int a=0; a<nj_s_gen; a++){
      TLorentzVector gen; gen.SetPxPyPzE(jmox_s_gen[a], jmoy_s_gen[a], jmoz_s_gen[a], jene_s_gen[a]); gen.SetTheta(TMath::ACos(jcost_s_gen[a]));
      genjets.push_back(gen);
      h_gj_Pt->Fill(gen.Pt());
      h_gj_Eta->Fill(gen.Eta());
      h_gj_Phi->Fill(gen.Phi());
      h_gj_E->Fill(gen.E());

      //acceptance cuts
      if(gen.Pt()>10 && abs(gen.Theta())>0.1745 && abs(gen.Theta())<2.9671){
	goodgenjets.push_back(gen);
	h_ggj_Pt->Fill(gen.Pt());
	h_ggj_Eta->Fill(gen.Eta());
	h_ggj_Phi->Fill(gen.Phi());
	h_ggj_E->Fill(gen.E());	

	h_Ndau_ggj->Fill(ndaughters_s_gen[a]);

	if(gen.Pt()>20 && abs(gen.Theta())>0.1745 && abs(gen.Theta())<2.9671){
	  tightgenjets.push_back(gen);
	}
      }
      //daughters
      h_Ndau_gj->Fill(ndaughters_s_gen[a]);      
    }
    h_Ngenjets->Fill(genjets.size());
    h_Ngoodgenjets->Fill(goodgenjets.size());
    h_Ntightgenjets->Fill(tightgenjets.size());
    
    //reco jets
    for(int a=0; a<nj_s_reco; a++){
      TLorentzVector reco; reco.SetPxPyPzE(jmox_s_reco[a], jmoy_s_reco[a], jmoz_s_reco[a], jene_s_reco[a]); reco.SetTheta(TMath::ACos(jcost_s_reco[a]));
      recojets.push_back(reco);
      h_rj_Pt->Fill(reco.Pt());
      h_rj_Eta->Fill(reco.Eta());
      h_rj_Phi->Fill(reco.Phi());
      h_rj_E->Fill(reco.E());
                 
      //daughters
      h_Ndau_rj->Fill(ndaughters_s_reco[a]);
      for(int b=0; b<ndaughters_s_reco[a]; b++){
	TLorentzVector dau; dau.SetPxPyPzE(daughters_PX_s_reco[a][b],daughters_PY_s_reco[a][b],daughters_PZ_s_reco[a][b],daughters_E_s_reco[a][b]);
	double pt_ratio = dau.Pt()/reco.Pt();
	pt_ratio_s_reco.push_back(pt_ratio);
	double e_ratio = dau.E()/reco.E();
	e_ratio_s_reco.push_back(e_ratio);	
      }

      //acceptance cuts
      if(reco.Pt()>10 && abs(reco.Theta())>0.1745 && abs(reco.Theta())<2.9671){
	goodrecojets.push_back(reco);
	h_grj_Pt->Fill(reco.Pt());
	h_grj_Eta->Fill(reco.Eta());
	h_grj_Phi->Fill(reco.Phi());
	h_grj_E->Fill(reco.E());

	h_Ndau_grj->Fill(ndaughters_s_reco[a]);

      //find max pt_ratio and e_ratio
	if(pt_ratio_s_reco.size()>0){
	  double max_pt_ratio = pt_ratio_s_reco.at(findmax(pt_ratio_s_reco));
	  double max_e_ratio = e_ratio_s_reco.at(findmax(e_ratio_s_reco));
	  h_maxPTratio_grj->Fill(max_pt_ratio);
	  h_maxEratio_grj->Fill(max_e_ratio);
	}
      }

    }
    h_Nrecojets->Fill(recojets.size());
    h_Ngoodrecojets->Fill(goodrecojets.size());

  
    int Ndau_mrj; double max_pt_ratio_mrj; double max_e_ratio_mrj;
    //matching recojets to genjets in min dR    
    if(genjets.size()>0 && recojets.size()>0){
      for(int j=0; j<genjets.size(); j++){
	int recoindex = FindRecoMatch(recojets, genjets.at(j));
	double dR = find_delR(genjets.at(j), recojets.at(recoindex));
	if(dR<0.5){
	  matchedrecojets.push_back(recojets.at(recoindex));	
	  
	  //daughters
	  Ndau_mrj = ndaughters_s_reco[recoindex];
	  h_Ndau_mrj->Fill(Ndau_mrj);
	  
	  //pt_ratio 
	  for(int k=0; k<ndaughters_s_reco[recoindex]; k++){
	    TLorentzVector dau; dau.SetPxPyPzE(daughters_PX_s_reco[recoindex][k],daughters_PY_s_reco[recoindex][k],daughters_PZ_s_reco[recoindex][k],daughters_E_s_reco[recoindex][k]);
	    double dau_pt_ratio_mrj = dau.Pt()/recojets.at(recoindex).Pt();
	    double dau_e_ratio_mrj = dau.E()/recojets.at(recoindex).E();
	    pt_ratio_mrj.push_back(dau_pt_ratio_mrj);
	    e_ratio_mrj.push_back(dau_e_ratio_mrj);
	  }
	  //findmax
	  max_pt_ratio_mrj = pt_ratio_mrj.at(findmax(pt_ratio_mrj));
	  max_e_ratio_mrj = e_ratio_mrj.at(findmax(e_ratio_mrj));
	  h_maxPTratio_mrj->Fill(max_pt_ratio_mrj);
	  h_maxEratio_mrj->Fill(max_e_ratio_mrj);

	  vector<TLorentzVector> vec; vec.push_back(genjets.at(j)); vec.push_back(recojets.at(recoindex));
	  gen_matchreco_jets.push_back(vec);
	}
      }
    }
  
    h_Nmatchedrecojets->Fill(matchedrecojets.size());
    for(int g=0; g<matchedrecojets.size(); g++){
      h_mrj_Pt->Fill(matchedrecojets.at(g).Pt());
      h_mrj_Eta->Fill(matchedrecojets.at(g).Eta());
      h_mrj_Phi->Fill(matchedrecojets.at(g).Phi());
      h_mrj_E->Fill(matchedrecojets.at(g).E());
      //acceptance cuts
      if(matchedrecojets.at(g).Pt()>10 && abs(matchedrecojets.at(g).Theta())>0.1745 && abs(matchedrecojets.at(g).Theta())<2.9671){
	matchedgoodrecojets.push_back(matchedrecojets.at(g));
	
	//daughters
	h_Ndau_mgrj->Fill(Ndau_mrj);
	h_maxPTratio_mgrj->Fill(max_pt_ratio_mrj);
	h_maxEratio_mgrj->Fill(max_e_ratio_mrj);	
	
      }
    }
    h_Nmatchedgoodrecojets->Fill(matchedgoodrecojets.size()); 


    //inv_mass of di-jet pair 1) two leading jets
    //reco jets
    if(recojets.size()>1){
      double inv_m1 = (recojets.at(0) + recojets.at(1)).M();
      h_M_j0j1_rj->Fill(inv_m1);
    }
    if(goodrecojets.size()>1){
      double inv_m2 = (goodrecojets.at(0) + goodrecojets.at(1)).M();
      h_M_j0j1_grj->Fill(inv_m2);
    }
    if(matchedrecojets.size()>1){
      double inv_m3 = (matchedrecojets.at(0) + matchedrecojets.at(1)).M();
      h_M_j0j1_mrj->Fill(inv_m3);
    }
    if(matchedgoodrecojets.size()>1){
      double inv_m4 = (matchedgoodrecojets.at(0) + matchedgoodrecojets.at(1)).M();
      h_M_j0j1_mgrj->Fill(inv_m4);
    }


    //gen jets
    if(genjets.size()>1){
      double inv_m_g = (genjets.at(0) + genjets.at(1)).M();
      h_M_j0j1_gj->Fill(inv_m_g);
    }
    if(goodgenjets.size()>1){
      double inv_m_gg = (goodgenjets.at(0) + goodgenjets.at(1)).M();
      h_M_j0j1_ggj->Fill(inv_m_gg);
    }
    if(tightgenjets.size()>1){
      double inv_m_tg = (tightgenjets.at(0) + tightgenjets.at(1)).M();
      h_M_j0j1_tgj->Fill(inv_m_tg);
    }
    

    //inv_mass of di-jet pair 2) closest to 125GeV
    if(recojets.size()>1){
      for(int c=0; c<recojets.size(); c++){
	for(int d=c+1; d<recojets.size(); d++){
	  double invm = (recojets.at(c) + recojets.at(d)).M(); 
	  double dM = abs(invm - 125.);
	  h_M_jj_rj->Fill(invm);

	  //delM.[];
	}
      }
    }

    if(goodrecojets.size()>1){
      for(int c=0; c<goodrecojets.size(); c++){
	for(int d=c+1; d<goodrecojets.size(); d++){
	  double invm2 = (goodrecojets.at(c) + goodrecojets.at(d)).M(); 
	  double dM2 = abs(invm2 - 125.);
	  h_M_jj_grj->Fill(invm2);

	  //delM.[];
	}
      }
    }



    //-----------------------------------    
    //Jet energy correction
    //find scale factor = PT(gen)/PT(reco) in eta-PT regions
    //let's define eta-PT regions
    //also loop over gen jets; with first have matchedrecojets to genjets in dR<0.5
  


    //transfer factors calculation here
    if(gen_matchreco_jets.size()>0){
      for(int k=0; k<gen_matchreco_jets.size(); k++){
	double genjEta = abs(gen_matchreco_jets[k][0].Eta()); double genjPt = gen_matchreco_jets[k][0].Pt();

	//eta [0, 0.5] & PT [10,50,100,>100]
	if(genjEta>=0 && genjEta<0.5){
	  if(genjPt>=10 && genjPt<50){
	    Tf_eta0To0p5_pt10To50.push_back(gen_matchreco_jets[k][0].Pt()/gen_matchreco_jets[k][1].Pt());
	  }
	  if(genjPt>=50 && genjPt<100){
	    Tf_eta0To0p5_pt50To100.push_back(gen_matchreco_jets[k][0].Pt()/gen_matchreco_jets[k][1].Pt());
	  }
	  if(genjPt>=100){
	    Tf_eta0To0p5_pt100.push_back(gen_matchreco_jets[k][0].Pt()/gen_matchreco_jets[k][1].Pt());
	  }	
	}

	//eta [0.5, 1.0] & PT [10,50,100,>100]
	if(genjEta>=0.5 && genjEta<1.0){
	  if(genjPt>=10 && genjPt<50){
	    Tf_eta0p5To1p0_pt10To50.push_back(gen_matchreco_jets[k][0].Pt()/gen_matchreco_jets[k][1].Pt());
	  }
	  if(genjPt>=50 && genjPt<100){
	    Tf_eta0p5To1p0_pt50To100.push_back(gen_matchreco_jets[k][0].Pt()/gen_matchreco_jets[k][1].Pt());
	  }
	  if(genjPt>=100){
	    Tf_eta0p5To1p0_pt100.push_back(gen_matchreco_jets[k][0].Pt()/gen_matchreco_jets[k][1].Pt());
	  }	
	}

	//eta [1.0, 2.5] & PT [10,50,100,>100]
	if(genjEta>=1.0 && genjEta<2.5){
	  if(genjPt>=10 && genjPt<50){
	    Tf_eta1p0To2p5_pt10To50.push_back(gen_matchreco_jets[k][0].Pt()/gen_matchreco_jets[k][1].Pt());
	  }
	  if(genjPt>=50 && genjPt<100){
	    Tf_eta1p0To2p5_pt50To100.push_back(gen_matchreco_jets[k][0].Pt()/gen_matchreco_jets[k][1].Pt());

	  }
	  if(genjPt>=100){
	    Tf_eta1p0To2p5_pt100.push_back(gen_matchreco_jets[k][0].Pt()/gen_matchreco_jets[k][1].Pt());
	  }
	}
      }
    }
    //transfer factors filled


  }//evt loop
  //-----------------------------------    
    
  /*cout<<"size of Tf_eta0To0p5_pt10To50=  "<<Tf_eta0To0p5_pt10To50.size()<<endl;
  cout<<"size of Tf_eta0To0p5_pt50To100=  "<<Tf_eta0To0p5_pt50To100.size()<<endl;
  cout<<"size of Tf_eta0To0p5_pt100=  "<<Tf_eta0To0p5_pt100.size()<<endl;
  cout<<"size of Tf_eta0p5To1p0_pt10To50=  "<<Tf_eta0p5To1p0_pt10To50.size()<<endl;
  cout<<"size of Tf_eta0p5To1p0_pt50To100=  "<<Tf_eta0p5To1p0_pt50To100.size()<<endl;
  cout<<"size of Tf_eta0p5To1p0_pt100=  "<<Tf_eta0p5To1p0_pt100.size()<<endl;
  cout<<"size of Tf_eta1p0To2p5_pt10To50=  "<<Tf_eta1p0To2p5_pt10To50.size()<<endl;
  cout<<"size of Tf_eta1p0To2p5_pt50To100=  "<<Tf_eta1p0To2p5_pt50To100.size()<<endl;
  cout<<"size of Tf_eta1p0To2p5_pt100=  "<<Tf_eta1p0To2p5_pt100.size()<<endl;
  */

  //find avg and std dev of Tfs
  double avg_Tf_eta0To0p5_pt10To50 = findAvg(Tf_eta0To0p5_pt10To50); double stddev_Tf_eta0To0p5_pt10To50 = findSD(Tf_eta0To0p5_pt10To50);
  double avg_Tf_eta0To0p5_pt50To100 = findAvg(Tf_eta0To0p5_pt50To100); double stddev_Tf_eta0To0p5_pt50To100 = findSD(Tf_eta0To0p5_pt50To100);
  double avg_Tf_eta0To0p5_pt100 = findAvg(Tf_eta0To0p5_pt100); double stddev_Tf_eta0To0p5_pt100 = findSD(Tf_eta0To0p5_pt100);

  double avg_Tf_eta0p5To1p0_pt10To50 = findAvg(Tf_eta0p5To1p0_pt10To50); double stddev_Tf_eta0p5To1p0_pt10To50 = findSD(Tf_eta0p5To1p0_pt10To50);
  double avg_Tf_eta0p5To1p0_pt50To100 = findAvg(Tf_eta0p5To1p0_pt50To100); double stddev_Tf_eta0p5To1p0_pt50To100 = findSD(Tf_eta0p5To1p0_pt50To100);
  double avg_Tf_eta0p5To1p0_pt100 = findAvg(Tf_eta0p5To1p0_pt100); double stddev_Tf_eta0p5To1p0_pt100 = findSD(Tf_eta0p5To1p0_pt100);

  double avg_Tf_eta1p0To2p5_pt10To50 = findAvg(Tf_eta1p0To2p5_pt10To50); double stddev_Tf_eta1p0To2p5_pt10To50 = findSD(Tf_eta1p0To2p5_pt10To50);
  double avg_Tf_eta1p0To2p5_pt50To100 = findAvg(Tf_eta1p0To2p5_pt50To100); double stddev_Tf_eta1p0To2p5_pt50To100 = findSD(Tf_eta1p0To2p5_pt50To100);
  double avg_Tf_eta1p0To2p5_pt100 = findAvg(Tf_eta1p0To2p5_pt100); double stddev_Tf_eta1p0To2p5_pt100 = findSD(Tf_eta1p0To2p5_pt100);
  
  
  cout<<"avg_Tf_eta0To0p5_pt10To50=  "<<avg_Tf_eta0To0p5_pt10To50<<"   stddev_Tf_eta0To0p5_pt10To50=  "<<stddev_Tf_eta0To0p5_pt10To50<<endl;
  cout<<"avg_Tf_eta0To0p5_pt50To100=  "<<avg_Tf_eta0To0p5_pt50To100<<"   stddev_Tf_eta0To0p5_pt50To100=  "<<stddev_Tf_eta0To0p5_pt50To100<<endl;
  cout<<"avg_Tf_eta0To0p5_pt100=  "<<avg_Tf_eta0To0p5_pt100<<"   stddev_Tf_eta0To0p5_pt100=  "<<stddev_Tf_eta0To0p5_pt100<<endl;

  cout<<"avg_Tf_eta0p5To1p0_pt10To50=  "<<avg_Tf_eta0p5To1p0_pt10To50<<"   stddev_Tf_eta0p5To1p0_pt10To50=  "<<stddev_Tf_eta0p5To1p0_pt10To50<<endl;
  cout<<"avg_Tf_eta0p5To1p0_pt50To100=  "<<avg_Tf_eta0p5To1p0_pt50To100<<"   stddev_Tf_eta0p5To1p0_pt50To100=  "<<stddev_Tf_eta0p5To1p0_pt50To100<<endl;
  cout<<"avg_Tf_eta0p5To1p0_pt100=  "<<avg_Tf_eta0p5To1p0_pt100<<"   stddev_Tf_eta0p5To1p0_pt100=  "<<stddev_Tf_eta0p5To1p0_pt100<<endl;

  cout<<"avg_Tf_eta1p0To2p5_pt10To50=  "<<avg_Tf_eta1p0To2p5_pt10To50<<"   stddev_Tf_eta1p0To2p5_pt10To50=  "<<stddev_Tf_eta1p0To2p5_pt10To50<<endl;
  cout<<"avg_Tf_eta1p0To2p5_pt50To100=  "<<avg_Tf_eta1p0To2p5_pt50To100<<"   stddev_Tf_eta1p0To2p5_pt50To100=  "<<stddev_Tf_eta1p0To2p5_pt50To100<<endl;
  cout<<"avg_Tf_eta1p0To2p5_pt100=  "<<avg_Tf_eta1p0To2p5_pt100<<"   stddev_Tf_eta1p0To2p5_pt100=  "<<stddev_Tf_eta1p0To2p5_pt100<<endl;







  //applying these corrections to sig-only jets just to see effect (of regions and avg) 
  for(unsigned int ientry=0; ientry<fChain_s_gen->GetEntries(); ++ientry){
    fChain_s_gen->GetEntry(ientry);
    fChain_s_reco->GetEntry(ientry);

    genjets.clear();
    goodgenjets.clear();
    recojets.clear();
    goodrecojets.clear();
    matchedrecojets.clear();
    matchedgoodrecojets.clear();
    gen_matchreco_jets.clear();
    goodcorrmatchrecojets.clear();

    for(int a=0; a<nj_s_gen; a++){
      TLorentzVector gen; gen.SetPxPyPzE(jmox_s_gen[a], jmoy_s_gen[a], jmoz_s_gen[a], jene_s_gen[a]); gen.SetTheta(TMath::ACos(jcost_s_gen[a]));
      genjets.push_back(gen);

      //acceptance cuts
      if(gen.Pt()>10 && abs(gen.Theta())>0.1745 && abs(gen.Theta())<2.9671){
	goodgenjets.push_back(gen);
      }
    }
   
    //reco jets
    for(int a=0; a<nj_s_reco; a++){
      TLorentzVector reco; reco.SetPxPyPzE(jmox_s_reco[a], jmoy_s_reco[a], jmoz_s_reco[a], jene_s_reco[a]); reco.SetTheta(TMath::ACos(jcost_s_reco[a]));
      recojets.push_back(reco);                 
      //acceptance cuts
      if(reco.Pt()>10 && abs(reco.Theta())>0.1745 && abs(reco.Theta())<2.9671){
	goodrecojets.push_back(reco);
      }
    }

    //matching recojets to genjets in min dR    
    if(genjets.size()>0 && recojets.size()>0){
      for(int j=0; j<genjets.size(); j++){
	int recoindex = FindRecoMatch(recojets, genjets.at(j));
	double dR = find_delR(genjets.at(j), recojets.at(recoindex));
	if(dR<0.5){
	  matchedrecojets.push_back(recojets.at(recoindex));	
	  //double jec_ptRatio = genjets.at(j).Pt()/recojets.at(recoindex).Pt() ;
	  vector<TLorentzVector> vec; vec.push_back(genjets.at(j)); vec.push_back(recojets.at(recoindex));
	  gen_matchreco_jets.push_back(vec);
	}
      }
    }
  
    h_Nmatchedrecojets->Fill(matchedrecojets.size());
    for(int g=0; g<matchedrecojets.size(); g++){
      //acceptance cuts
      if(matchedrecojets.at(g).Pt()>10 && abs(matchedrecojets.at(g).Theta())>0.1745 && abs(matchedrecojets.at(g).Theta())<2.9671){
	matchedgoodrecojets.push_back(matchedrecojets.at(g));
      }
    }

    //applying JEC to matched reco jets  
    if(gen_matchreco_jets.size()>0){
      for(int k=0; k<gen_matchreco_jets.size(); k++){
	double genjEta = abs(gen_matchreco_jets[k][0].Eta()); double genjPt = gen_matchreco_jets[k][0].Pt();

	//eta [0, 0.5] & PT [10,50,100,>100]
	if(genjEta>=0 && genjEta<0.5){
	  if(genjPt>=10 && genjPt<50){	  
	    gen_matchreco_jets[k][1] = avg_Tf_eta0To0p5_pt10To50*gen_matchreco_jets[k][1] ;
	  }
	  if(genjPt>=50 && genjPt<100){
	    gen_matchreco_jets[k][1] = avg_Tf_eta0To0p5_pt50To100*gen_matchreco_jets[k][1];
	  }
	  if(genjPt>=100){
	    gen_matchreco_jets[k][1] = avg_Tf_eta0To0p5_pt100*gen_matchreco_jets[k][1];
	  }	
	}

	//eta [0.5, 1.0] & PT [10,50,100,>100]
	if(genjEta>=0.5 && genjEta<1.0){
	  if(genjPt>=10 && genjPt<50){
	    gen_matchreco_jets[k][1] = avg_Tf_eta0p5To1p0_pt10To50*gen_matchreco_jets[k][1];
	  }
	  if(genjPt>=50 && genjPt<100){
	    gen_matchreco_jets[k][1] = avg_Tf_eta0p5To1p0_pt50To100*gen_matchreco_jets[k][1];
	  }
	  if(genjPt>=100){
	    gen_matchreco_jets[k][1] = avg_Tf_eta0p5To1p0_pt100*gen_matchreco_jets[k][1];
	  }	
	}

	//eta [1.0, 2.5] & PT [10,50,100,>100]
	if(genjEta>=1.0 && genjEta<2.5){
	  if(genjPt>=10 && genjPt<50){
	    gen_matchreco_jets[k][1] = avg_Tf_eta1p0To2p5_pt10To50*gen_matchreco_jets[k][1] ;
	  }
	  if(genjPt>=50 && genjPt<100){
	    gen_matchreco_jets[k][1] = avg_Tf_eta1p0To2p5_pt50To100*gen_matchreco_jets[k][1] ;
	  }
	  if(genjPt>=100){
	    gen_matchreco_jets[k][1] = avg_Tf_eta1p0To2p5_pt100*gen_matchreco_jets[k][1] ;
	  }
	}
      }
    }

    //kinematics after JEC
    if(gen_matchreco_jets.size()>0){
      for(int a=0; a<gen_matchreco_jets.size(); a++){
	TLorentzVector vv = gen_matchreco_jets[a][1];
	h_cmrj_Pt->Fill(vv.Pt());
	h_cmrj_Eta->Fill(vv.Eta());
	h_cmrj_Phi->Fill(vv.Phi());
	h_cmrj_E->Fill(vv.E());     

	//acceptance cuts
	if(vv.Pt()>10 && abs(vv.Theta())>0.1745 && abs(vv.Theta())<2.9671){	
	  goodcorrmatchrecojets.push_back(vv);
	  h_cmgrj_Pt->Fill(vv.Pt());
	  h_cmgrj_Eta->Fill(vv.Eta());
	  h_cmgrj_Phi->Fill(vv.Phi());
	  h_cmgrj_E->Fill(vv.E());     

	}
      }
    }    

    //now plot inv mass after jec applied
    //also need to select good reco jets
    if(gen_matchreco_jets.size()>1){
      double mass_cmrj = (gen_matchreco_jets[0][1] + gen_matchreco_jets[1][1]).M();
      h_M_cmrj->Fill(mass_cmrj);
    }
    if(goodcorrmatchrecojets.size()>1){
      double mass_cmgrj = (goodcorrmatchrecojets.at(0) + goodcorrmatchrecojets.at(1)).M();
      h_M_cmgrj->Fill(mass_cmgrj);
    }




  }//event loop
  

  //------------------------------------------------------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------SIG+BIB--------------------------------------------------------------------------------------
  //----------------------------------------------------------------------event---------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //Loop over sig+BIB events and apply JEC to reco jets.

  Int_t           ndaughters_sb_gen[30];
  Float_t         daughters_PX_sb_gen[30][60];   
  Float_t         daughters_PY_sb_gen[30][60];   
  Float_t         daughters_PZ_sb_gen[30][60];   
  Float_t         daughters_M_sb_gen[30][60];   
  Float_t         daughters_E_sb_gen[30][60];   

  Int_t           ndaughters_sb_reco[30];
  Float_t         daughters_PX_sb_reco[30][60];   
  Float_t         daughters_PY_sb_reco[30][60];   
  Float_t         daughters_PZ_sb_reco[30][60];   
  Float_t         daughters_M_sb_reco[30][60];   
  Float_t         daughters_E_sb_reco[30][60];   

  Int_t           nj_sb_gen;
  Float_t         jmox_sb_gen[30];   //[njet]                                                                                                                                                              
  Float_t         jmoy_sb_gen[30];   //[njet]                                                                                                                                                               
  Float_t         jmoz_sb_gen[30];   //[njet]                                                                                                                                                               
  Float_t         jmas_sb_gen[30];   //[njet]                                                                                                                                                               
  Float_t         jene_sb_gen[30];   //[njet]                                                                                                                                                               
  Float_t         jcost_sb_gen[30];   //[njet]                                                                                                                                                              
  Float_t         jmom_sb_gen[30];   //[njet]                                                                                                                                                               

  Int_t           nj_sb_reco;
  Float_t         jmox_sb_reco[30];   //[njet]                                                                                                                                                              
  Float_t         jmoy_sb_reco[30];   //[njet]                                                                                                                                                              
  Float_t         jmoz_sb_reco[30];   //[njet]                                                                                                                                                              
  Float_t         jmas_sb_reco[30];   //[njet]                                                                                                                                                              
  Float_t         jene_sb_reco[30];   //[njet]                                                                                                                                                              
  Float_t         jcost_sb_reco[30];   //[njet]                                                                                                                                                             
  Float_t         jmom_sb_reco[30];   //[njet]                                                                                                                                                             

  TChain* fChain_sb_gen = new TChain("fChain_sb_gen");
  TChain* fChain_sb_reco = new TChain("fChain_sb_reco");
  fChain_sb_gen->Add("JetHistogramsHbb_sig_and_BIB_bib100_2mev_250ps_ak5_10kEvents.root/JetHistogramGenJetTuple");
  fChain_sb_reco->Add("JetHistogramsHbb_sig_and_BIB_bib100_2mev_250ps_ak5_10kEvents.root/JetHistogramRecoJetTuple");

  fChain_sb_gen->SetBranchAddress("ndaughters", &ndaughters_sb_gen, &b_ndaughters);
  fChain_sb_gen->SetBranchAddress("daughters_PX", daughters_PX_sb_gen, &b_daughters_PX);
  fChain_sb_gen->SetBranchAddress("daughters_PY", daughters_PY_sb_gen, &b_daughters_PY);
  fChain_sb_gen->SetBranchAddress("daughters_PZ", daughters_PZ_sb_gen, &b_daughters_PZ);
  fChain_sb_gen->SetBranchAddress("daughters_M", daughters_M_sb_gen, &b_daughters_M);
  fChain_sb_gen->SetBranchAddress("daughters_E", daughters_E_sb_gen, &b_daughters_E);

  fChain_sb_reco->SetBranchAddress("ndaughters", &ndaughters_sb_reco, &b_ndaughters);
  fChain_sb_reco->SetBranchAddress("daughters_PX", daughters_PX_sb_reco, &b_daughters_PX);
  fChain_sb_reco->SetBranchAddress("daughters_PY", daughters_PY_sb_reco, &b_daughters_PY);
  fChain_sb_reco->SetBranchAddress("daughters_PZ", daughters_PZ_sb_reco, &b_daughters_PZ);
  fChain_sb_reco->SetBranchAddress("daughters_M", daughters_M_sb_reco, &b_daughters_M);
  fChain_sb_reco->SetBranchAddress("daughters_E", daughters_E_sb_reco, &b_daughters_E);
    
  fChain_sb_gen->SetBranchAddress("nj", &nj_sb_gen, &b_njet);
  fChain_sb_gen->SetBranchAddress("jmox", jmox_sb_gen, &b_jmox);
  fChain_sb_gen->SetBranchAddress("jmoy", jmoy_sb_gen, &b_jmoy);
  fChain_sb_gen->SetBranchAddress("jmoz", jmoz_sb_gen, &b_jmoz);
  fChain_sb_gen->SetBranchAddress("jmas", jmas_sb_gen, &b_jmas);
  fChain_sb_gen->SetBranchAddress("jene", jene_sb_gen, &b_jene);
  fChain_sb_gen->SetBranchAddress("jcost", jcost_sb_gen, &b_jcost);
  fChain_sb_gen->SetBranchAddress("jmom", jmom_sb_gen, &b_jmom);
  

  fChain_sb_reco->SetBranchAddress("nj", &nj_sb_reco, &b_njet);
  fChain_sb_reco->SetBranchAddress("jmox", jmox_sb_reco, &b_jmox);
  fChain_sb_reco->SetBranchAddress("jmoy", jmoy_sb_reco, &b_jmoy);
  fChain_sb_reco->SetBranchAddress("jmoz", jmoz_sb_reco, &b_jmoz);
  fChain_sb_reco->SetBranchAddress("jmas", jmas_sb_reco, &b_jmas);
  fChain_sb_reco->SetBranchAddress("jene", jene_sb_reco, &b_jene);
  fChain_sb_reco->SetBranchAddress("jcost", jcost_sb_reco, &b_jcost);
  fChain_sb_reco->SetBranchAddress("jmom", jmom_sb_reco, &b_jmom);


  vector <TLorentzVector> genjets_sb;
  vector <TLorentzVector> goodgenjets_sb;
  vector <TLorentzVector> tightgenjets_sb;
  vector <TLorentzVector> recojets_sb;
  vector <TLorentzVector> goodrecojets_sb;
  vector <TLorentzVector> matchedrecojets_sb;
  vector <TLorentzVector> matchedgoodrecojets_sb;
  vector <TLorentzVector> goodcorrmatchrecojets_sb;
  vector <vector<TLorentzVector>> gen_matchreco_jets_sb;

  vector <double> pt_ratio_sb_reco;
  vector <double> e_ratio_sb_reco;
  vector <double> pt_ratio_mrj_sb;
  vector <double> e_ratio_mrj_sb;


  for(unsigned int ientry=0; ientry<fChain_sb_gen->GetEntries(); ++ientry){
    fChain_sb_gen->GetEntry(ientry);
    fChain_sb_reco->GetEntry(ientry);

    genjets_sb.clear();
    goodgenjets_sb.clear();
    tightgenjets_sb.clear();
    recojets_sb.clear();
    goodrecojets_sb.clear();
    matchedrecojets_sb.clear();
    matchedgoodrecojets_sb.clear();
    gen_matchreco_jets_sb.clear();
    goodcorrmatchrecojets_sb.clear();

    pt_ratio_sb_reco.clear();
    e_ratio_sb_reco.clear();
    pt_ratio_mrj_sb.clear();
    e_ratio_mrj_sb.clear();


    for(int a=0; a<nj_sb_gen; a++){
      TLorentzVector gen; gen.SetPxPyPzE(jmox_sb_gen[a], jmoy_sb_gen[a], jmoz_sb_gen[a], jene_sb_gen[a]); gen.SetTheta(TMath::ACos(jcost_sb_gen[a]));
      genjets_sb.push_back(gen);
      h_sb_gj_Pt->Fill(gen.Pt());
      h_sb_gj_Eta->Fill(gen.Eta());
      h_sb_gj_Phi->Fill(gen.Phi());
      h_sb_gj_E->Fill(gen.E());

      h_sb_Ndau_gj->Fill(ndaughters_sb_gen[a]);

      //acceptance cuts
      if(gen.Pt()>10 && abs(gen.Theta())>0.1745 && abs(gen.Theta())<2.9671){
	goodgenjets_sb.push_back(gen);
	h_sb_ggj_Pt->Fill(gen.Pt());
	h_sb_ggj_Eta->Fill(gen.Eta());
	h_sb_ggj_Phi->Fill(gen.Phi());
	h_sb_ggj_E->Fill(gen.E());	

	h_sb_Ndau_ggj->Fill(ndaughters_sb_gen[a]);
	if(gen.Pt()>20 && abs(gen.Theta())>0.1745 && abs(gen.Theta())<2.9671){
	  tightgenjets_sb.push_back(gen);
	}
      }
    }
    h_sb_Ngenjets->Fill(genjets_sb.size());
    h_sb_Ngoodgenjets->Fill(goodgenjets_sb.size());
    h_sb_Ntightgenjets->Fill(tightgenjets_sb.size());
    
    for(int a=0; a<nj_sb_reco; a++){
      TLorentzVector reco; reco.SetPxPyPzE(jmox_sb_reco[a], jmoy_sb_reco[a], jmoz_sb_reco[a], jene_sb_reco[a]); reco.SetTheta(TMath::ACos(jcost_sb_reco[a]));
      recojets_sb.push_back(reco);
      h_sb_rj_Pt->Fill(reco.Pt());
      h_sb_rj_Eta->Fill(reco.Eta());
      h_sb_rj_Phi->Fill(reco.Phi());
      h_sb_rj_E->Fill(reco.E());
      
      h_sb_Ndau_rj->Fill(ndaughters_sb_reco[a]);

      //daughters
      for(int b=0; b<ndaughters_sb_reco[a]; b++){
	TLorentzVector dau; dau.SetPxPyPzE(daughters_PX_sb_reco[a][b],daughters_PY_sb_reco[a][b],daughters_PZ_sb_reco[a][b],daughters_E_sb_reco[a][b]);
	double pt_ratio = dau.Pt()/reco.Pt();
	pt_ratio_sb_reco.push_back(pt_ratio);
	double e_ratio = dau.E()/reco.E();
	e_ratio_sb_reco.push_back(e_ratio);	
      }

      //acceptance cuts
      if(reco.Pt()>10 && abs(reco.Theta())>0.1745 && abs(reco.Theta())<2.9671){
	goodrecojets_sb.push_back(reco);
	h_sb_grj_Pt->Fill(reco.Pt());
	h_sb_grj_Eta->Fill(reco.Eta());
	h_sb_grj_Phi->Fill(reco.Phi());
	h_sb_grj_E->Fill(reco.E());

	h_sb_Ndau_grj->Fill(ndaughters_sb_reco[a]);	
      
	//find max pt_ratio and e_ratio
	if(pt_ratio_sb_reco.size()>0){
	  double max_pt_ratio = pt_ratio_sb_reco.at(findmax(pt_ratio_sb_reco));
	  double max_e_ratio = e_ratio_sb_reco.at(findmax(e_ratio_sb_reco));
	  h_sb_maxPTratio_grj->Fill(max_pt_ratio);
	  h_sb_maxEratio_grj->Fill(max_e_ratio);
	}
      }  
             
    }    

    h_sb_Nrecojets->Fill(recojets_sb.size());
    h_sb_Ngoodrecojets->Fill(goodrecojets_sb.size());


    int Ndau_mrj_sb; double max_pt_ratio_mrj_sb; double max_e_ratio_mrj_sb;
    //matching recojets to genjets in min dR    
    if(genjets_sb.size()>0 && recojets_sb.size()>0){
      for(int j=0; j<genjets_sb.size(); j++){
	int recoindex = FindRecoMatch(recojets_sb, genjets_sb.at(j));
	double dR = find_delR(genjets_sb.at(j), recojets_sb.at(recoindex));
	if(dR<0.5){
	  matchedrecojets_sb.push_back(recojets_sb.at(recoindex));	

	  //daughters                                                                                                                                                                                       
          Ndau_mrj_sb = ndaughters_sb_reco[recoindex];
          h_sb_Ndau_mrj->Fill(Ndau_mrj_sb);

	  //pt_ratio 
	  for(int k=0; k<ndaughters_sb_reco[recoindex]; k++){
	    TLorentzVector dau; dau.SetPxPyPzE(daughters_PX_sb_reco[recoindex][k],daughters_PY_sb_reco[recoindex][k],daughters_PZ_sb_reco[recoindex][k],daughters_E_sb_reco[recoindex][k]);
	    double dau_pt_ratio_mrj = dau.Pt()/recojets_sb.at(recoindex).Pt();
	    double dau_e_ratio_mrj = dau.E()/recojets_sb.at(recoindex).E();
	    pt_ratio_mrj_sb.push_back(dau_pt_ratio_mrj);
	    e_ratio_mrj_sb.push_back(dau_e_ratio_mrj);
	  }
	  //findmax
	  max_pt_ratio_mrj_sb = pt_ratio_mrj_sb.at(findmax(pt_ratio_mrj_sb));
	  max_e_ratio_mrj_sb = e_ratio_mrj_sb.at(findmax(e_ratio_mrj_sb));
	  h_sb_maxPTratio_mrj->Fill(max_pt_ratio_mrj_sb);
	  h_sb_maxEratio_mrj->Fill(max_e_ratio_mrj_sb);

	  vector<TLorentzVector> vec; vec.push_back(genjets_sb.at(j)); vec.push_back(recojets_sb.at(recoindex));
	  gen_matchreco_jets_sb.push_back(vec);
	}
      }
    }

    h_sb_Ngen_matchreco_jets->Fill(gen_matchreco_jets_sb.size());
    h_sb_Nmatchedrecojets->Fill(matchedrecojets_sb.size());
    for(int g=0; g<matchedrecojets_sb.size(); g++){
      h_sb_mrj_Pt->Fill(matchedrecojets_sb.at(g).Pt());
      h_sb_mrj_Eta->Fill(matchedrecojets_sb.at(g).Eta());
      h_sb_mrj_Phi->Fill(matchedrecojets_sb.at(g).Phi());
      h_sb_mrj_E->Fill(matchedrecojets_sb.at(g).E());
      //acceptance cuts
      if(matchedrecojets_sb.at(g).Pt()>10 && abs(matchedrecojets_sb.at(g).Theta())>0.1745 && abs(matchedrecojets_sb.at(g).Theta())<2.9671){
	matchedgoodrecojets_sb.push_back(matchedrecojets_sb.at(g));

	h_sb_Ndau_mgrj->Fill(Ndau_mrj_sb);
	h_sb_maxPTratio_mgrj->Fill(max_pt_ratio_mrj_sb);
	h_sb_maxEratio_mgrj->Fill(max_e_ratio_mrj_sb);
      }
    }
    h_sb_Nmatchedgoodrecojets->Fill(matchedgoodrecojets_sb.size()); 


    //inv_mass of di-jet pair 1) two leading jets
    //reco jets
    if(recojets_sb.size()>1){
      double inv_m1 = (recojets_sb.at(0) + recojets_sb.at(1)).M();
      h_sb_M_j0j1_rj->Fill(inv_m1);
    }
    if(goodrecojets_sb.size()>1){
      double inv_m2 = (goodrecojets_sb.at(0) + goodrecojets_sb.at(1)).M();
      h_sb_M_j0j1_grj->Fill(inv_m2);
    }
    if(matchedrecojets_sb.size()>1){
      double inv_m3 = (matchedrecojets_sb.at(0) + matchedrecojets_sb.at(1)).M();
      h_sb_M_j0j1_mrj->Fill(inv_m3);
    }
    if(matchedgoodrecojets_sb.size()>1){
      double inv_m4 = (matchedgoodrecojets_sb.at(0) + matchedgoodrecojets_sb.at(1)).M();
      h_sb_M_j0j1_mgrj->Fill(inv_m4);
    }


    //gen jets
    if(genjets_sb.size()>1){
      double inv_m_g = (genjets_sb.at(0) + genjets_sb.at(1)).M();
      h_sb_M_j0j1_gj->Fill(inv_m_g);
    }
    if(goodgenjets_sb.size()>1){
      double inv_m_gg = (goodgenjets_sb.at(0) + goodgenjets_sb.at(1)).M();
      h_sb_M_j0j1_ggj->Fill(inv_m_gg);
    }
    if(tightgenjets_sb.size()>1){
      double inv_m_tg = (tightgenjets_sb.at(0) + tightgenjets_sb.at(1)).M();
      h_sb_M_j0j1_tgj->Fill(inv_m_tg);
    }
  

    //applying JEC to matched reco jets  
    if(gen_matchreco_jets_sb.size()>0){
      for(int k=0; k<gen_matchreco_jets_sb.size(); k++){
	double genjEta = abs(gen_matchreco_jets_sb[k][0].Eta()); double genjPt = gen_matchreco_jets_sb[k][0].Pt();

	//eta [0, 0.5] & PT [10,50,100,>100]
	if(genjEta>=0 && genjEta<0.5){
	  if(genjPt>=10 && genjPt<50){	  
	    gen_matchreco_jets_sb[k][1] = avg_Tf_eta0To0p5_pt10To50*gen_matchreco_jets_sb[k][1] ;
	  }
	  if(genjPt>=50 && genjPt<100){
	    gen_matchreco_jets_sb[k][1] = avg_Tf_eta0To0p5_pt50To100*gen_matchreco_jets_sb[k][1];
	  }
	  if(genjPt>=100){
	    gen_matchreco_jets_sb[k][1] = avg_Tf_eta0To0p5_pt100*gen_matchreco_jets_sb[k][1];
	  }	
	}

	//eta [0.5, 1.0] & PT [10,50,100,>100]
	if(genjEta>=0.5 && genjEta<1.0){
	  if(genjPt>=10 && genjPt<50){
	    gen_matchreco_jets_sb[k][1] = avg_Tf_eta0p5To1p0_pt10To50*gen_matchreco_jets_sb[k][1];
	  }
	  if(genjPt>=50 && genjPt<100){
	    gen_matchreco_jets_sb[k][1] = avg_Tf_eta0p5To1p0_pt50To100*gen_matchreco_jets_sb[k][1];
	  }
	  if(genjPt>=100){
	    gen_matchreco_jets_sb[k][1] = avg_Tf_eta0p5To1p0_pt100*gen_matchreco_jets_sb[k][1];
	  }	
	}

	//eta [1.0, 2.5] & PT [10,50,100,>100]
	if(genjEta>=1.0 && genjEta<2.5){
	  if(genjPt>=10 && genjPt<50){
	    gen_matchreco_jets_sb[k][1] = avg_Tf_eta1p0To2p5_pt10To50*gen_matchreco_jets_sb[k][1] ;
	  }
	  if(genjPt>=50 && genjPt<100){
	    gen_matchreco_jets_sb[k][1] = avg_Tf_eta1p0To2p5_pt50To100*gen_matchreco_jets_sb[k][1] ;
	  }
	  if(genjPt>=100){
	    gen_matchreco_jets_sb[k][1] = avg_Tf_eta1p0To2p5_pt100*gen_matchreco_jets_sb[k][1] ;
	  }
	}
      }
    }

    //kinematics after JEC
    if(gen_matchreco_jets_sb.size()>0){
      for(int a=0; a<gen_matchreco_jets_sb.size(); a++){
	TLorentzVector vv = gen_matchreco_jets_sb[a][1];
	h_sb_cmrj_Pt->Fill(vv.Pt());
	h_sb_cmrj_Eta->Fill(vv.Eta());
	h_sb_cmrj_Phi->Fill(vv.Phi());
	h_sb_cmrj_E->Fill(vv.E());     

	//acceptance cuts
	if(vv.Pt()>10 && abs(vv.Theta())>0.1745 && abs(vv.Theta())<2.9671){	
	  goodcorrmatchrecojets_sb.push_back(vv);
	  h_sb_cmgrj_Pt->Fill(vv.Pt());
	  h_sb_cmgrj_Eta->Fill(vv.Eta());
	  h_sb_cmgrj_Phi->Fill(vv.Phi());
	  h_sb_cmgrj_E->Fill(vv.E());     

	}
      }
    }
    h_sb_Ncmrj->Fill(gen_matchreco_jets_sb.size());
    h_sb_Ncmgrj->Fill(goodcorrmatchrecojets_sb.size());


    //now plot inv mass after jec applied
    //also need to select good reco jets

    if(gen_matchreco_jets_sb.size()>1){
      double mass_cmrj = (gen_matchreco_jets_sb[0][1] + gen_matchreco_jets_sb[1][1]).M();
      h_sb_M_cmrj->Fill(mass_cmrj);
    }

    if(goodcorrmatchrecojets_sb.size()>1){
      double mass_cmgrj = (goodcorrmatchrecojets_sb.at(0) + goodcorrmatchrecojets_sb.at(1)).M();
      h_sb_M_cmgrj->Fill(mass_cmgrj);
    }


  }//sb events loop






  //---------------------------------------------------------------------------------
  /*
  TCanvas *c1 = new TCanvas("c1","c1",800,600);


  //Jet
  h3_theta_phi_ene_j->Draw("lego");
  c1->SaveAs("plots_2mev_250ps_jpar/Theta_Phi_E_Jet.pdf");
  h_z_jet->Draw();
  c1->SaveAs("plots_2mev_250ps_jpar/Z_Jet.pdf");
  h_R_jet->Draw();
  c1->SaveAs("plots_2mev_250ps_jpar/R_Jet.pdf");
  h2_z_R_jet->Draw("colz");
  c1->SaveAs("plots_2mev_250ps_jpar/Z_R_Jet.pdf");
  h3_z_R_E_jet->Draw("lego");
  c1->SaveAs("plots_2mev_250ps_jpar/Z_R_E_Jet.pdf");

  //jpar
  h3_theta_phi_ene_jpar->Draw("lego");
  c1->SaveAs("plots_2mev_250ps_jpar/Theta_Phi_E_jpar.pdf");
  h2_z_R_jpar->Draw("colz");
  c1->SaveAs("plots_2mev_250ps_jpar/Z_R_jpar.pdf");
  h3_z_R_E_jpar->Draw("lego");
  c1->SaveAs("plots_2mev_250ps_jpar/Z_R_E_jpar.pdf");
  h_pt_ratio_max->Draw();
  c1->SaveAs("plots_2mev_250ps_jpar/PTmax_ratio_jpar.pdf");
  h_ene_ratio_max->Draw();
  c1->SaveAs("plots_2mev_250ps_jpar/Emax_ratio_jpar.pdf");
  h_min_z_jpar->Draw();
  c1->SaveAs("plots_2mev_250ps_jpar/Zmin_jpar.pdf");
  h_max_z_jpar->Draw();
  c1->SaveAs("plots_2mev_250ps_jpar/Zmax_jpar.pdf");
  h_min_R_jpar->Draw();
  c1->SaveAs("plots_2mev_250ps_jpar/Rmin_jpar.pdf");
  h_max_R_jpar->Draw();
  c1->SaveAs("plots_2mev_250ps_jpar/Rmin_jpar.pdf");
  h_delZ_jpar->Draw();
  c1->SaveAs("plots_2mev_250ps_jpar/delZ_jpar.pdf");
  h_delR_jpar->Draw();
  c1->SaveAs("plots_2mev_250ps_jpar/delR_jpar.pdf");
*/



  TFile* f=new TFile("plots_Hbb_bib100_2mev_250ps_ak5_10kEvents_s.root","RECREATE");
  if ( f->IsOpen() ) cout << "File opened successfully" << endl;

  //sig only
  h_Ngenjets->Write();
  h_Nrecojets->Write();
  h_Ngoodgenjets->Write();
  h_Ntightgenjets->Write();
  h_Ngoodrecojets->Write();
  h_Nmatchedrecojets->Write();
  h_Nmatchedgoodrecojets->Write();
  h_Ndau_gj->Write();
  h_Ndau_ggj->Write();
  h_Ndau_rj->Write();
  h_Ndau_grj->Write();
  h_Ndau_mrj->Write();
  h_Ndau_mgrj->Write();

  h_gj_Pt->Write();
  h_gj_Eta->Write();
  h_gj_Phi->Write();
  h_gj_E->Write();
  h_rj_Pt->Write();
  h_rj_Eta->Write();
  h_rj_Phi->Write();
  h_rj_E->Write();
  h_ggj_Pt->Write();
  h_ggj_Eta->Write();
  h_ggj_Phi->Write();
  h_ggj_E->Write();
  h_grj_Pt->Write();
  h_grj_Eta->Write();
  h_grj_Phi->Write();
  h_grj_E->Write();
  h_mrj_Pt->Write();
  h_mrj_Eta->Write();
  h_mrj_Phi->Write();
  h_mrj_E->Write();

  h_M_j0j1_gj->Write();
  h_M_j0j1_ggj->Write();
  h_M_j0j1_tgj->Write();
  h_M_j0j1_rj->Write();
  h_M_j0j1_grj->Write();
  h_M_j0j1_mrj->Write();
  h_M_j0j1_mgrj->Write();

  h_M_jj_rj->Write();
  h_M_jj_grj->Write();

  h_cmrj_Pt->Write();
  h_cmrj_Eta->Write();
  h_cmrj_Phi->Write();
  h_cmrj_E->Write();
  h_cmgrj_Pt->Write();
  h_cmgrj_Eta->Write();
  h_cmgrj_Phi->Write();
  h_cmgrj_E->Write();

  h_M_cmrj->Write();
  h_M_cmgrj->Write();

  h_maxPTratio_mrj->Write();
  h_maxEratio_mrj->Write();
  h_maxPTratio_mgrj->Write();
  h_maxEratio_mgrj->Write();
  
  h_maxPTratio_grj->Write();
  h_maxEratio_grj->Write();


  //sig and BIB
  h_sb_Ngenjets->Write();
  h_sb_Nrecojets->Write();
  h_sb_Ngoodgenjets->Write();
  h_sb_Ntightgenjets->Write();
  h_sb_Ngoodrecojets->Write();
  h_sb_Nmatchedrecojets->Write();
  h_sb_Nmatchedgoodrecojets->Write();
  h_sb_Ngen_matchreco_jets->Write();
  h_sb_Ndau_gj->Write();
  h_sb_Ndau_ggj->Write();
  h_sb_Ndau_rj->Write();
  h_sb_Ndau_grj->Write();
  h_sb_Ndau_mrj->Write();
  h_sb_Ndau_mgrj->Write();

  h_sb_gj_Pt->Write();
  h_sb_gj_Eta->Write();
  h_sb_gj_Phi->Write();
  h_sb_gj_E->Write();
  h_sb_rj_Pt->Write();
  h_sb_rj_Eta->Write();
  h_sb_rj_Phi->Write();
  h_sb_rj_E->Write();
  h_sb_ggj_Pt->Write();
  h_sb_ggj_Eta->Write();
  h_sb_ggj_Phi->Write();
  h_sb_ggj_E->Write();
  h_sb_grj_Pt->Write();
  h_sb_grj_Eta->Write();
  h_sb_grj_Phi->Write();
  h_sb_grj_E->Write();
  h_sb_mrj_Pt->Write();
  h_sb_mrj_Eta->Write();
  h_sb_mrj_Phi->Write();
  h_sb_mrj_E->Write();

  h_sb_Ncmrj->Write();
  h_sb_Ncmgrj->Write();

  h_sb_cmrj_Pt->Write();
  h_sb_cmrj_Eta->Write();
  h_sb_cmrj_Phi->Write();
  h_sb_cmrj_E->Write();
  h_sb_cmgrj_Pt->Write();
  h_sb_cmgrj_Eta->Write();
  h_sb_cmgrj_Phi->Write();
  h_sb_cmgrj_E->Write();

  h_sb_M_j0j1_gj->Write();
  h_sb_M_j0j1_ggj->Write();
  h_sb_M_j0j1_tgj->Write();
  h_sb_M_j0j1_rj->Write();
  h_sb_M_j0j1_grj->Write();
  h_sb_M_j0j1_mrj->Write();
  h_sb_M_j0j1_mgrj->Write();

  h_sb_M_cmrj->Write();
  h_sb_M_cmgrj->Write();

  h_sb_maxPTratio_mrj->Write();
  h_sb_maxEratio_mrj->Write();
  h_sb_maxPTratio_mgrj->Write();
  h_sb_maxEratio_mgrj->Write();

  h_sb_maxPTratio_grj->Write();
  h_sb_maxEratio_grj->Write();

  //Jet
  /*  h2_theta_phi_j->Write();
  h3_theta_phi_ene_j->Write();
  h_t_jet->Write();
  h_z_jet->Write();
  h_R_jet->Write();
  h2_z_R_jet->Write();
  h3_z_R_E_jet->Write();
  h3_z_R_t_jet->Write();
  h_Ndau_jet->Write();

  //jpar
  h2_theta_phi_jpar->Write();
  h3_theta_phi_ene_jpar->Write();
  h2_z_R_jpar->Write();
  h3_z_R_E_jpar->Write();
  h_pt_ratio_max->Write();
  h_ene_ratio_max->Write();
  h_min_z_jpar->Write();
  h_max_z_jpar->Write();
  h_min_R_jpar->Write();
  h_max_R_jpar->Write();
  h_delZ_jpar->Write();
  h_delR_jpar->Write();
  h2_delZ_delR_jpar->Write();
  h_delT_z_jpar->Write();
  h_delT_R_jpar->Write();
  h_delR_delZg05->Write();
  h_delR_delZl05->Write();
  h_Ndau_delZg05->Write();
  h_Ndau_delZl05->Write();

  h_Nj->Write();
  h_jj_mass->Write();
  h_jj_p->Write();
  h_Ngj->Write();
  h_Ngjs->Write();
  h_dijet_M->Write();
  h_dijet_Pt->Write();
  */

}
