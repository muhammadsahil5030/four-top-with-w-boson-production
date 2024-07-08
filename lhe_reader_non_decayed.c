//&>/dev/null;x="${0%.*}";[ ! "$x" -ot "$0" ]||(rm -f "$x";g++ -o "$x" "$0" -I`root-config --incdir` `root-config --libs`);exit

// Build: g++ lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
// OR, g++ -std=c++11 lhe_reader_non_decayed.c -o lhe_reader_non_decayed -I`root-config --incdir` `root-config --libs`
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
//#include </home/muhammad/root6/build/include/TLorentzVector.h>
#include <TLorentzVector.h>
using namespace std;

// Pour une description du format leshouches
// hep-ph/0609017
//
// pour chaque evenement
// une ligne générale : NbPart idprocess poids scale alpha_em alpha_s
// pour chaque particule : id status mere1 mere2 couleur1 couleur2 px py pz E m lifetime spin  
int main(int argc, char **argv) {

  if (argc != 2) {
    cout << "Missing argument. Expected LHE filename without '.lhe'"<< endl;
    exit(-5);
  }
  float DeltaR(float eta1, float eta2, float phi1, float phi2);
  float DeltaPhi(float phi1, float phi2);

  string basename=argv[1];

  string lhefname = basename+".lhe";
  string rootfname = basename+".root";

  string tt;
  int event=0;
  int npart,idprocess;
  double weight,scale,alpha_em,alpha_s;

  TH1::SetDefaultSumw2(true);
  
//--------------------------*insert here histogram definition*---------------------------//
//

TH1F* hPt_up = new TH1F("hPt_up", "Pt up", 50, 0, 400);
TH1F* hEta_up = new TH1F("hEta_up","Eta up", 50, -6, 6);
TH1F* hPhi_up = new TH1F("hPhi_up", "Phi up", 50, -3.14, 3.14);

TH1F* hPt_down = new TH1F("hPt_down", "Pt down", 50, 0, 400);
TH1F* hEta_down = new TH1F("hEta_down","Eta down", 50, -6, 6);
TH1F* hPhi_down = new TH1F("hPhi_down", "Phi down", 50, -3.14, 3.14);

TH1F* hPt_charm = new TH1F("hPt_charm", "Pt charm", 50, 0, 400);
TH1F* hEta_charm = new TH1F("hEta_charm","Eta charm", 50, -6, 6);
TH1F* hPhi_charm = new TH1F("hPhi_charm", "Phi charm", 50, -3.14, 3.14);

TH1F* hPt_strange = new TH1F("hPt_strange", "Pt strange", 50, 0, 400);
TH1F* hEta_strange = new TH1F("hEta_strange","Eta strange", 50, -6, 6);
TH1F* hPhi_strange = new TH1F("hPhi_strange", "Phi strange", 50, -3.14, 3.14);

TH1F* hPt_b = new TH1F("hPt_b", "Pt b-quark", 50, 0, 400);
TH1F* hEta_b = new TH1F("hEta_b","Eta b-quark", 50, -6, 6);
TH1F* hPhi_b = new TH1F("hPhi_b", "Phi b-quark", 50, -3.14, 3.14);

  TH1F** hPt_bq = new TH1F*[4];
  TH1F** hEta_bq = new TH1F*[4];
  TH1F** hPhi_bq = new TH1F*[4];
  for (int bquark=0 ; bquark<4 ; bquark++)
  {
    ostringstream oss;
    oss << "bquark" << bquark ;
    string ptstr = oss.str()+"_pt";
    hPt_bq[bquark] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,400) ;
    string etastr=oss.str()+"_eta";
    hEta_bq[bquark] = new TH1F(etastr.c_str(),etastr.c_str(),50,-2.5, 2.5) ;
    string phistr=oss.str()+"_phi";
    hPhi_bq[bquark] = new TH1F(phistr.c_str(),phistr.c_str(),50,-3.14, 3.14) ;
  }
  
  TH1F** hPt_lq1 = new TH1F*[5];
  TH1F** hEta_lq1 = new TH1F*[5];
  TH1F** hPhi_lq1 = new TH1F*[5];
  for (int lquark1=0 ; lquark1<5 ; lquark1++)
  {
    ostringstream oss;
    oss << "up_quark" << lquark1 ;
    string ptstr = oss.str()+"_pt";
    hPt_lq1[lquark1] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,400) ;
    string etastr=oss.str()+"_eta";
    hEta_lq1[lquark1] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_lq1[lquark1] = new TH1F(phistr.c_str(),phistr.c_str(),50,-5,5) ;
  }
  
  TH1F** hPt_lq2 = new TH1F*[5];
  TH1F** hEta_lq2 = new TH1F*[5];
  TH1F** hPhi_lq2 = new TH1F*[5];
  for (int lquark2=0 ; lquark2<5 ; lquark2++)
  {
    ostringstream oss;
    oss << "down_quark" << lquark2 ;
    string ptstr = oss.str()+"_pt";
    hPt_lq2[lquark2] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,400) ;
    string etastr=oss.str()+"_eta";
    hEta_lq2[lquark2] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_lq2[lquark2] = new TH1F(phistr.c_str(),phistr.c_str(),50,-5,5) ;
  }
  
  TH1F** hPt_lq3 = new TH1F*[5];
  TH1F** hEta_lq3 = new TH1F*[5];
  TH1F** hPhi_lq3 = new TH1F*[5];
  for (int lquark3=0 ; lquark3<5 ; lquark3++)
  {
    ostringstream oss;
    oss << "charm_quark" << lquark3;
    string ptstr = oss.str()+"_pt";
    hPt_lq3[lquark3] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,400) ;
    string etastr=oss.str()+"_eta";
    hEta_lq3[lquark3] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_lq3[lquark3] = new TH1F(phistr.c_str(),phistr.c_str(),50,-5,5) ;
  }
  
  TH1F** hPt_lq4 = new TH1F*[5];
  TH1F** hEta_lq4 = new TH1F*[5];
  TH1F** hPhi_lq4 = new TH1F*[5];
  for (int lquark4=0 ; lquark4<5 ; lquark4++)
  {
    ostringstream oss;
    oss << "strange_quark" << lquark4 ;
    string ptstr = oss.str()+"_pt";
    hPt_lq4[lquark4] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,400) ;
    string etastr=oss.str()+"_eta";
    hEta_lq4[lquark4] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_lq4[lquark4] = new TH1F(phistr.c_str(),phistr.c_str(),50,-5,5) ;
  }
    
  TH1F** hPt_wboson = new TH1F*[5];
  TH1F** hEta_wboson = new TH1F*[5];
  TH1F** hPhi_wboson = new TH1F*[5];
  TH1F** hM_wboson = new TH1F*[5];
  for (int wboson=0 ; wboson<5 ; wboson++)
  {
    ostringstream oss;
    oss << "wboson" << wboson ;
    string ptstr = oss.str()+"_pt";
    hPt_wboson[wboson] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,300) ;
    string etastr=oss.str()+"_eta";
    hEta_wboson[wboson] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_wboson[wboson] = new TH1F(phistr.c_str(),phistr.c_str(),50,-3.16,3.16) ;
    string massstr=oss.str()+"_mass";
    hM_wboson[wboson] = new TH1F(massstr.c_str(),massstr.c_str(),50,60.0,100.0) ;
  }
  
  
  TH1F** hPt_top = new TH1F*[5];
  TH1F** hEta_top = new TH1F*[5];
  TH1F** hPhi_top = new TH1F*[5];
  TH1F** hM_top = new TH1F*[5];
  for (int t=0 ; t<5 ; t++) {
    ostringstream oss;
    oss << "top_quark" << t ;
    string ptstr = oss.str()+"_pt";
    hPt_top[t] = new TH1F(ptstr.c_str(),ptstr.c_str(),50,0,600) ;
    string etastr=oss.str()+"_eta";
    hEta_top[t] = new TH1F(etastr.c_str(),etastr.c_str(),50,-6,6) ;
    string phistr=oss.str()+"_phi";
    hPhi_top[t] = new TH1F(phistr.c_str(),phistr.c_str(),50,-3.5,3.5) ;
    string massstr=oss.str()+"_mass";
    hM_top[t] = new TH1F(massstr.c_str(),massstr.c_str(),50,150.0,190.0) ;
  } 
//-------------------------------*end histogram definition*-----------------------------------//
//
cout<<"I am here: "<<endl;

  int nlept=0, nsemi=0, nhadr=0;
  ifstream ff(lhefname.c_str(),ios::in); //ouverture du fichier .lhe
  //ifstream ff("test.lhe",ios::in); //ouverture du fichier .lhe
  //ifstream ff("/home/cms/perries/madgraph-newphysics/s1/madevent/Events/zp4000_unweighted_events.lhe",ios::in); //ouverture du fichier .lhe
  //ifstream ff("/home/cms/perries/madgraph-newphysics/QCD/madevent/Events/qcd_unweighted_events.lhe",ios::in); //ouverture du fichier .lhe
  int negativeWeight = 0 , n_bq=0;
  long long line = 0;
  while(!ff.eof()) {
    std::stringstream buffer;
    ff>>tt;
    buffer << tt << std::endl;
    line++;

    if(tt=="<event>") {
      ff>>npart>>idprocess>>weight>>scale>>alpha_em>>alpha_s; //event definition
      buffer << npart << " " << idprocess << " " << weight << " " << scale << " " << alpha_em << " " << alpha_s << std::endl;
      line++;
      event++;
      if (weight < 0) {
        negativeWeight++;
        weight = -1;
      } else {
        weight = 1;
      }
      /*weight = 1.;*/

      if (event%1==0) cout << "reading event "<< event << endl;
      int lmin=-1, lmax=-1, bmin=-1, met1=-1, met2=-1, bmax=-1, met=-1, topLep,topbJet, tbarLep, tbarbJet, top_jet, tbar_jet;
      int Part_Id, Moth_Id, Part_absId, Moth_absId;
int n_topjets=0, n_tbarjets=0, n_topbjets=0, n_bJets=0, n_tbar_nu=0, n_top_nu=0, n_wlep=0;
int n_bq=0, n_upq=0, n_downq=0, n_charmq=0, n_strangeq=0;

int  bj[2]={-1,-1};
int lp[2]={-1,-1};

      int lq1[5]={-1,-1,-1,-1,-1};
      int lq2[5]={-1,-1,-1,-1,-1};
      int lq3[5]={-1,-1,-1,-1,-1};
      int lq4[5]={-1,-1,-1,-1,-1};
      int bq[4]={-1,-1,-1,-1};
      
      int top=-1,topbar=-1,zprime=-1;
      int *Id      = new int[npart+1];
      int *Status  = new int[npart+1];
      int *Mother1 = new int[npart+1];
      int *Mother2 = new int[npart+1];
      int *Color1  = new int[npart+1];
      int *Color2  = new int[npart+1];
      double *px = new double[npart+1];
      double *py = new double[npart+1];
      double *pz = new double[npart+1];
      double *E = new double[npart+1];
      double *m = new double[npart+1];
      double *lifetime = new double[npart+1];
      double *spin = new double[npart+1];
      
      TLorentzVector **v = new TLorentzVector*[npart+1];
      TLorentzVector v_upq, v_up1, v_up2, v_up3, v_up4, v_up5, v_charmq, v_charm1, v_charm2, v_charm3, v_charm4, v_charm5;
      TLorentzVector v_downq, v_down1, v_down2, v_down3, v_down4, v_down5, v_strangeq, v_strange1, v_strange2, v_strange3, v_strange4, v_strange5;
      TLorentzVector v_lep1, v_lep2, v_lep3, v_lep4, v_lep5;
      TLorentzVector v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4, v_w_boson5;
      TLorentzVector v_bq, v_bq1, v_bq2, v_bq3, v_bq4, v_bqtotal;
      TLorentzVector v_top1, v_top2, v_top3, v_top4, v_top5;
      TLorentzVector v_top_lep, v_top_alep, v_lep_nu, v_top_bJet, v_tbar_lep, v_top_jet, v_tbar_jet, v_tbar_nu, v_top_nu, v_wminus_lep;
      
      // in lhe first line is number 1, so fill unused array [0] with a crazy value;
      Id[0]= -99999;
      Status[0]= -99999;
      Mother1[0]= -99999;
      Mother2[0]= -99999;
      Color1[0]= -99999;
      Color2[0]= -99999;
      px[0]= -99999;
      py[0]= -99999;
      pz[0]= -99999;
      E[0]= -99999;
      m[0]= -99999;
      lifetime[0]= -99999;
      spin[0]= -99999;
     for (int i=1 ; i<npart+1 ; i++) { //start at one
        ff >> Id[i] >> Status[i] >> Mother1[i] >> Mother2[i] >> Color1[i] >> Color2[i]
           >> px[i] >> py[i] >> pz[i] >> E[i] >> m[i] >> lifetime[i] >> spin[i] ;
        buffer << Id[i] << " " << Status[i] << " " << std::endl;
        line++;
        v[i] = new TLorentzVector(px[i], py[i], pz[i], E[i]);
        if (Status[i]==-1) continue; // status -1 = initial quark ==> skip
        if (Id[i]==6)  top=i;
        if (Id[i]==-6) topbar=i;
        if (Id[i]>6000000) zprime=i;


//------------------------*Light Quarks and B Quarks*-----------------------------

	//b qaurks:
	if (abs(Id[Mother1[i]])==6)
	{ //Mother top
	
		if (abs(Id[i])==5)
		{// bquarks
		v_bq.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		hPt_b->Fill(v_bq.Pt());
		hEta_b->Fill(v_bq.Eta());
		hPhi_b->Fill(v_bq.Phi());
		
		     if (bq[0]==-1) bq[0]=i;
		else if (bq[1]==-1) bq[1]=i;
		else if (bq[2]==-1) bq[2]=i;
		else if (bq[3]==-1) bq[3]=i;
		else cout << "ERROR : more than 4 b quarks" << endl;
		n_bq++;
		}
	}
	
	//Light Quarks:
	if ( abs(Id[Mother1[i]])==24)
	{// Mother W boson
	
		if (abs(Id[i]) == 1)
		{// up quarks
		v_upq.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		hPt_up->Fill(v_upq.Pt());
		hEta_up->Fill(v_upq.Eta());
		hPhi_up->Fill(v_upq.Phi());
		
		     if (lq1[0] == -1 ) lq1[0]=i;
		else if (lq1[1] == -1 ) lq1[1]=i;
		else if (lq1[2] == -1 ) lq1[2]=i;
		else if (lq1[3] == -1 ) lq1[3]=i;
		else if (lq1[4] == -1 ) lq1[4]=i;
		else cout<<"more than 4 up quarks"<<endl;
		n_upq++;
		}
		
		if (abs(Id[i]) == 2)
		{// down quarks
		v_downq.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		hPt_down->Fill(v_downq.Pt());
		hEta_down->Fill(v_downq.Eta());
		hPhi_down->Fill(v_downq.Phi());
		
		     if (lq2[0] == -1 ) lq2[0]=i;
		else if (lq2[1] == -1 ) lq2[1]=i;
		else if (lq2[2] == -1 ) lq2[2]=i;
		else if (lq2[3] == -1 ) lq2[3]=i;
		else if (lq2[4] == -1 ) lq2[4]=i;
		else cout<<"more than 4 down quarks"<<endl;
		n_downq++;
		}
		
		if (abs(Id[i]) == 3)
		{// charm quarks
		v_charmq.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		hPt_charm->Fill(v_charmq.Pt());
		hEta_charm->Fill(v_charmq.Eta());
		hPhi_charm->Fill(v_charmq.Phi());
		
		     if (lq3[0] == -1 ) lq3[0]=i;
		else if (lq3[1] == -1 ) lq3[1]=i;
		else if (lq3[2] == -1 ) lq3[2]=i;
		else if (lq3[3] == -1 ) lq3[3]=i;
		else if (lq3[4] == -1 ) lq3[4]=i;
		else cout<<"more than 4 charm quarks"<<endl;
		n_charmq++;
		}
		
		if (abs(Id[i]) == 4)
		{// strange quarks
		v_strangeq.SetPxPyPzE(v[i]->Px(),v[i]->Py(),v[i]->Pz(),v[i]->E());
		hPt_strange->Fill(v_strangeq.Pt());
		hEta_strange->Fill(v_strangeq.Eta());
		hPhi_strange->Fill(v_strangeq.Phi());
		
		     if (lq4[0] == -1 ) lq4[0]=i;
		else if (lq4[1] == -1 ) lq4[1]=i;
		else if (lq4[2] == -1 ) lq4[2]=i;
		else if (lq4[3] == -1 ) lq4[3]=i;
		else if (lq4[4] == -1 ) lq4[4]=i;
		else cout<<"more than 4 strange quarks"<<endl;
		n_strangeq++;
		}
	}	  	  
} //loop of i

cout<<"no. b quarks: "<<n_bq<<endl;
cout<<"no. up quarks: "<<n_upq<<endl;
cout<<"no. down quarks: "<<n_downq<<endl;
cout<<"no. charm quarks: "<<n_charmq<<endl;
cout<<"no. strange quarks: "<<n_strangeq<<endl;

//--------------------------------*Filling Histograms*----------------------------------//

double bq_pt[4] = {v[bq[0]]->Pt(),v[bq[1]]->Pt(),v[bq[2]]->Pt(),v[bq[3]]->Pt()};	
        std::sort(bq_pt,bq_pt+4);
double bq_eta[4] = {v[bq[0]]->Eta(),v[bq[1]]->Eta(),v[bq[2]]->Eta(),v[bq[3]]->Eta()};	
        std::sort(bq_eta,bq_eta+4);
double bq_phi[4] = {v[bq[0]]->Phi(),v[bq[1]]->Phi(),v[bq[2]]->Phi(),v[bq[3]]->Phi()};	
        std::sort(bq_phi,bq_phi+4);
        
      v_bq1.SetPxPyPzE(v[bq[0]]->Px(),v[bq[0]]->Py(),v[bq[0]]->Pz(),v[bq[0]]->E());
      v_bq2.SetPxPyPzE(v[bq[1]]->Px(),v[bq[1]]->Py(),v[bq[1]]->Pz(),v[bq[1]]->E());
      v_bq3.SetPxPyPzE(v[bq[2]]->Px(),v[bq[2]]->Py(),v[bq[2]]->Pz(),v[bq[2]]->E());
      v_bq4.SetPxPyPzE(v[bq[3]]->Px(),v[bq[3]]->Py(),v[bq[3]]->Pz(),v[bq[3]]->E());
      
      //b quarks:
      for (int nbq=0; nbq<4 ; nbq++)
      {
	  hPt_bq[nbq]->Fill( bq_pt[nbq] );
          hEta_bq[nbq]->Fill( bq_eta[nbq] );
          hPhi_bq[nbq]->Fill( bq_phi[nbq] );
      }    

if (n_upq == 5 || n_downq == 5)
{
double up_quark_pt[5] = {v[lq1[0]]->Pt(),v[lq1[1]]->Pt(),v[lq1[2]]->Pt(),v[lq1[3]]->Pt(),v[lq1[4]]->Pt()};	
        std::sort(up_quark_pt, up_quark_pt+5);
double up_quark_eta[5] = {v[lq1[0]]->Eta(),v[lq1[1]]->Eta(),v[lq1[2]]->Eta(),v[lq1[3]]->Eta(),v[lq1[4]]->Eta()};	
        std::sort(up_quark_eta, up_quark_eta+5);
double up_quark_phi[5] = {v[lq1[0]]->Phi(),v[lq1[1]]->Phi(),v[lq1[2]]->Phi(),v[lq1[3]]->Phi(),v[lq1[4]]->Phi()};	
        std::sort(up_quark_phi, up_quark_phi+5);

double down_quark_pt[5] = {v[lq2[0]]->Pt(),v[lq2[1]]->Pt(),v[lq2[2]]->Pt(),v[lq2[3]]->Pt(),v[lq2[4]]->Pt()};	
        std::sort(down_quark_pt, down_quark_pt+5);
double down_quark_eta[5] = {v[lq2[0]]->Eta(),v[lq2[1]]->Eta(),v[lq2[2]]->Eta(),v[lq2[3]]->Eta(),v[lq2[4]]->Eta()};	
        std::sort(down_quark_eta, down_quark_eta+5);
double down_quark_phi[5] = {v[lq2[0]]->Phi(),v[lq2[1]]->Phi(),v[lq2[2]]->Phi(),v[lq2[3]]->Phi(),v[lq2[4]]->Phi()};	
        std::sort(down_quark_phi, down_quark_phi+5);
                      
      v_up1.SetPxPyPzE(v[lq1[0]]->Px(),v[lq1[0]]->Py(),v[lq1[0]]->Pz(),v[lq1[0]]->E());
      v_up2.SetPxPyPzE(v[lq1[1]]->Px(),v[lq1[1]]->Py(),v[lq1[1]]->Pz(),v[lq1[1]]->E());
      v_up3.SetPxPyPzE(v[lq1[2]]->Px(),v[lq1[2]]->Py(),v[lq1[2]]->Pz(),v[lq1[2]]->E());
      v_up4.SetPxPyPzE(v[lq1[3]]->Px(),v[lq1[3]]->Py(),v[lq1[3]]->Pz(),v[lq1[3]]->E());
      v_up5.SetPxPyPzE(v[lq1[4]]->Px(),v[lq1[4]]->Py(),v[lq1[4]]->Pz(),v[lq1[4]]->E());
      for (int n_up=0; n_up<5 ; n_up++)
      {
	  hPt_lq1[n_up]->Fill( up_quark_pt[n_up] );
          hEta_lq1[n_up]->Fill( up_quark_eta[n_up] );
          hPhi_lq1[n_up]->Fill( up_quark_phi[n_up] );
      }
      cout<<"up_quark Pt :  "<<up_quark_pt[0]<<",  "<<up_quark_pt[1]<<",  "<<up_quark_pt[2]<<",  "<<up_quark_pt[3]<<",  "<<up_quark_pt[4]<<endl;
      
           
      v_down1.SetPxPyPzE(v[lq2[0]]->Px(),v[lq2[0]]->Py(),v[lq2[0]]->Pz(),v[lq2[0]]->E());
      v_down2.SetPxPyPzE(v[lq2[1]]->Px(),v[lq2[1]]->Py(),v[lq2[1]]->Pz(),v[lq2[1]]->E());
      v_down3.SetPxPyPzE(v[lq2[2]]->Px(),v[lq2[2]]->Py(),v[lq2[2]]->Pz(),v[lq2[2]]->E());
      v_down4.SetPxPyPzE(v[lq2[3]]->Px(),v[lq2[3]]->Py(),v[lq2[3]]->Pz(),v[lq2[3]]->E());
      v_down5.SetPxPyPzE(v[lq2[4]]->Px(),v[lq2[4]]->Py(),v[lq2[4]]->Pz(),v[lq2[4]]->E());
      for (int n_down=0; n_down<5 ; n_down++)
      {
      hPt_lq2[n_down]->Fill( down_quark_pt[n_down] );
      hEta_lq2[n_down]->Fill( down_quark_eta[n_down] );
      hPhi_lq2[n_down]->Fill( down_quark_phi[n_down] );
      }     
      cout<<"down_quark Pt :  "<<down_quark_pt[0]<<",  "<<down_quark_pt[1]<<",  "<<down_quark_pt[2]<<",  "<<down_quark_pt[3]<<",  "<<down_quark_pt[4]<<endl;
             
        // w boson kinematics:
        v_w_boson1 = v_up1+v_down1;
        v_w_boson2 = v_up2+v_down2;
        v_w_boson3 = v_up3+v_down3;
        v_w_boson4 = v_up4+v_down4;
        v_w_boson5 = v_up5+v_down5;
        
TLorentzVector v_wbosons[5]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4,  v_w_boson5};
	for (int wb=0; wb<5 ; wb++)
	{
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
        hM_wboson[wb]->Fill( v_wbosons[wb].M() );
        }
       
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
}

if (n_upq==4 || n_downq == 4)
{
double up_quark_pt[5] = {v[lq1[0]]->Pt(),v[lq1[1]]->Pt(),v[lq1[2]]->Pt(),v[lq1[3]]->Pt()};	
        std::sort(up_quark_pt, up_quark_pt+4);
double up_quark_eta[5] = {v[lq1[0]]->Eta(),v[lq1[1]]->Eta(),v[lq1[2]]->Eta(),v[lq1[3]]->Eta()};	
        std::sort(up_quark_eta, up_quark_eta+4);
double up_quark_phi[5] = {v[lq1[0]]->Phi(),v[lq1[1]]->Phi(),v[lq1[2]]->Phi(),v[lq1[3]]->Phi()};	
        std::sort(up_quark_phi, up_quark_phi+4);
        
double down_quark_pt[5] = {v[lq2[0]]->Pt(),v[lq2[1]]->Pt(),v[lq2[2]]->Pt(),v[lq2[3]]->Pt()};	
        std::sort(down_quark_pt, down_quark_pt+4);
double down_quark_eta[5] = {v[lq2[0]]->Eta(),v[lq2[1]]->Eta(),v[lq2[2]]->Eta(),v[lq2[3]]->Eta()};	
        std::sort(down_quark_eta, down_quark_eta+4);
double down_quark_phi[5] = {v[lq2[0]]->Phi(),v[lq2[1]]->Phi(),v[lq2[2]]->Phi(),v[lq2[3]]->Phi()};	
        std::sort(down_quark_phi, down_quark_phi+4);

double charm_quark_pt[5] = {v[lq3[0]]->Pt()};	
        std::sort(charm_quark_pt, charm_quark_pt+1);
double charm_quark_eta[5] = {v[lq3[0]]->Eta()};	
        std::sort(charm_quark_eta, charm_quark_eta+1);
double charm_quark_phi[5] = {v[lq3[0]]->Phi()};	
        std::sort(charm_quark_phi, charm_quark_phi+1);
        
double strange_quark_pt[5] = {v[lq4[0]]->Pt()};	
        std::sort(strange_quark_pt, strange_quark_pt+1);
double strange_quark_eta[5] = {v[lq4[0]]->Eta()};	
        std::sort(strange_quark_eta, strange_quark_eta+1);
double strange_quark_phi[5] = {v[lq4[0]]->Phi()};	
        std::sort(strange_quark_phi, strange_quark_phi+1);
                
      v_up1.SetPxPyPzE(v[lq1[0]]->Px(),v[lq1[0]]->Py(),v[lq1[0]]->Pz(),v[lq1[0]]->E());
      v_up2.SetPxPyPzE(v[lq1[1]]->Px(),v[lq1[1]]->Py(),v[lq1[1]]->Pz(),v[lq1[1]]->E());
      v_up3.SetPxPyPzE(v[lq1[2]]->Px(),v[lq1[2]]->Py(),v[lq1[2]]->Pz(),v[lq1[2]]->E());
      v_up4.SetPxPyPzE(v[lq1[3]]->Px(),v[lq1[3]]->Py(),v[lq1[3]]->Pz(),v[lq1[3]]->E());
      for (int n_up=0; n_up<4 ; n_up++)
      {
	  hPt_lq1[n_up]->Fill( up_quark_pt[n_up] );
          hEta_lq1[n_up]->Fill( up_quark_eta[n_up] );
          hPhi_lq1[n_up]->Fill( up_quark_phi[n_up] );
      }
      
      v_down1.SetPxPyPzE(v[lq2[0]]->Px(),v[lq2[0]]->Py(),v[lq2[0]]->Pz(),v[lq2[0]]->E());
      v_down2.SetPxPyPzE(v[lq2[1]]->Px(),v[lq2[1]]->Py(),v[lq2[1]]->Pz(),v[lq2[1]]->E());
      v_down3.SetPxPyPzE(v[lq2[2]]->Px(),v[lq2[2]]->Py(),v[lq2[2]]->Pz(),v[lq2[2]]->E());
      v_down4.SetPxPyPzE(v[lq2[3]]->Px(),v[lq2[3]]->Py(),v[lq2[3]]->Pz(),v[lq2[3]]->E());      
      for (int n_down=0; n_down<4 ; n_down++)
      {
      hPt_lq2[n_down]->Fill( down_quark_pt[n_down] );
      hEta_lq2[n_down]->Fill( down_quark_eta[n_down] );
      hPhi_lq2[n_down]->Fill( down_quark_phi[n_down] );
      }    
           
      v_charm1.SetPxPyPzE(v[lq3[0]]->Px(),v[lq3[0]]->Py(),v[lq3[0]]->Pz(),v[lq3[0]]->E());
      for (int n_charm=0; n_charm<1 ; n_charm++)
      {
	  hPt_lq3[n_charm]->Fill( charm_quark_pt[n_charm] );
          hEta_lq3[n_charm]->Fill( charm_quark_eta[n_charm] );
          hPhi_lq3[n_charm]->Fill( charm_quark_phi[n_charm] );
      }
      
      v_strange1.SetPxPyPzE(v[lq4[0]]->Px(),v[lq4[0]]->Py(),v[lq4[0]]->Pz(),v[lq4[0]]->E());
      for (int n_strange=0; n_strange<1 ; n_strange++)
      {
	  hPt_lq4[n_strange]->Fill( strange_quark_pt[n_strange] );
          hEta_lq4[n_strange]->Fill( strange_quark_eta[n_strange] );
          hPhi_lq4[n_strange]->Fill( strange_quark_phi[n_strange] );
      }
      
      cout<<"up_quark Pt :  "<<up_quark_pt[0]<<",  "<<up_quark_pt[1]<<",  "<<up_quark_pt[2]<<",  "<<up_quark_pt[3]<<",  "<<endl;
      cout<<"down_quark Pt :  "<<down_quark_pt[0]<<",  "<<down_quark_pt[1]<<",  "<<down_quark_pt[2]<<",  "<<down_quark_pt[3]<<",  "<<endl; 
      cout<<"charm_quark Pt :  "<<charm_quark_pt[0]<<endl;
      cout<<"strange_quark Pt :  "<<strange_quark_pt[0]<<endl;
      
        // w boson kinematics:
        v_w_boson1 = v_up1+v_down1;
        v_w_boson2 = v_up2+v_down2;
        v_w_boson3 = v_up3+v_down3;
        v_w_boson4 = v_up4+v_down4;
        v_w_boson5 = v_charm1+v_strange1;
        
TLorentzVector v_wbosons[5]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4,  v_w_boson5};
	for (int wb=0; wb<5 ; wb++)
	{
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
        hM_wboson[wb]->Fill( v_wbosons[wb].M() );
        }
       
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
}

if (n_upq==3 || n_downq == 3)
{
double up_quark_pt[5] = {v[lq1[0]]->Pt(),v[lq1[1]]->Pt(),v[lq1[2]]->Pt()};	
        std::sort(up_quark_pt, up_quark_pt+3);
double up_quark_eta[5] = {v[lq1[0]]->Eta(),v[lq1[1]]->Eta(),v[lq1[2]]->Eta()};	
        std::sort(up_quark_eta, up_quark_eta+3);
double up_quark_phi[5] = {v[lq1[0]]->Phi(),v[lq1[1]]->Phi(),v[lq1[2]]->Phi()};	
        std::sort(up_quark_phi, up_quark_phi+3);
        
double down_quark_pt[5] = {v[lq2[0]]->Pt(),v[lq2[1]]->Pt(),v[lq2[2]]->Pt()};	
        std::sort(down_quark_pt, down_quark_pt+3);
double down_quark_eta[5] = {v[lq2[0]]->Eta(),v[lq2[1]]->Eta(),v[lq2[2]]->Eta()};	
        std::sort(down_quark_eta, down_quark_eta+3);
double down_quark_phi[5] = {v[lq2[0]]->Phi(),v[lq2[1]]->Phi(),v[lq2[2]]->Phi()};	
        std::sort(down_quark_phi, down_quark_phi+3);

double charm_quark_pt[5] = {v[lq3[0]]->Pt(),v[lq3[1]]->Pt()};	
        std::sort(charm_quark_pt, charm_quark_pt+2);
double charm_quark_eta[5] = {v[lq3[0]]->Eta(),v[lq3[1]]->Eta()};	
        std::sort(charm_quark_eta, charm_quark_eta+2);
double charm_quark_phi[5] = {v[lq3[0]]->Phi(),v[lq3[1]]->Phi()};	
        std::sort(charm_quark_phi, charm_quark_phi+2);
        
double strange_quark_pt[5] = {v[lq4[0]]->Pt(),v[lq4[1]]->Pt()};	
        std::sort(strange_quark_pt, strange_quark_pt+2);
double strange_quark_eta[5] = {v[lq4[0]]->Eta(),v[lq4[1]]->Eta()};	
        std::sort(strange_quark_eta, strange_quark_eta+2);
double strange_quark_phi[5] = {v[lq4[0]]->Phi(),v[lq4[1]]->Phi()};	
        std::sort(strange_quark_phi, strange_quark_phi+2);
                
      v_up1.SetPxPyPzE(v[lq1[0]]->Px(),v[lq1[0]]->Py(),v[lq1[0]]->Pz(),v[lq1[0]]->E());
      v_up2.SetPxPyPzE(v[lq1[1]]->Px(),v[lq1[1]]->Py(),v[lq1[1]]->Pz(),v[lq1[1]]->E());
      v_up3.SetPxPyPzE(v[lq1[2]]->Px(),v[lq1[2]]->Py(),v[lq1[2]]->Pz(),v[lq1[2]]->E());
      for (int n_up=0; n_up<3 ; n_up++)
      {
	  hPt_lq1[n_up]->Fill( up_quark_pt[n_up] );
          hEta_lq1[n_up]->Fill( up_quark_eta[n_up] );
          hPhi_lq1[n_up]->Fill( up_quark_phi[n_up] );
      }
      
      v_down1.SetPxPyPzE(v[lq2[0]]->Px(),v[lq2[0]]->Py(),v[lq2[0]]->Pz(),v[lq2[0]]->E());
      v_down2.SetPxPyPzE(v[lq2[1]]->Px(),v[lq2[1]]->Py(),v[lq2[1]]->Pz(),v[lq2[1]]->E());
      v_down3.SetPxPyPzE(v[lq2[2]]->Px(),v[lq2[2]]->Py(),v[lq2[2]]->Pz(),v[lq2[2]]->E());    
      for (int n_down=0; n_down<3 ; n_down++)
      {
      hPt_lq2[n_down]->Fill( down_quark_pt[n_down] );
      hEta_lq2[n_down]->Fill( down_quark_eta[n_down] );
      hPhi_lq2[n_down]->Fill( down_quark_phi[n_down] );
      }    
           
      v_charm1.SetPxPyPzE(v[lq3[0]]->Px(),v[lq3[0]]->Py(),v[lq3[0]]->Pz(),v[lq3[0]]->E());
      v_charm2.SetPxPyPzE(v[lq3[1]]->Px(),v[lq3[1]]->Py(),v[lq3[1]]->Pz(),v[lq3[1]]->E());
      for (int n_charm=0; n_charm<2 ; n_charm++)
      {
	  hPt_lq3[n_charm]->Fill( charm_quark_pt[n_charm] );
          hEta_lq3[n_charm]->Fill( charm_quark_eta[n_charm] );
          hPhi_lq3[n_charm]->Fill( charm_quark_phi[n_charm] );
      }
      
      v_strange1.SetPxPyPzE(v[lq4[0]]->Px(),v[lq4[0]]->Py(),v[lq4[0]]->Pz(),v[lq4[0]]->E());
      v_strange2.SetPxPyPzE(v[lq4[1]]->Px(),v[lq4[1]]->Py(),v[lq4[1]]->Pz(),v[lq4[1]]->E());
      for (int n_strange=0; n_strange<2 ; n_strange++)
      {
	  hPt_lq4[n_strange]->Fill( strange_quark_pt[n_strange] );
          hEta_lq4[n_strange]->Fill( strange_quark_eta[n_strange] );
          hPhi_lq4[n_strange]->Fill( strange_quark_phi[n_strange] );
      }
      
      cout<<"up_quark Pt :  "<<up_quark_pt[0]<<",  "<<up_quark_pt[1]<<",  "<<up_quark_pt[2]<<",  "<<endl;
      cout<<"down_quark Pt :  "<<down_quark_pt[0]<<",  "<<down_quark_pt[1]<<",  "<<down_quark_pt[2]<<",  "<<endl; 
      cout<<"charm_quark Pt :  "<<charm_quark_pt[0]<<", "<<charm_quark_pt[1]<<endl;
      cout<<"strange_quark Pt :  "<<strange_quark_pt[0]<<", "<<strange_quark_pt[1]<<endl;
      
        // w boson kinematics:
        v_w_boson1 = v_up1+v_down1;
        v_w_boson2 = v_up2+v_down2;
        v_w_boson3 = v_up3+v_down3;
        v_w_boson4 = v_charm1+v_strange1;
        v_w_boson5 = v_charm2+v_strange2;
        
TLorentzVector v_wbosons[5]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4,  v_w_boson5};
	for (int wb=0; wb<5 ; wb++)
	{
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
        hM_wboson[wb]->Fill( v_wbosons[wb].M() );
        }
       
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
}

if (n_upq==2 || n_downq == 2)
{
double up_quark_pt[5] = {v[lq1[0]]->Pt(),v[lq1[1]]->Pt()};	
        std::sort(up_quark_pt, up_quark_pt+2);
double up_quark_eta[5] = {v[lq1[0]]->Eta(),v[lq1[1]]->Eta()};	
        std::sort(up_quark_eta, up_quark_eta+2);
double up_quark_phi[5] = {v[lq1[0]]->Phi(),v[lq1[1]]->Phi()};	
        std::sort(up_quark_phi, up_quark_phi+2);
        
double down_quark_pt[5] = {v[lq2[0]]->Pt(),v[lq2[1]]->Pt()};	
        std::sort(down_quark_pt, down_quark_pt+2);
double down_quark_eta[5] = {v[lq2[0]]->Eta(),v[lq2[1]]->Eta()};	
        std::sort(down_quark_eta, down_quark_eta+2);
double down_quark_phi[5] = {v[lq2[0]]->Phi(),v[lq2[1]]->Phi()};	
        std::sort(down_quark_phi, down_quark_phi+2);

double charm_quark_pt[5] = {v[lq3[0]]->Pt(),v[lq3[1]]->Pt(),v[lq3[2]]->Pt()};	
        std::sort(charm_quark_pt, charm_quark_pt+3);
double charm_quark_eta[5] = {v[lq3[0]]->Eta(),v[lq3[1]]->Eta(),v[lq3[2]]->Eta()};	
        std::sort(charm_quark_eta, charm_quark_eta+3);
double charm_quark_phi[5] = {v[lq3[0]]->Phi(),v[lq3[1]]->Phi(),v[lq3[2]]->Phi()};	
        std::sort(charm_quark_phi, charm_quark_phi+3);
        
double strange_quark_pt[5] = {v[lq4[0]]->Pt(),v[lq4[1]]->Pt(),v[lq4[2]]->Pt()};	
        std::sort(strange_quark_pt, strange_quark_pt+3);
double strange_quark_eta[5] = {v[lq4[0]]->Eta(),v[lq4[1]]->Eta(),v[lq4[2]]->Eta()};	
        std::sort(strange_quark_eta, strange_quark_eta+3);
double strange_quark_phi[5] = {v[lq4[0]]->Phi(),v[lq4[1]]->Phi(),v[lq4[2]]->Phi()};	
        std::sort(strange_quark_phi, strange_quark_phi+3);
                
      v_up1.SetPxPyPzE(v[lq1[0]]->Px(),v[lq1[0]]->Py(),v[lq1[0]]->Pz(),v[lq1[0]]->E());
      v_up2.SetPxPyPzE(v[lq1[1]]->Px(),v[lq1[1]]->Py(),v[lq1[1]]->Pz(),v[lq1[1]]->E());
      for (int n_up=0; n_up<2 ; n_up++)
      {
	  hPt_lq1[n_up]->Fill( up_quark_pt[n_up] );
          hEta_lq1[n_up]->Fill( up_quark_eta[n_up] );
          hPhi_lq1[n_up]->Fill( up_quark_phi[n_up] );
      }
      
      v_down1.SetPxPyPzE(v[lq2[0]]->Px(),v[lq2[0]]->Py(),v[lq2[0]]->Pz(),v[lq2[0]]->E());
      v_down2.SetPxPyPzE(v[lq2[1]]->Px(),v[lq2[1]]->Py(),v[lq2[1]]->Pz(),v[lq2[1]]->E());   
      for (int n_down=0; n_down<2 ; n_down++)
      {
      hPt_lq2[n_down]->Fill( down_quark_pt[n_down] );
      hEta_lq2[n_down]->Fill( down_quark_eta[n_down] );
      hPhi_lq2[n_down]->Fill( down_quark_phi[n_down] );
      }    
           
      v_charm1.SetPxPyPzE(v[lq3[0]]->Px(),v[lq3[0]]->Py(),v[lq3[0]]->Pz(),v[lq3[0]]->E());
      v_charm2.SetPxPyPzE(v[lq3[1]]->Px(),v[lq3[1]]->Py(),v[lq3[1]]->Pz(),v[lq3[1]]->E());
      v_charm3.SetPxPyPzE(v[lq3[2]]->Px(),v[lq3[2]]->Py(),v[lq3[2]]->Pz(),v[lq3[2]]->E());
      for (int n_charm=0; n_charm<3 ; n_charm++)
      {
	  hPt_lq3[n_charm]->Fill( charm_quark_pt[n_charm] );
          hEta_lq3[n_charm]->Fill( charm_quark_eta[n_charm] );
          hPhi_lq3[n_charm]->Fill( charm_quark_phi[n_charm] );
      }
      
      v_strange1.SetPxPyPzE(v[lq4[0]]->Px(),v[lq4[0]]->Py(),v[lq4[0]]->Pz(),v[lq4[0]]->E());
      v_strange2.SetPxPyPzE(v[lq4[1]]->Px(),v[lq4[1]]->Py(),v[lq4[1]]->Pz(),v[lq4[1]]->E());
      v_strange3.SetPxPyPzE(v[lq4[2]]->Px(),v[lq4[2]]->Py(),v[lq4[2]]->Pz(),v[lq4[2]]->E());
      
      for (int n_strange=0; n_strange<3 ; n_strange++)
      {
	  hPt_lq4[n_strange]->Fill( strange_quark_pt[n_strange] );
          hEta_lq4[n_strange]->Fill( strange_quark_eta[n_strange] );
          hPhi_lq4[n_strange]->Fill( strange_quark_phi[n_strange] );
      }
      
      cout<<"up_quark Pt :  "<<up_quark_pt[0]<<",  "<<up_quark_pt[1]<<",  "<<endl;
      cout<<"down_quark Pt :  "<<down_quark_pt[0]<<",  "<<down_quark_pt[1]<<",  "<<endl; 
      cout<<"charm_quark Pt :  "<<charm_quark_pt[0]<<", "<<charm_quark_pt[1]<<", "<<charm_quark_pt[2]<<endl;
      cout<<"strange_quark Pt :  "<<strange_quark_pt[0]<<", "<<strange_quark_pt[1]<<", "<<strange_quark_pt[2]<<endl;
      
        // w boson kinematics:
        v_w_boson1 = v_up1+v_down1;
        v_w_boson2 = v_up2+v_down2;
        v_w_boson3 = v_charm1+v_strange1;
        v_w_boson4 = v_charm2+v_strange2;
        v_w_boson5 = v_charm3+v_strange3;
        
TLorentzVector v_wbosons[5]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4,  v_w_boson5};
	for (int wb=0; wb<5 ; wb++)
	{
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
        hM_wboson[wb]->Fill( v_wbosons[wb].M() );
        }
       
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
}

if (n_upq==1 || n_downq == 1)
{
double up_quark_pt[5] = {v[lq1[0]]->Pt()};	
        std::sort(up_quark_pt, up_quark_pt+1);
double up_quark_eta[5] = {v[lq1[0]]->Eta()};	
        std::sort(up_quark_eta, up_quark_eta+1);
double up_quark_phi[5] = {v[lq1[0]]->Phi()};	
        std::sort(up_quark_phi, up_quark_phi+1);
        
double down_quark_pt[5] = {v[lq2[0]]->Pt()};	
        std::sort(down_quark_pt, down_quark_pt+1);
double down_quark_eta[5] = {v[lq2[0]]->Eta()};	
        std::sort(down_quark_eta, down_quark_eta+1);
double down_quark_phi[5] = {v[lq2[0]]->Phi()};	
        std::sort(down_quark_phi, down_quark_phi+1);

double charm_quark_pt[5] = {v[lq3[0]]->Pt(),v[lq3[1]]->Pt(),v[lq3[2]]->Pt(),v[lq3[3]]->Pt()};	
        std::sort(charm_quark_pt, charm_quark_pt+4);
double charm_quark_eta[5] = {v[lq3[0]]->Eta(),v[lq3[1]]->Eta(),v[lq3[2]]->Eta(),v[lq3[3]]->Eta()};	
        std::sort(charm_quark_eta, charm_quark_eta+4);
double charm_quark_phi[5] = {v[lq3[0]]->Phi(),v[lq3[1]]->Phi(),v[lq3[2]]->Phi(),v[lq3[3]]->Phi()};	
        std::sort(charm_quark_phi, charm_quark_phi+4);
        
double strange_quark_pt[5] = {v[lq4[0]]->Pt(),v[lq4[1]]->Pt(),v[lq4[2]]->Pt(),v[lq4[3]]->Pt()};	
        std::sort(strange_quark_pt, strange_quark_pt+4);
double strange_quark_eta[5] = {v[lq4[0]]->Eta(),v[lq4[1]]->Eta(),v[lq4[2]]->Eta(),v[lq4[3]]->Eta()};	
        std::sort(strange_quark_eta, strange_quark_eta+4);
double strange_quark_phi[5] = {v[lq4[0]]->Phi(),v[lq4[1]]->Phi(),v[lq4[2]]->Phi(),v[lq4[3]]->Phi()};	
        std::sort(strange_quark_phi, strange_quark_phi+4);
                
      v_up1.SetPxPyPzE(v[lq1[0]]->Px(),v[lq1[0]]->Py(),v[lq1[0]]->Pz(),v[lq1[0]]->E());
      for (int n_up=0; n_up<1 ; n_up++)
      {
	  hPt_lq1[n_up]->Fill( up_quark_pt[n_up] );
          hEta_lq1[n_up]->Fill( up_quark_eta[n_up] );
          hPhi_lq1[n_up]->Fill( up_quark_phi[n_up] );
      }
      
      v_down1.SetPxPyPzE(v[lq2[0]]->Px(),v[lq2[0]]->Py(),v[lq2[0]]->Pz(),v[lq2[0]]->E());  
      for (int n_down=0; n_down<1 ; n_down++)
      {
      hPt_lq2[n_down]->Fill( down_quark_pt[n_down] );
      hEta_lq2[n_down]->Fill( down_quark_eta[n_down] );
      hPhi_lq2[n_down]->Fill( down_quark_phi[n_down] );
      }    
           
      v_charm1.SetPxPyPzE(v[lq3[0]]->Px(),v[lq3[0]]->Py(),v[lq3[0]]->Pz(),v[lq3[0]]->E());
      v_charm2.SetPxPyPzE(v[lq3[1]]->Px(),v[lq3[1]]->Py(),v[lq3[1]]->Pz(),v[lq3[1]]->E());
      v_charm3.SetPxPyPzE(v[lq3[2]]->Px(),v[lq3[2]]->Py(),v[lq3[2]]->Pz(),v[lq3[2]]->E());
      v_charm4.SetPxPyPzE(v[lq3[3]]->Px(),v[lq3[3]]->Py(),v[lq3[3]]->Pz(),v[lq3[3]]->E());
      for (int n_charm=0; n_charm<4 ; n_charm++)
      {
	  hPt_lq3[n_charm]->Fill( charm_quark_pt[n_charm] );
          hEta_lq3[n_charm]->Fill( charm_quark_eta[n_charm] );
          hPhi_lq3[n_charm]->Fill( charm_quark_phi[n_charm] );
      }
      
      v_strange1.SetPxPyPzE(v[lq4[0]]->Px(),v[lq4[0]]->Py(),v[lq4[0]]->Pz(),v[lq4[0]]->E());
      v_strange2.SetPxPyPzE(v[lq4[1]]->Px(),v[lq4[1]]->Py(),v[lq4[1]]->Pz(),v[lq4[1]]->E());
      v_strange3.SetPxPyPzE(v[lq4[2]]->Px(),v[lq4[2]]->Py(),v[lq4[2]]->Pz(),v[lq4[2]]->E());
      v_strange4.SetPxPyPzE(v[lq4[3]]->Px(),v[lq4[3]]->Py(),v[lq4[3]]->Pz(),v[lq4[3]]->E());
      for (int n_strange=0; n_strange<4 ; n_strange++)
      {
	  hPt_lq4[n_strange]->Fill( strange_quark_pt[n_strange] );
          hEta_lq4[n_strange]->Fill( strange_quark_eta[n_strange] );
          hPhi_lq4[n_strange]->Fill( strange_quark_phi[n_strange] );
      }
      
      cout<<"up_quark Pt :  "<<up_quark_pt[0]<<endl;
      cout<<"down_quark Pt :  "<<down_quark_pt[0]<<endl; 
      cout<<"charm_quark Pt :  "<<charm_quark_pt[0]<<", "<<charm_quark_pt[1]<<", "<<charm_quark_pt[2]<<", "<<charm_quark_pt[3]<<endl;
      cout<<"strange_quark Pt :  "<<strange_quark_pt[0]<<", "<<strange_quark_pt[1]<<", "<<strange_quark_pt[2]<<", "<<strange_quark_pt[3]<<endl;
      
        // w boson kinematics:
        v_w_boson1 = v_up1+v_down1;
        v_w_boson2 = v_charm1+v_strange1;
        v_w_boson3 = v_charm2+v_strange2;
        v_w_boson4 = v_charm3+v_strange3;
        v_w_boson5 = v_charm4+v_strange4;
        
TLorentzVector v_wbosons[5]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4,  v_w_boson5};
	for (int wb=0; wb<5 ; wb++)
	{
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
        hM_wboson[wb]->Fill( v_wbosons[wb].M() );
        }
       
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
}

if (n_charmq == 5)
{
double charm_quark_pt[5] = {v[lq3[0]]->Pt(),v[lq3[1]]->Pt(),v[lq3[2]]->Pt(),v[lq3[3]]->Pt(),v[lq3[4]]->Pt()};	
        std::sort(charm_quark_pt, charm_quark_pt+5);
double charm_quark_eta[5] = {v[lq3[0]]->Eta(),v[lq3[1]]->Eta(),v[lq3[2]]->Eta(),v[lq3[3]]->Eta(),v[lq3[4]]->Eta()};	
        std::sort(charm_quark_eta, charm_quark_eta+5);
double charm_quark_phi[5] = {v[lq3[0]]->Phi(),v[lq3[1]]->Phi(),v[lq3[2]]->Phi(),v[lq3[3]]->Phi(),v[lq3[4]]->Phi()};	
        std::sort(charm_quark_phi, charm_quark_phi+5);

double strange_quark_pt[5] = {v[lq4[0]]->Pt(),v[lq4[1]]->Pt(),v[lq4[2]]->Pt(),v[lq4[3]]->Pt(),v[lq4[4]]->Pt()};	
        std::sort(strange_quark_pt, strange_quark_pt+5);
double strange_quark_eta[5] = {v[lq4[0]]->Eta(),v[lq4[1]]->Eta(),v[lq4[2]]->Eta(),v[lq4[3]]->Eta(),v[lq4[4]]->Eta()};	
        std::sort(strange_quark_eta, strange_quark_eta+5);
double strange_quark_phi[5] = {v[lq4[0]]->Phi(),v[lq4[1]]->Phi(),v[lq4[2]]->Phi(),v[lq4[3]]->Phi(),v[lq4[4]]->Phi()};	
        std::sort(strange_quark_phi, strange_quark_phi+5);
                      
      v_charm1.SetPxPyPzE(v[lq3[0]]->Px(),v[lq3[0]]->Py(),v[lq3[0]]->Pz(),v[lq3[0]]->E());
      v_charm2.SetPxPyPzE(v[lq3[1]]->Px(),v[lq3[1]]->Py(),v[lq3[1]]->Pz(),v[lq3[1]]->E());
      v_charm3.SetPxPyPzE(v[lq3[2]]->Px(),v[lq3[2]]->Py(),v[lq3[2]]->Pz(),v[lq3[2]]->E());
      v_charm4.SetPxPyPzE(v[lq3[3]]->Px(),v[lq3[3]]->Py(),v[lq3[3]]->Pz(),v[lq3[3]]->E());
      v_charm5.SetPxPyPzE(v[lq3[4]]->Px(),v[lq3[4]]->Py(),v[lq3[4]]->Pz(),v[lq3[4]]->E());
       for (int n_charm=0; n_charm<5 ; n_charm++)
      {
	  hPt_lq3[n_charm]->Fill( charm_quark_pt[n_charm] );
          hEta_lq3[n_charm]->Fill( charm_quark_eta[n_charm] );
          hPhi_lq3[n_charm]->Fill( charm_quark_phi[n_charm] );
      }
      cout<<"charm_quark Pt :  "<<charm_quark_pt[0]<<",  "<<charm_quark_pt[1]<<",  "<<charm_quark_pt[2]<<",  "<<charm_quark_pt[3]<<",  "<<charm_quark_pt[4]<<endl;
      
           
      v_strange1.SetPxPyPzE(v[lq4[0]]->Px(),v[lq4[0]]->Py(),v[lq4[0]]->Pz(),v[lq4[0]]->E());
      v_strange2.SetPxPyPzE(v[lq4[1]]->Px(),v[lq4[1]]->Py(),v[lq4[1]]->Pz(),v[lq4[1]]->E());
      v_strange3.SetPxPyPzE(v[lq4[2]]->Px(),v[lq4[2]]->Py(),v[lq4[2]]->Pz(),v[lq4[2]]->E());
      v_strange4.SetPxPyPzE(v[lq4[3]]->Px(),v[lq4[3]]->Py(),v[lq4[3]]->Pz(),v[lq4[3]]->E());
      v_strange4.SetPxPyPzE(v[lq4[4]]->Px(),v[lq4[4]]->Py(),v[lq4[4]]->Pz(),v[lq4[4]]->E());
      for (int n_strange=0; n_strange<5 ; n_strange++)
      {
	  hPt_lq4[n_strange]->Fill( strange_quark_pt[n_strange] );
          hEta_lq4[n_strange]->Fill( strange_quark_eta[n_strange] );
          hPhi_lq4[n_strange]->Fill( strange_quark_phi[n_strange] );
      }  
      cout<<"strange_quark Pt :  "<<strange_quark_pt[0]<<",  "<<strange_quark_pt[1]<<",  "<<strange_quark_pt[2]<<",  "<<strange_quark_pt[3]<<",  "<<strange_quark_pt[4]<<endl;
             
        // w boson kinematics:
        v_w_boson1 = v_charm1+v_strange1;
        v_w_boson2 = v_charm2+v_strange2;
        v_w_boson3 = v_charm3+v_strange3;
        v_w_boson4 = v_charm4+v_strange4;
        v_w_boson5 = v_charm5+v_strange5;
        
TLorentzVector v_wbosons[5]={v_w_boson1, v_w_boson2, v_w_boson3, v_w_boson4,  v_w_boson5};
	for (int wb=0; wb<5 ; wb++)
	{
        hPt_wboson[wb]->Fill( v_wbosons[wb].Pt() );
        hEta_wboson[wb]->Fill( v_wbosons[wb].Eta() );
        hPhi_wboson[wb]->Fill( v_wbosons[wb].Phi() );
        hM_wboson[wb]->Fill( v_wbosons[wb].M() );
        }
       
        // top quark kinematics
        v_top1= v_w_boson1+v_bq1;
        v_top2= v_w_boson2+v_bq2;
        v_top3= v_w_boson3+v_bq3;
        v_top4= v_w_boson4+v_bq4;
}
cout<<"\n";
//
TLorentzVector v_tops[5]={v_top1, v_top2, v_top3, v_top4};
	for (int t=0; t<5 ; t++)
	{
        hPt_top[t]->Fill( v_tops[t].Pt() );
        hEta_top[t]->Fill( v_tops[t].Eta() );
        hPhi_top[t]->Fill( v_tops[t].Phi() );
        hM_top[t]->Fill( v_tops[t].M() );
        }

//----------------------------*end filling the histograms*-------------------------------

      ff>>tt;
      line++;
      //if (event==100)  break;
      delete Id;
      delete Status;
      delete Mother1;
      delete Mother2;
      delete Color1;
      delete Color2;
      delete px;
      delete py;
      delete pz;
      delete E;
      delete m;
      delete lifetime;
      delete spin;
      for (int k=1 ; k<npart+1 ; delete v[k++]);

    }
    
}

  cout <<"Total number of events --> "<<event<<endl;
  TFile *rootfile = new TFile(rootfname.c_str(),"recreate");

hPt_b->Write();
hEta_b->Write();
hPhi_b->Write();

hPt_up->Write();
hEta_up->Write();
hPhi_up->Write();

hPt_down->Write();
hEta_down->Write();
hPhi_down->Write();

hPt_charm->Write();
hEta_charm->Write();
hPhi_charm->Write();

hPt_strange->Write();
hEta_strange->Write();
hPhi_strange->Write();
  
for (int i=0;i<4;i++)
    {
    hPt_bq[i]->Write();
    hEta_bq[i]->Write();
    hPhi_bq[i]->Write();
    }

for (int i=0;i<5;i++)
    {
    hPt_lq1[i]->Write();
    hEta_lq1[i]->Write();
    hPhi_lq1[i]->Write();
    }
    
for (int i=0;i<5;i++)
    {
    hPt_lq2[i]->Write();
    hEta_lq2[i]->Write();
    hPhi_lq2[i]->Write();
    }
    
for (int i=0;i<5;i++)
    {
    hPt_lq3[i]->Write();
    hEta_lq3[i]->Write();
    hPhi_lq3[i]->Write();
    }
    
for (int i=0;i<5;i++)
    {
    hPt_lq4[i]->Write();
    hEta_lq4[i]->Write();
    hPhi_lq4[i]->Write();
    }    
      
for (int i=0;i<5;i++)
  {
    hPt_wboson[i]->Write();
    hEta_wboson[i]->Write();
    hPhi_wboson[i]->Write();
    hM_wboson[i]->Write();
  }
  
for (int i=0;i<5;i++)
  {
    hPt_top[i]->Write();
    hEta_top[i]->Write();
    hPhi_top[i]->Write();
    hM_top[i]->Write();
  }

  rootfile->Close();

  cout << "Events with negative weight: " << negativeWeight << endl;
  cout<<"\n";
  exit(0);
}

float DeltaPhi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  return dphi;
}
float DeltaR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta2 - eta1;
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  float DELTAR = sqrt(pow(dphi,2)+pow(deta,2))*1.0;
  return DELTAR;
}

