//root -l examples/ttbarH_bg_codes/fullyHadronic_kinematics_BDT.C'("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxH_bg_tttxtxwm/Events/run_01_decayed_1/tag_1_delphes_CMS_events.root")'

//root -l examples/ttbarH_bg_codes/fullyHadronic_kinematics_BDT.C'("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxwm_bg_tttxtxwm/Events/run_01_decayed_1/tag_1_delphes_CMS_events.root")'

//root -l examples/ttbarH_bg_codes/fullyHadronic_kinematics_BDT.C'("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxz_bg_tttxtxwm/Events/run_01_decayed_1/tag_1_delphes_CMS_events.root")'

#include "TH1.h"
#include "TSystem.h"
#include "TCanvas.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

// Function to calculate mass difference score for a given candidate mass and true mass
double massScore(double candidateMass, double trueMass, double tolerance) 
{
    return (abs(candidateMass - trueMass) < tolerance) ? abs(candidateMass - trueMass) : tolerance;
}

void fullyHadronic_kinematics_BDT(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create a chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create an object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchWeight = treeReader->UseBranch("Weight");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");

  // Book histograms for PT distributions
  Jet *jet[6]; // 4 leading Bjets
  Electron *electron[4];
  Muon *muon[4];
  
//____________________________________Creating BackGround Tree______________________________________________ 

float jet1_pt, jet2_pt, jet3_pt, jet4_pt, jet5_pt, jet6_pt, bjet1_pt, sumbjet_pt, Scalar_HT, missing_ET, rat_MET_HT, top_mass;
int numJets, numBJets;
float deltaEta_jets, deltaPhi_jets, deltaEta_bjets, deltaPhi_bjets, deltaR_jets, deltaR_bjets, bkg_weight;

//TTree *backgroundTree1 = new TTree("ttbarHBackground", "Tree for ttbarH background");
//TTree *backgroundTree1 = new TTree("ttbarWBackground", "Tree for ttbarW background");
TTree *backgroundTree1 = new TTree("ttbarZBackground", "Tree for ttbarZ background");

backgroundTree1->Branch("jet1_pt", &jet1_pt, "jet1_pt/F");
backgroundTree1->Branch("jet2_pt", &jet2_pt, "jet2_pt/F");
backgroundTree1->Branch("jet3_pt", &jet3_pt, "jet3_pt/F");
backgroundTree1->Branch("jet4_pt", &jet4_pt, "jet4_pt/F");
backgroundTree1->Branch("jet5_pt", &jet5_pt, "jet5_pt/F");
backgroundTree1->Branch("jet6_pt", &jet6_pt, "jet6_pt/F");

backgroundTree1->Branch("bjet1_pt", &bjet1_pt, "bjet1_pt/F");
backgroundTree1->Branch("sumbjet_pt", &sumbjet_pt, "sumbjet_pt/F");
backgroundTree1->Branch("Scalar_HT", &Scalar_HT, "Scalar_HT/F");
backgroundTree1->Branch("missing_ET", &missing_ET, "missing_ET/F");
backgroundTree1->Branch("rat_MET_HT", &rat_MET_HT, "rat_MET_HT/F");

backgroundTree1->Branch("numJets", &numJets, "numJets/I");
backgroundTree1->Branch("numBJets", &numBJets, "numBJets/I");

backgroundTree1->Branch("deltaEta_jets", &deltaEta_jets, "deltaEta_jets/F");
backgroundTree1->Branch("deltaPhi_jets", &deltaPhi_jets, "deltaPhi_jets/F");
backgroundTree1->Branch("deltaEta_bjets", &deltaEta_bjets, "deltaEta_bjets/F");
backgroundTree1->Branch("deltaPhi_bjets", &deltaPhi_bjets, "deltaPhi_bjets/F");
backgroundTree1->Branch("deltaR_jets", &deltaR_jets, "deltaR_jets/F");
backgroundTree1->Branch("deltaR_bjets", &deltaR_bjets, "deltaR_jets/F");

backgroundTree1->Branch("top_mass", &top_mass, "top_mass/F");

//backgroundTree1->Branch("ttbarH_weight", &bkg_weight, "ttbarHbkg_weight/F");
//backgroundTree1->Branch("ttbarW_weight", &bkg_weight, "ttbarWbkg_weight/F");
backgroundTree1->Branch("ttbarZ_weight", &bkg_weight, "ttbarZbkg_weight/F");

//_______________________________________end of the Tree_______________________________________________

  // Create histograms for B-Jets and Jets
  TH1F *hPt_lBJet[1];
  hPt_lBJet[0] = new TH1F("b_jet_pt_0", "leading b-jets P_{T}", 50, 0.0, 500.0);
 
  //BDT variables:  
  TH1F *hM_wboson = new TH1F("hM_wboson", "W Boson Mass", 50, 0, 160);
  TH1F *hM_top = new TH1F("hM_top", "top Mass", 50, 50, 300);
  TH1F *hScalar_HT = new TH1F("Scalar_HT", "Scalar HT", 50, 0, 4000);
  TH1F *hN_jets = new TH1F("hN_jets", "Number of Jets", 12, 3, 16);
  TH1F *hN_bjet = new TH1F("hN_bjet", "Number of b-jets", 9, 0, 8);
  TH1F *hbjet_HT = new TH1F("bjets_HT", "Sum(Pt) of b-jets", 50, 0, 1500);
  TH1F *hMissing_ET = new TH1F("Missing_ET", "Missing ET", 50, 0, 400);
  TH1F *hRatio_MET_HT = new TH1F("Ration_MET_HT", "Ratio MET ET", 30, 0, 10);
  
  TH1F *hPt_lJet[6];
  hPt_lJet[0] = new TH1F("jet_pt_0", "leading jets P_{T}", 50, 0.0, 1000.0);
  hPt_lJet[1] = new TH1F("jet_pt_1", "2nd leading jet P_{T}", 50, 0.0, 1000.0);
  hPt_lJet[2] = new TH1F("jet_pt_2", "3rd leading jet P_{T}", 50, 0.0, 500.0);
  hPt_lJet[3] = new TH1F("jet_pt_3", "4th leading jet P_{T}", 50, 0.0, 500.0);
  hPt_lJet[4] = new TH1F("jet_pt_4", "5th leading jet P_{T}", 50, 0.0, 500.0);
  hPt_lJet[5] = new TH1F("jet_pt_5", "6th leading jet P_{T}", 50, 0.0, 500.0);
  
  TH1F *hdEta_jets = new TH1F("hdEta_jets", "Jets deltaEta", 20, -3.0, 3.0);
  TH1F *hdPhi_jets = new TH1F("hdphi_jets", "Jets deltaphi", 20, -4.0, 4.0);
  TH1F *hdEta_bjets = new TH1F("hdEta_bjets", "bJets deltaEta", 10, -3.0, 3.0);
  TH1F *hdPhi_bjets = new TH1F("hdphi_bjets", "bJets deltaphi", 10, -4.0, 4.0);
  
  TH1F *hdEta_lj_lbj = new TH1F("hdEta_lj_lbj", "Leading Jets B-Jets #Delta#eta", 20, -3.0, 3.0);
  TH1F *hdEta_l1j_l1bj = new TH1F("hdEta_l1j_l1bj", "Sub-Leading Jets B-Jets #Delta#eta", 20, -3.0, 3.0);
  TH1F *hdPhi_lj_lbj = new TH1F("hdPhi_lj_lbj", "Leading Jets B-Jets #Delta#phi", 20, -4.0, 4.0);
  TH1F *hdPhi_l1j_l1bj = new TH1F("hdPhi_l1j_l1bj", "Sub-Leading Jets B-Jets #Delta#phi", 20, -4.0, 4.0);
  
  TH1F *hdR_jets = new TH1F("hdR_jets", "Jets #DeltaR", 20, 0.0, 5.);
  TH1F *hdR_bjets = new TH1F("hdR_bjets", "B-Jets #DeltaR", 20, 0.0, 5.);
  TH1F *hdR_lj_lbj = new TH1F("dR_lj_lbj", "#DeltaR", 20, 0.0, 5.);
  TH1F *hdR_l1j_l1bj = new TH1F("hdR_l1j_l1bj", "#DeltaR", 20, 0.0, 5.);
  
  TH2F *hPt_jetcorr = new TH2F("hPt_jetcorr", "Jet correlation", 50, 0, 500, 50, 0, 500);
  TH2F *hPt_bjetcorr = new TH2F("hPt_bjetcorr", "B-Jet correlation", 50, 0, 500, 50, 0, 500);
  
//-------------------------------------------------------------------------------------
    TLorentzVector v_lepton, v_lep1, v_lep2;
    Double_t leadingLeptonPT = 0.0;
    Double_t subleadingLeptonPT = 0.0;
    int n_events=0, n_cutjet=0, jetmul_cut=0, bjetmul_cut=0, n_HTcut=0, n_cutbjet=0;

//--------------------------------Loop over all events----------------------------------//
  for (Long64_t entry = 0; entry < numberOfEntries; ++entry)
  {
    int n_jets=0, n_lbj=0, n_mbj=0, n_tbj=0;
    int numPlottedBJets = 0;
    double HT_bjets=0, deltaPhi, ratio_MET_sqrtHT, HT, met;
    float dEta_ljets, dEta_lbjets, dEta_l1jets, dEta_l1bjets, dPhi_ljets, dPhi_lbjets, dPhi_l1jets, dPhi_l1bjets;
  
  // Create a vector to store combined b-jets
  std::vector<Jet*> bjets;
  std::vector<Jet*> jets;
  
    // Load selected branches with data from the specified event
    treeReader->ReadEntry(entry);
    
    //weight
    Weight *weight = (Weight*)branchWeight->At(0);
    if (!weight) continue;
    //hweight->Fill(weight->Weight);
    bkg_weight = (weight->Weight);
    
    ScalarHT *scalarHT = (ScalarHT*)branchScalarHT->At(0);
    if (!scalarHT) continue;
    hScalar_HT->Fill(scalarHT->HT);
    HT = scalarHT->HT;
    Scalar_HT = scalarHT->HT;
    
    //MET
    MissingET *missingET = (MissingET*)branchMissingET->At(0);
    if(!missingET) continue;
    hMissing_ET->Fill(missingET->MET);
    met = missingET->MET;
    missing_ET = missingET->MET;
    
    //Ratio of the Missing ET with sqrt of HT
    ratio_MET_sqrtHT = met/TMath::Sqrt(HT);
    hRatio_MET_HT->Fill(ratio_MET_sqrtHT);
    rat_MET_HT = ratio_MET_sqrtHT;
    
   //finding the correlation between the leading and subleading jets:
   if (branchJet->GetEntries()>1)
   {
   Jet *jet1 = (Jet*)branchJet->At(0);
   Jet *jet2 = (Jet*)branchJet->At(1);
   hPt_jetcorr->Fill(jet1->PT, jet2->PT);
   }
        
    //first 6 leading Jets:
    if (branchJet->GetEntriesFast() >= 6)
    {
    	for (Int_t j = 0; j < 6; ++j)
     	{
     		//if (j == 0) continue;
     		//{ 
     		//	jet[j] = (Jet *)branchJet->At(j);
     		//	hPt_lJet[j]->Fill(jet[j]->PT);
     		//}
     		jet[j] = (Jet *)branchJet->At(j);
     		hPt_lJet[j]->Fill(jet[j]->PT);
     	}
    }
    
//-----------------------reconstruction of W bosons and top quarks----------------------//
    double bestScore = std::numeric_limits<double>::max();
    
    // Loop over all possible jet combinations to form top quarks
    for (Int_t i = 0; i < branchJet->GetEntries(); ++i) 
    {
          for (Int_t j = i + 1; j < branchJet->GetEntries(); ++j) 
          {
                for (Int_t k = 0; k < branchJet->GetEntries(); ++k) 
                {
                    if (k == i || k == j) continue; // Ensure b-jet is not one of the W boson jets
                    Jet *jet1 = (Jet*)branchJet->At(i);
                    Jet *jet2 = (Jet*)branchJet->At(j);
                    Jet *bjet = (Jet*)branchJet->At(k);

                    // Form W boson and top quark candidates
                    TLorentzVector W_candidate = jet1->P4() + jet2->P4();
                    TLorentzVector top_candidate = W_candidate + bjet->P4();

                    // Calculate scores based on mass proximity to known values
                    double W_score = massScore(W_candidate.M(), 80.4, 25);
                    double top_score = massScore(top_candidate.M(), 172.5, 30);

                    // Consider each triplet independently
                    if (W_score + top_score < bestScore && bjet->BTag) 
                    {
                        bestScore = W_score + top_score;
                        hM_wboson->Fill(W_candidate.M());
                        hM_top->Fill(top_candidate.M());
                        top_mass = top_candidate.M();
                    }
                }
            }
        }
     
//-------------------Loop over all jets in the event to find the b-jets-------------------//
    for (Int_t jetIndex = 0; jetIndex < branchJet->GetEntries(); ++jetIndex)
    {
      // Access the jet
      Jet *jet = (Jet *)branchJet->At(jetIndex);
      jets.push_back(jet);
      n_jets++;     
      //loose b jets
      if (jet->BTag == 1)
      {
        n_lbj++;
        bjets.push_back(jet);

        // Plot pT of the b-jet and store it in the array
        if (numPlottedBJets < 1)
        {
          hPt_lBJet[numPlottedBJets]->Fill(jet->PT);
          ++numPlottedBJets;
        }
        //Sum of Pt of all b-jets
        HT_bjets = HT_bjets+(jet->PT);
        sumbjet_pt = HT_bjets;
      }
    }
    
//-------------------------------sort jets by Pt------------------------------//
    auto sortByPt= [](Jet* a, Jet* b){ return a->PT > b->PT; };
    std::sort(jets.begin(), jets.end(), sortByPt); 
    std::sort(jets.begin(), jets.end(), sortByPt);      
    
    if (jets.size() >=2)
    {
    float dEta_jets = (jets[0]->Eta - jets[1]->Eta);
    float dPhi_jets = (jets[0]->Phi - jets[1]->Phi);
    float dR_jets = TMath::Sqrt(pow(dEta_jets,2) + pow(dPhi_jets,2));
    hdEta_jets->Fill(dEta_jets);
    hdPhi_jets->Fill(dPhi_jets);
    hdR_jets->Fill(dR_jets);
    
    deltaEta_jets = dEta_jets;
    deltaPhi_jets = dPhi_jets;
    deltaR_jets = dR_jets;
    }
    if (bjets.size() >=2)
    {
    float dEta_bjets = (bjets[0]->Eta - bjets[1]->Eta);
    float dPhi_bjets = (bjets[0]->Phi - bjets[1]->Phi);
    float dR_bjets = TMath::Sqrt(pow(dEta_bjets,2) + pow(dPhi_bjets,2));
    hPt_bjetcorr->Fill(bjets[0]->PT, bjets[1]->PT);
    hdEta_bjets->Fill(dEta_bjets);
    hdPhi_bjets->Fill(dPhi_bjets);
    hdR_bjets->Fill(dR_bjets);
    
    deltaEta_bjets = dEta_bjets;
    deltaPhi_bjets = dPhi_bjets;
    deltaR_bjets = dR_bjets;
    
    bjet1_pt = bjets[0]->PT;
    }
    
    //leading jets and b-jets
    if (!jets.empty() && !bjets.empty())
    {
    float dEta_lj_lbj=(jets[0]->Eta - bjets[0]->Eta);
    hdEta_lj_lbj->Fill(dEta_lj_lbj);
    float dPhi_lj_lbj = (jets[0]->Phi - bjets[0]->Phi);
    float dR_lj_lbj = TMath::Sqrt(pow(dEta_lj_lbj, 2) + pow(dPhi_lj_lbj, 2));
    hdPhi_lj_lbj->Fill(dPhi_lj_lbj);
    hdR_lj_lbj->Fill(dR_lj_lbj);
    }
    
    //first subleading jets and b-jets
    if(jets.size()>=2 && bjets.size() >=2)
    {
    float dEta_l1j_l1bj = jets[1]->Eta - bjets[1]->Eta;
    hdEta_l1j_l1bj->Fill(dEta_l1j_l1bj);
    float dPhi_l1j_l1bj = (jets[1]->Phi - bjets[1]->Phi);
    float dR_l1j_l1bj = TMath::Sqrt(pow(dEta_l1j_l1bj, 2) + pow(dPhi_l1j_l1bj, 2));
    hdPhi_l1j_l1bj->Fill(dPhi_l1j_l1bj);
    hdR_l1j_l1bj->Fill(dR_l1j_l1bj);
    }
      
      if(jets.size()>0 ) jet1_pt = jet[0]->PT;
      if(jets.size()>1 ) jet2_pt = jets[1]->PT;
      if(jets.size()>2 ) jet3_pt = jets[2]->PT;
      if(jets.size()>3 ) jet4_pt = jets[3]->PT;
      if(jets.size()>4 ) jet5_pt = jets[4]->PT;
      if(jets.size()>5 ) jet6_pt = jets[5]->PT;
      
//---------------------------------applying cuts:---------------------------------//
/*      if (scalarHT->HT > 1000)
      {
      hScalar_HT->Fill(scalarHT->HT);
      Scalar_HT = (scalarHT->HT);
      n_HTcut++;
      }
      
      if(!jets.empty() && jets[0]->PT > 150)
      {
      hPt_lJet[0]->Fill(jets[0]->PT);
      jet1_pt = (jets[0]->PT);
      n_cutjet++;
      }
      
      if(bjets.size()>0)
      {
      if(bjets[0]->PT > 100)
      {
      hPt_lBJet[0]->Fill(bjets[0]->PT);
      bjet1_pt = (bjets[0]->PT);
      n_cutbjet++;
      }
      }
      
      if (n_jets >= 10)
      {
      hN_jets->Fill(n_jets);
      numJets = n_jets;
      jetmul_cut++;
      }
      
      if (n_lbj >= 2)
      {
      hN_bjet->Fill(n_lbj);
      numBJets = n_lbj;
      bjetmul_cut++;
      }
*/      
      hN_jets->Fill(n_jets);
      numJets = n_jets;
      hN_bjet->Fill(n_lbj);
      numBJets = n_lbj;
      hbjet_HT->Fill(HT_bjets);
      
      //cout<<"HT_bjets: "<<HT_bjets<<endl;
      //cout<<"loose bjets: "<<n_lbj<<endl;
      //cout<<"medium bjets: "<<n_mbj<<endl;
      //cout<<"tight bjets: "<<n_tbj<<endl;
      //cout<<"deltaPhi: "<<deltaPhi<<endl;
      cout<<"Evenets: "<<entry<<endl;
            
      backgroundTree1->Fill();
         
  }//end of for loop
  
  //-------------------------------output----------------------------------------//
  
    //cout<<"No. of events after cut on Jet Pt: "<<n_cutjet<<endl;
    //cout<<"No. of events after cut on bJet Pt: "<<n_cutbjet<<endl;
    //cout<<"No. of events after cut on Scalar HT: "<<n_HTcut<<endl;
    //cout<<"Jet multiplicity after cut: "<<jetmul_cut<<endl;
    //cout<<"bJet multiplicity after cut: "<<bjetmul_cut<<endl;
  
/* 
//--------------------------------------------------------------------------------------------      
  //TFile *rootFile = new TFile("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxH_bg_tttxtxwm/Events/run_01_decayed_1/ttbarH_bg_bdt_variables.root", "RECREATE");
  //TFile *rootFile = new TFile("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxH_bg_tttxtxwm/Events/run_01_decayed_1/ttbarH_bg_kvariables_cut.root", "RECREATE");
  
  //TFile *rootFile = new TFile("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxwm_bg_tttxtxwm/Events/run_01_decayed_1/ttbarwm_bg_bdt_variables.root", "RECREATE");
  //TFile *rootFile = new TFile("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxwm_bg_tttxtxwm/Events/run_01_decayed_1/ttbarwm_bg_kvariables_cut.root", "RECREATE");
  
  //TFile *rootFile = new TFile("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxz_bg_tttxtxwm/Events/run_01_decayed_1/ttbarz_bg_bdt_variables.root", "RECREATE");
  TFile *rootFile = new TFile("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxz_bg_tttxtxwm/Events/run_01_decayed_1/ttbarz_bg_kvariables_cut.root", "RECREATE");


//---------------------------------*Histograms Plotting*----------------------------------------

  TCanvas *canvas1 = new TCanvas("myCanvas1", "BDT Variables leading jets", 800, 600);
  canvas1->Divide(2,3);
    
  canvas1->cd(1);
  hPt_lJet[0]->GetXaxis()->SetTitle("Pt_{j1} [GeV]");
  hPt_lJet[0]->GetYaxis()->SetTitle("Events");
  hPt_lJet[0]->Draw("hist");
  hPt_lJet[0]->SetStats(0);
  
  canvas1->cd(2);
  hPt_lJet[1]->GetXaxis()->SetTitle("Pt_{j2} [GeV]");
  hPt_lJet[1]->GetYaxis()->SetTitle("Events");
  hPt_lJet[1]->Draw("hist same");
  hPt_lJet[1]->SetStats(0);
  
  canvas1->cd(3);
  hPt_lJet[2]->GetXaxis()->SetTitle("Pt_{j3} [GeV]");
  hPt_lJet[2]->GetYaxis()->SetTitle("Events");
  hPt_lJet[2]->Draw("hist same");
  hPt_lJet[2]->SetStats(0);
  
  canvas1->cd(4);
  hPt_lJet[3]->GetXaxis()->SetTitle("Pt_{j4} [GeV]");
  hPt_lJet[3]->GetYaxis()->SetTitle("Events");
  hPt_lJet[3]->Draw("hist same");
  hPt_lJet[3]->SetStats(0);
  
  canvas1->cd(5);
  hPt_lJet[4]->GetXaxis()->SetTitle("Pt_{j5} [GeV]");
  hPt_lJet[4]->GetYaxis()->SetTitle("Events");
  hPt_lJet[4]->Draw("hist same");
  hPt_lJet[4]->SetStats(0);
  
  canvas1->cd(6);
  hPt_lJet[5]->GetXaxis()->SetTitle("Pt_{j6} [GeV]");
  hPt_lJet[5]->GetYaxis()->SetTitle("Events");
  hPt_lJet[5]->Draw("hist same");
  hPt_lJet[5]->SetStats(0);
  canvas1->Draw();
  //canvas1->SaveAs("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxH_bg_tttxtxwm/analysis_code/bdt_plots/six_leading_subleading_jets.png");
  canvas1->Update();
  
  TCanvas *canvas2 = new TCanvas("myCanvas2", "BDT Variables N-jets and N-bjets", 800, 600);
  canvas2->Divide(1,2);
  canvas2->cd(1);
  hN_jets->GetXaxis()->SetTitle("N_{j}");
  hN_jets->GetYaxis()->SetTitle("Events");
  hN_jets->Draw("hist");
  hN_jets->SetStats(0);
  canvas2->cd(2);
  hN_bjet->GetXaxis()->SetTitle("N_{bj}");
  hN_bjet->GetYaxis()->SetTitle("Events");
  hN_bjet->Draw("hist");
  hN_bjet->SetStats(0);
  canvas2->Draw();
  //canvas2->SaveAs("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxH_bg_tttxtxwm/analysis_code/bdt_plots/Njets_Nbjets.png");
  canvas2->Update();
  
  TCanvas *canvas3 = new TCanvas("myCanvas3", "BDT Variables", 800, 600);
  canvas3->Divide(2,2);
  canvas3->cd(1);
  hScalar_HT->GetXaxis()->SetTitle("ScalarHT [GeV]");
  hScalar_HT->GetYaxis()->SetTitle("Normalized");
  hScalar_HT->Draw("hist");
  hScalar_HT->SetStats(0);
  canvas3->cd(2);
  hbjet_HT->GetXaxis()->SetTitle("#sum P_{t} [GeV]");
  hbjet_HT->GetYaxis()->SetTitle("Normalized");
  hbjet_HT->Draw("hist");
  hbjet_HT->SetStats(0);
  canvas3->cd(3);
  hM_wboson->GetXaxis()->SetTitle("m_{w} [GeV]");
  hM_wboson->GetYaxis()->SetTitle("Normalized");
  hM_wboson->Draw("hist");
  hM_wboson->SetStats(0);
  canvas3->cd(4);
  hM_top->GetXaxis()->SetTitle("m_{t} [GeV]");
  hM_top->GetYaxis()->SetTitle("Normalized");
  hM_top->Draw("hist");
  hM_top->SetStats(0);
  canvas3->Draw();
  //canvas3->SaveAs("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxH_bg_tttxtxwm/analysis_code/bdt_plots/reco_masses_scalarHT.png");
  canvas3->Update();

  TCanvas *canvas4 = new TCanvas("myCanvas4", "BDT Variables bjets", 800, 600);
  hPt_lBJet[0]->GetXaxis()->SetTitle("leading bj P_{t} [GeV]");
  hPt_lBJet[0]->GetYaxis()->SetTitle("Events");
  hPt_lBJet[0]->Draw("hist");
  hPt_lBJet[0]->SetStats(0); 
  canvas4->Draw();
  //canvas4->SaveAs("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxH_bg_tttxtxwm/analysis_code/bdt_plots/leading_bjets.png");
  canvas4->Update();
  
  TCanvas *canvas5 = new TCanvas("myCanvas5", "BDT Variables jets", 800, 600);
  canvas5->Divide(2,2);
  canvas5->cd(1);
  hdEta_jets->GetXaxis()->SetTitle("#Delta#eta");
  hdEta_jets->GetYaxis()->SetTitle("Events");
  hdEta_jets->Draw("Hist");
  hdEta_jets->SetStats(0);
  canvas5->cd(2);
  hdPhi_jets->GetXaxis()->SetTitle("#Delta#phi");
  hdPhi_jets->GetYaxis()->SetTitle("Events");
  hdPhi_jets->Draw("hist");
  hdPhi_jets->SetStats(0);
  canvas5->Draw();
  canvas5->cd(3);
  hdEta_bjets->GetXaxis()->SetTitle("#Delta#eta");
  hdEta_bjets->GetYaxis()->SetTitle("Events");
  hdEta_bjets->Draw();
  hdEta_bjets->SetStats(0);
  canvas5->cd(4);
  hdPhi_bjets->GetXaxis()->SetTitle("#Delta#phi");
  hdPhi_bjets->GetYaxis()->SetTitle("Events");
  hdPhi_bjets->Draw();
  hdPhi_bjets->SetStats(0);
  canvas5->Draw();
  //canvas5->SaveAs("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxH_bg_tttxtxwm/analysis_code/bdt_plots/delta_eta_phi_jets_bjets.png");
  canvas5->Update();
    
  TCanvas *canvas6 = new TCanvas("myCanvas6", "Correlation of Jets and BJets", 800, 600);
  canvas6->Divide(2,1);
  canvas6->cd(1);
  hPt_jetcorr->GetXaxis()->SetTitle("P_{t} [GeV]");
  hPt_jetcorr->GetYaxis()->SetTitle("P_{t} [GeV]");
  hPt_jetcorr->Draw("COLZ");
  hPt_jetcorr->SetStats(0);
  canvas6->cd(2);
  hPt_bjetcorr->GetXaxis()->SetTitle("P_{t} [GeV]");
  hPt_bjetcorr->GetYaxis()->SetTitle("P_{t} [GeV]");
  hPt_bjetcorr->Draw("COLZ");
  hPt_bjetcorr->SetStats(0);
  //canvas6->SaveAs("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_pp_ttxH_bg_tttxtxwm/analysis_code/bdt_plots/correlation_jets_bjets.png"); 
  canvas6->Update();
  
  TCanvas *canvas7 = new TCanvas("myCanvas7", "#DeltaR Jets", 800, 600);
  hdR_jets->GetXaxis()->SetTitle("#DeltaR_{j,j1}");
  hdR_jets->GetYaxis()->SetTitle("Events");
  hdR_jets->Draw("hist");
  hdR_jets->SetStats(0);
  
  TCanvas *canvas8 = new TCanvas("myCanvas8", "#DeltaR B-Jets", 800, 600);
  hdR_bjets->GetXaxis()->SetTitle("#DeltaR_{bj,bj1}");
  hdR_bjets->GetYaxis()->SetTitle("Events");
  hdR_bjets->Draw("hist");
  hdR_bjets->SetStats(0);
  
  TCanvas *canvas9 = new TCanvas("myCanvas9", "#DeltaEta", 800, 600);
  hdEta_lj_lbj->GetXaxis()->SetTitle("#DeltaEta_{lj,lbj}");
  hdEta_lj_lbj->GetYaxis()->SetTitle("Events");
  hdEta_lj_lbj->Draw("hist");
  hdEta_lj_lbj->SetStats(0);
  
  TCanvas *canvas10 = new TCanvas("myCanvas10", "#DeltaEta", 800, 600);
  hdEta_l1j_l1bj->GetXaxis()->SetTitle("#DeltaEta_{l1j,l1bj}");
  hdEta_l1j_l1bj->GetYaxis()->SetTitle("Events");
  hdEta_l1j_l1bj->Draw("hist");
  hdEta_l1j_l1bj->SetStats(0);
  
  TCanvas *canvas11 = new TCanvas("myCanvas11", "#DeltaPhi", 800, 600);
  hdPhi_lj_lbj->GetXaxis()->SetTitle("#DeltaPhi_{lj,lbj}");
  hdPhi_lj_lbj->GetYaxis()->SetTitle("Events");
  hdPhi_lj_lbj->Draw("hist");
  hdPhi_lj_lbj->SetStats(0);
  
  TCanvas *canvas12 = new TCanvas("myCanvas12", "#DeltaPhi", 800, 600);
  hdPhi_l1j_l1bj->GetXaxis()->SetTitle("#DeltaPhi_{l1j,l1bj}");
  hdPhi_l1j_l1bj->GetYaxis()->SetTitle("Events");
  hdPhi_l1j_l1bj->Draw("hist");
  hdPhi_l1j_l1bj->SetStats(0);
  
  TCanvas *canvas13 = new TCanvas("myCanvas13", "#DeltaR", 800, 600);
  hdR_lj_lbj->GetXaxis()->SetTitle("#DeltaR_{lj,lbj}");
  hdR_lj_lbj->GetYaxis()->SetTitle("Events");
  hdR_lj_lbj->Draw("hist");
  hdR_lj_lbj->SetStats(0);
  
  TCanvas *canvas14 = new TCanvas("myCanvas14", "#DeltaR", 800, 600);
  hdR_l1j_l1bj->GetXaxis()->SetTitle("#DeltaR_{l1j,l1bj}");
  hdR_l1j_l1bj->GetYaxis()->SetTitle("Events");
  hdR_l1j_l1bj->Draw("hist");
  hdR_l1j_l1bj->SetStats(0);
  
  TCanvas *canvas15 = new TCanvas("myCanvas15", "#MissingET", 800, 600);
  hMissing_ET->GetXaxis()->SetTitle("MET");
  hMissing_ET->GetYaxis()->SetTitle("Events");
  hMissing_ET->Draw("hist");
  hMissing_ET->SetStats(0);
  
  TCanvas *canvas16 = new TCanvas("myCanvas16", "#Ratio MET HT", 800, 600);
  hRatio_MET_HT->GetXaxis()->SetTitle("MET/\\sqrt{HT}");
  hRatio_MET_HT->GetYaxis()->SetTitle("Events");
  hRatio_MET_HT->Draw("hist");
  hRatio_MET_HT->SetStats(0);
 
  hPt_lJet[0]->Write();
  hPt_lJet[1]->Write();
  hPt_lJet[2]->Write();
  hPt_lJet[3]->Write();
  hPt_lJet[4]->Write();
  hPt_lJet[5]->Write();
  
  hN_jets->Write();
  hN_bjet->Write();
  hScalar_HT->Write();
  hbjet_HT->Write();
  hPt_lBJet[0]->Write();
  
  hdEta_jets->Write();
  hdPhi_jets->Write();
  hdEta_bjets->Write();
  hdPhi_bjets->Write();
  
  hdEta_lj_lbj->Write();
  hdEta_l1j_l1bj->Write();
  hdPhi_lj_lbj->Write();
  hdPhi_l1j_l1bj->Write();
  
  hdR_jets->Write();
  hdR_bjets->Write();
  hdR_lj_lbj->Write();
  hdR_l1j_l1bj->Write();
  
  hMissing_ET->Write();
  hRatio_MET_HT->Write();
  
  hPt_jetcorr->Write();
  hPt_bjetcorr->Write();
  hM_wboson->Write();
  hM_top->Write();
     
  rootFile->Close(); 
*/
  TFile *outputFile1 = new TFile("/home/msahil/work/root-6.20.08/tutorials/tmva/files/ttbarZ_tree.root", "RECREATE");
  backgroundTree1->Write();
  outputFile1->Close();

}


