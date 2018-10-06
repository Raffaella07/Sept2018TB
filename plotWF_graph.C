
using namespace std;
void plotWF_graph(const char * filename){
  
  
  TFile *  file= TFile::Open(filename);
  TTree * WFTree = (TTree*)file->Get("wf");
  TTree * digiTree = (TTree*)file->Get("digi");

  Float_t amp_max[28];
  Int_t k;
  Double_t max=0;
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1200,550);
  wf_c->Divide(2,1);

  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",500,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",500,0.0,1);

  digiTree->SetBranchAddress("amp_max",&amp_max);
  
  /*  for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    if(amp_max[6]>max) {max=amp_max[6];}
    if(amp_max[7]>max) {max=amp_max[7];}
  }//chiudo for k
  */
   max=4096;
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[6]/max);
    hl_amp->Fill(amp_max[7]/max);
  }//chiudo for k   
  
  //hl_amp->SetLineColor();
  
  cout<< max << endl;

  wf_c->cd(1)->SetLogy();
  hr_amp->GetXaxis()->SetTitle("max.amplitude [mV]");
  hr_amp->GetYaxis()->SetTitle("counts");
  hr_amp->DrawNormalized();
  wf_c->cd(2)->SetLogy();
  hl_amp->GetXaxis()->SetTitle("max.amplitude [mV]");
  hl_amp->GetYaxis()->SetTitle("counts");
  hl_amp->DrawNormalized();


  TString histoname = "";
  histoname.Append("histo");
  histoname.Append(filename);

  
 TFile* f2 = new TFile(histoname.Data(),"recreate");
  f2->cd();
  hr_amp->Write();
  hl_amp->Write();
  f2->Close();
  // wf_c->cd(2)->SetLogy();
  //hr_amp->Draw();
  //hl_amp->Draw("same");
  
}

