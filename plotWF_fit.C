
void plotWF_fit(const char * filename){
  
  
  TFile *  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree * digiTree = (TTree*)file->Get("digi");

  Float_t amp_max[54];
  Int_t k;
  Double_t max=0;
  
 
 

  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",500,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",500,0.0,1);
  TF1 *fit_r = new TF1("f_r","landau",0.137,2);
  TF1 *fit_l = new TF1("f_l","landau",0.13,2);
  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  
  for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    if(amp_max[3]>max) {max=amp_max[3];}
    if(amp_max[4]>max) {max=amp_max[4];}
  }//chiudo for k
  
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%3000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }//chiudo for k
   
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));

  TCanvas* wf_c =new TCanvas("wf","Plot wf",1200,550);
  wf_c->Clear();
  cout<< max << endl;

  gStyle->SetOptFit();
  wf_c->Divide(2,1);
  
  wf_c->cd(1)->SetLogy();

  hr_amp->SetLineColor(kBlue);
  hr_amp->GetXaxis()->SetTitle("max.amplitude [mV]");
  hr_amp->GetYaxis()->SetTitle("counts");
  hr_amp->Draw("HISTO");
  hr_amp->Fit("f_r","R");
  fit_r->DrawF1(0,1,"same");

  wf_c->cd(2)->SetLogy();

  hl_amp->SetLineColor(kBlue);
  hl_amp->GetXaxis()->SetTitle("max.amplitude [mV]");
  hl_amp->GetYaxis()->SetTitle("counts");
  hl_amp->Draw("HISTO");
  hl_amp->Fit("f_l","R");
  fit_l->DrawF1(0,1,"same");
  /*
  TString histoname = "";
  histoname.Append("histo");
  histoname.Append(filename);

  
 TFile* f2 = new TFile(histoname.Data(),"RECREATE");
  f2->cd();
  hr_amp->Write();
  hl_amp->Write();
  f2->Close();
  // wf_c->cd(2)->SetLogy();
  //hr_amp->Draw();
  //hl_amp->Draw("same");
  */
}

