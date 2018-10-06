//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_cut(const char * filename){
  
  
  TFile *  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree * digiTree = (TTree*)file->Get("digi");

  Float_t amp_max[225];
  Int_t k;
  Double_t max=0;
  
  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",500,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",500,0.0,1);
  TF1 *fit_r = new TF1("f_r","landau",0.15,1);
  TF1 *fit_l = new TF1("f_l","landau",0.03,1);
  TH1F *hr_cut =new TH1F("hr_cut","histos_cut",500,0.0,1);
  TH1F *hl_cut =new TH1F("hl_cut","histos_cut ",500,0.0,1);
  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  

  max=4096;
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%3000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[6]/max);
    hl_amp->Fill(amp_max[5]/max);
  }//chiudo for k
   
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  cout<< max << endl;
  
  hr_amp->Fit("f_r","R0");
  hl_amp->Fit("f_l","R0");
  cout << 0.8*fit_l->GetParameter(1) << "____" << 3*fit_l\
	   ->GetParameter(1) << endl;
  for(k=0;k<digiTree->GetEntries();k++){
    
    digiTree->GetEntry(k);
    
    if (0.8*fit_r->GetParameter(1) < amp_max[6]/max && amp_max[6]/max < 3*fit_r->GetParameter(1)){
      
      hr_cut->Fill(amp_max[6]/max);
    }
    if (0.8*fit_l->GetParameter(1) < amp_max[5]/max && amp_max[5]/max < 3*fit_l->GetParameter(1)){
      hl_cut->Fill(amp_max[5]/max);}
  }
//chiudo for k
   
   hr_cut->Scale(1/(hr_cut->Integral()));
  hl_cut->Scale(1/(hl_cut->Integral()));
  hr_cut->Scale(hr_amp->Integral(0.8*fit_r->GetParameter(1)*500, 3*fit_r->GetParameter(1)*500)/hr_amp->Integral());
  hl_cut->Scale(hl_amp->Integral(0.8*fit_l->GetParameter(1)*500, 3*fit_l->GetParameter(1)*500)/hl_amp->Integral());
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1200,550);
  wf_c->Clear();
  
  wf_c->Divide(2,1);

  wf_c->cd(1)->SetLogy();
  gStyle->SetOptFit();
  hr_amp->SetLineColor(4);
  hr_amp->GetXaxis()->SetTitle("max.amplitude [mV]");
  hr_amp->GetYaxis()->SetTitle("counts");
  hr_amp->Draw("HISTO"); 
  fit_r->DrawF1(0,1,"same");
  hr_cut->SetLineColor(kBlack);
  hr_cut->Draw("HISTO same");
  
  wf_c->cd(2)->SetLogy();
  hl_amp->SetLineColor(4);
  hl_amp->GetXaxis()->SetTitle("max.amplitude [mV]");
  hl_amp->GetYaxis()->SetTitle("counts");
  //gStyle->SetOptFit();
  hl_amp->Draw("HISTO");
  hl_cut->SetLineColor(kBlack);
  fit_l->DrawF1(0,1,"same");
  hl_cut->Draw("HISTO same"); 
  
}


