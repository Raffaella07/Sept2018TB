//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_cut(const char * filename, TH1D* amp_l, TH1D* amp_r,Double_t* Mipl,Double_t* sMipl, Double_t* Mipr, Double_t* sMipr){
  
  
  TFile *  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree * digiTree = (TTree*)file->Get("digi");

  Float_t amp_max[54];
  Int_t k;
  Double_t max=4096;
  
  TH1D *hr_amp =new TH1D("hr_amp","histos_ampr",500,0.0,1);
  TH1D *hl_amp =new TH1D("hl_amp","histos_ampl",500,0.0,1);
  TF1 *fit_r = new TF1("f_r","landau",0.06,1);
  TF1 *fit_l = new TF1("f_l","landau",0.05,1);
  TH1D *hr_cut =new TH1D("hr_cut","histos_cut",500,0.0,1);
  TH1D *hl_cut =new TH1D("hl_cut","histos_cut ",500,0.0,1);
  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  

  max=4096;
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%3000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }//chiudo for k
   
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  cout<< max << endl;
  
  hr_amp->Fit("f_r","R0");
  hl_amp->Fit("f_l","R0");
  
  for(k=0;k<digiTree->GetEntries();k++){
    
    digiTree->GetEntry(k);
    
    if (0.8*fit_l->GetParameter(1) < amp_max[4]/max && amp_max[4]/max < 3*fit_l->GetParameter(1)){

      hr_cut->Fill(amp_max[3]/max);
      hl_cut->Fill(amp_max[4]/max);}
    
  }//chiudo for k
   
  hr_cut->Scale(1/(hr_cut->Integral()));
  hl_cut->Scale(1/(hl_cut->Integral()));
  hr_cut->Scale(hr_amp->Integral(0.8*fit_r->GetParameter(1)*500, 3*fit_r->GetParameter(1)*500)/hr_amp->Integral());
  hl_cut->Scale(hl_amp->Integral(0.8*fit_l->GetParameter(1)*500, 3*fit_l->GetParameter(1)*500)/hl_amp->Integral());

  TCanvas* wf_c =new TCanvas("wf","Plot wf",1200,550);
  wf_c->Clear();
  
  wf_c->Divide(2,1);

  wf_c->cd(1)->SetLogy();
  //gStyle->SetOptFit();
  hr_amp->SetLineColor(4);
  hr_amp->GetXaxis()->SetTitle("max.amplitude [mV]");
  hr_amp->GetYaxis()->SetTitle("counts(normalized)");
  hr_amp->Draw("HISTO"); 
  // fit_r->DrawF1(0,1,"same");
  hr_cut->SetLineColor(kBlack);
  // hr_cut->Draw("HISTO same");
  
  wf_c->cd(2)->SetLogy();
  hl_amp->SetLineColor(4);
  hl_amp->GetXaxis()->SetTitle("max.amplitude [mV]");
  hl_amp->GetYaxis()->SetTitle("counts(normalized)");
  //gStyle->SetOptFit();
  hl_amp->Draw("HISTO");
  hl_cut->SetLineColor(kBlack);
  // fit_l->DrawF1(0,1,"same");
  // hl_cut->Draw("HISTO same"); 
 
  hl_amp->Copy(*amp_l);
  hr_amp->Copy(*amp_r);
  *Mipl=fit_l->GetParameter(1);
  *sMipl=fit_l->GetParError(1);
 *Mipr=fit_r->GetParameter(1);
  *sMipr=fit_r->GetParError(1);
}

void Darkbkg(){
  Int_t ncurrents=4,i;
  Double_t DCR[ncurrents],erry_l[ncurrents],erry_r[ncurrents];
  Double_t Mipl[2][ncurrents], Mipr[2][ncurrents],sMipl[2][ncurrents],sMipr[2][ncurrents];
  TCanvas* super[2][ncurrents];
  TPad* wf[2][ncurrents];
  TPad* bkg[2][ncurrents];
  DCR[0]=50;
  DCR[1]=500;
  DCR[2]=1000;
  DCR[3]=1940;
   
  TH1D *amp_l[ncurrents];
  TH1D *amp_r[ncurrents];
  TH1D *amp_ldcr[ncurrents];
  TH1D *amp_rdcr[ncurrents];
  TH1D *bgk_l[ncurrents];
  TH1D *bgk_r[ncurrents];
  TF1 *fitr;
  TF1 *fitl;
  TF1 *fitr_dcr[ncurrents];
  TF1 *fitl_dcr[ncurrents];
 
 
  
  string filename = "Pd/ConfT100-B72-1.2.root";
  cout << "hereee" << endl;
   for(i=0;i<ncurrents;i++){
     amp_ldcr[i] =new TH1D(((string)"amp_L_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),((string)"amp_L"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),500,0.0,1);
     amp_rdcr[i] =new TH1D(((string)"amp_R_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),((string)"amp_R_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),500,0.0,1);
     amp_l[i] =new TH1D(((string)"amp_L_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),((string)"amp_L_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),500,0.0,1);
     amp_r[i] =new TH1D(((string)"amp_R_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),((string)"amp_R_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),500,0.0,1);
     bgk_l[i] =new TH1D(((string)"bgk_L_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),((string)"bgk_L_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),500,0.0,1);
     bgk_r[i] =new TH1D(((string)"bgk_R_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),((string)"bgk_R_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),500,0.0,1);

     fitl= new TF1(((string)"fitl"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),"landau",0.03,1);
     fitr= new TF1(((string)"fitr"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),"landau",0.05,1);
     fitl_dcr[i]= new TF1(((string)"fitl_dcr"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),"landau",0.05,1);
     fitr_dcr[i]= new TF1(((string)"fitr_dcr"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),"landau",0.05,1);

    
     gROOT->SetBatch(kTRUE);
     
     super[0][i] =new TCanvas(((string)"superL"+to_string((Int_t)DCR[i])).c_str(),((string)"Plot super"+to_string((Int_t)DCR[i])).c_str(),800,1200);
    
     wf[0][i] =new TPad(((string)"superwf"+to_string((Int_t)DCR[i])).c_str(),((string)"Plot super"+to_string((Int_t)DCR[i])).c_str(),0.0,0.2,1,1);
     bkg[0][i] =new TPad(((string)"superb"+to_string((Int_t)DCR[i])).c_str(),((string)"Plot super"+to_string((Int_t)DCR[i])).c_str(),0.0,0.0,1,0.2);
     wf[1][i] =new TPad(((string)"superwf"+to_string((Int_t)DCR[i])).c_str(),((string)"Plot super"+to_string((Int_t)DCR[i])).c_str(),0.0,0.2,1,1);
     bkg[1][i] =new TPad(((string)"superb"+to_string((Int_t)DCR[i])).c_str(),((string)"Plot super"+to_string((Int_t)DCR[i])).c_str(),0.0,0.0,1,0.2);
     wf[0][i]->Draw();
     bkg[0][i]->Draw();
    

     plotWF_cut(filename.c_str(),amp_l[i],amp_r[i], &Mipl[0][i],&sMipl[0][i],&Mipr[0][i],&sMipr[0][i]);

     filename = "DCR/ConfT500-B72-1.2DCR"+to_string((Int_t)DCR[i])+".root";
     
     plotWF_cut(filename.c_str(),amp_ldcr[i],amp_rdcr[i],&Mipl[1][i],&sMipl[1][i],&Mipr[1][i],&sMipr[1][i]);
     bgk_l[i]->Add(amp_ldcr[i],amp_l[i],1,-1);
     bgk_r[i]->Add(amp_rdcr[i],amp_r[i],1,-1);
    
    

     amp_l[i]->Fit(((string)"fitl"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),"0R");
     amp_r[i]->Fit(((string)"fitr"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),"0R");
     amp_ldcr[i]->Fit(((string)"fitl_dcr"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),"0R");
     amp_rdcr[i]->Fit(((string)"fitr_dcr"+to_string((Int_t)DCR[i])+(string)"#muA").c_str(),"0R");

     wf[0][i]->cd()->SetLogy();
     
     gPad->SetFrameBorderSize(0);
     gPad->SetBottomMargin(0.000000000000);
     gPad->SetTopMargin(0.1);
     gPad->SetLeftMargin(0.15);
     gPad->SetRightMargin(0.01);
     gStyle->SetOptStat(0);
     TLegend* l1 = new TLegend();
     l1->AddEntry(amp_ldcr[i],"ampL_dcr","l");
     amp_ldcr[i]->SetTitle(((string)"amp_L_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str());
    
     
     amp_l[i]->SetLineColor(kRed);
     fitl->SetLineColor(kRed);
     fitl_dcr[i]->SetLineColor(kBlue);
     l1->AddEntry(amp_l[i],"ampL","l");
     amp_ldcr[i]->Draw("hist");
     fitl_dcr[i]->Draw("same");
     fitl->Draw("same");
     amp_l[i]->Draw("HISTsame");
     
     l1->Draw();
     
     cout << Mipl[0][i] << endl;
     erry_l[i]=sqrt(pow(sMipl[0][i],2)+ pow(sMipl[1][i],2));
     erry_r[i]=sqrt(pow(sMipr[0][i],2)+ pow(sMipr[1][i],2));
     // Mipl[0][i]-=Mipl[1][i];
     //Mipr[0][i]-=Mipr[1][i];
    
     bkg[0][i]->cd();
     gPad->SetTopMargin(0.0000000000000000);
     gPad->SetRightMargin(0.01);
     gPad->SetBottomMargin(0.2);
     gPad->SetLeftMargin(0.15);
      
   
     bgk_l[i]->GetYaxis()->SetLabelSize(0.1);
     bgk_l[i]->GetXaxis()->SetLabelSize(0.1);
     bgk_l[i]->GetXaxis()->SetTitleSize(0.1);
     bgk_l[i]->GetXaxis()->SetTitle("amp_max(mV)");
     bgk_l[i]->Draw();

     super[0][i]->SaveAs(((string)"controlplots/superL"+to_string((Int_t)DCR[i])+(string)".eps").c_str());
     super[1][i] =new TCanvas(((string)"superR"+to_string((Int_t)DCR[i])).c_str(),((string)"Plot super"+to_string((Int_t)DCR[i])).c_str(),800,1200);
     wf[1][i]->Draw();
     bkg[1][i]->Draw(); 
     wf[1][i]->cd()->SetLogy();
     
     gPad->SetFrameBorderSize(0);
     gPad->SetBottomMargin(0.000000000000);
     gPad->SetTopMargin(0.1);
     gPad->SetRightMargin(0.01);
     gPad->SetLeftMargin(0.15);
     gStyle->SetOptStat(0);
     TLegend* l6 = new TLegend();
     l6->AddEntry(amp_ldcr[i],"ampR_dcr","l");
     amp_rdcr[i]->SetTitle(((string)"amp_R_"+to_string((Int_t)DCR[i])+(string)"#muA").c_str());
    
     
     amp_r[i]->SetLineColor(kRed);
     fitr->SetLineColor(kRed);
     fitr_dcr[i]->SetLineColor(kBlue);
     l1->AddEntry(amp_r[i],"ampr","l");
     amp_rdcr[i]->Draw("hist");
     fitr_dcr[i]->Draw("same");
     fitr->Draw("same");
     amp_r[i]->Draw("HISTsame");
     
     l6->Draw();
     
        
     bkg[1][i]->cd();
     gPad->SetTopMargin(0.0000000000000000);
     gPad->SetRightMargin(0.01);
     gPad->SetBottomMargin(0.2);
      gPad->SetLeftMargin(0.15);
   
     bgk_r[i]->GetYaxis()->SetLabelSize(0.1);
     bgk_r[i]->GetXaxis()->SetLabelSize(0.1);
     bgk_r[i]->GetXaxis()->SetTitleSize(0.1);
     bgk_r[i]->GetXaxis()->SetTitle("amp_max(mV)");
     bgk_r[i]->Draw();

     super[1][i]->SaveAs(((string)"controlplots/superR"+to_string((Int_t)DCR[i])+(string)".eps").c_str()); 
     gROOT->SetBatch(kFALSE);

     
   }
   
 
  
   TGraphErrors* leftcurrents = new TGraphErrors(4,DCR,Mipl[1],0,erry_l);
   TGraphErrors* rightcurrents = new TGraphErrors(4,DCR,Mipr[1],0,erry_r);

   TCanvas* currents = new TCanvas("currents","plot_currents",600,550);
   TLegend* l2= new TLegend();
   
    leftcurrents->GetYaxis()->SetRangeUser(0.0,0.3);
   leftcurrents->GetXaxis()->SetTitle("SIPM (dark) current (#muA)");
   leftcurrents->GetYaxis()->SetTitle("Peak_{Mip}^{DCR}(mV)");
   leftcurrents->SetMarkerSize(.8);
   rightcurrents->SetMarkerSize(.8);
   leftcurrents->SetMarkerStyle(8);
   rightcurrents->SetMarkerStyle(8);
   leftcurrents->SetMarkerColor(kBlue);
   rightcurrents->SetMarkerColor(kRed);

   l2->AddEntry(leftcurrents, " SIPM Left", "P");
   l2->AddEntry(rightcurrents, " SIPM Right", "P");
   
   leftcurrents->Draw("AP");
   rightcurrents->Draw("P");
   l2->Draw();
  
  
  

}

