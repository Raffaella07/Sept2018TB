//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void ResolutionUncorr(const char * filename){


  TFile*  file= TFile::Open(filename);
  TTree* digiTree = (TTree*)file->Get("digi");



  Float_t amp_max[54], time[54];
  Int_t k,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax;
  bool debug=false;
  Double_t max=0;
  rxmin=0;
  rxmax=0.5;


  const Int_t  nbinx=200,nbiny=300;

  rymin_l=0;
  rymax_l=200;
  rymin_r=0;
  rymax_r=200;
  tymin=0;
  tymax=100;



  Float_t x_r[nbinx],y_r[nbiny], x_l[nbinx],y_l[nbiny],rmsy_l[nbiny],rmsy_r[nbiny];
  Float_t xt[nbinx],yt[nbinx],rmsyt[nbinx];


  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);

  TF1 *fit_r = new TF1("f_r","landau",0.10,1);
  TF1 *fit_l = new TF1("f_l","landau",0.10,1);


  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  

  max=4096;
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);

    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }/*chiudo for */

  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));

  
  cout<< max << endl;

  hr_amp->Fit("f_r","R0");
  hl_amp->Fit("f_l","R0");


  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2F* h2_m= new TH2F("h2_m", "histo h2_m",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);

  for(k=0;k<digiTree->GetEntries();k++){
    
    digiTree->GetEntry(k);
    
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	h2_l->Fill(amp_max[3]/max,time[1+32]-time[0]);
	h2_r->Fill(amp_max[4]/max,time[2+32]-time[0]);
	h2_m->Fill(time[1+32]-time[2+32],(time[1+32]+time[2+32])/2-time[0]);
	if(k%500==0) cout<<time[1+32]-time[0]<<endl;
      }
    
  }//chiudo for k
  
  
  
  TH1D* histotemp_l;
  TH1D* histotemp_r;
  TH1D* histotemp_m;
  
  histotemp_l=h2_l->ProjectionY("h2_lprojY",0,nbinx);
  histotemp_r=h2_r->ProjectionY("h2_rprojY",0,nbinx);
  histotemp_m=h2_m->ProjectionY("h2_tprojY",0,nbinx);
  
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",600,550);
  TF1* g_r = new TF1("g_r","gaus",0,21);
  TF1* g_l = new TF1("g_l","gaus",0,21);
  TF1* g_m = new TF1("g_m","gaus",0,21);
  

  histotemp_m->Fit("g_m","0R");
  histotemp_l->Fit("g_l","0R");
  histotemp_r->Fit("g_r","0R");

  // hyp_r->SetParLimits(0,1,8);
  gStyle->SetOptStat("");
  gStyle->SetOptFit();
  
  histotemp_r->SetLineColor(kRed);
  histotemp_l->SetLineColor(kBlue);
  

  histotemp_m->GetXaxis()->SetTitle("t_ave-t_MCP");
  histotemp_m->GetYaxis()->SetTitle("counts");
  
  histotemp_m->Draw();
  histotemp_r->Draw("same");
  histotemp_l->Draw("same");
  
  g_l->Draw("same");
  g_r->Draw("same");
  g_m->Draw("same");

  cout<<histotemp_l->GetEntries()<<"    "<<histotemp_r->GetEntries()<<"    "<< histotemp_m->GetEntries()<<endl; 

}

