//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_su(const char * filename){
  
  
  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");


  
  Float_t amp_max[54], time[54];
  Int_t k,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin,rymax;
  bool debug=false;
  Double_t max=0,tmax=0;
  rxmin=0.1;
  rxmax=0.2;
  rymin=13;
  rymax=21;
  const Int_t  nbinx=30,nbiny=100;

   Float_t x_r[nbinx],y_r[nbiny], x_l[nbinx],y_l[nbiny],rmsy_l[nbiny],rmsy_r[nbiny];
  Float_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
 

  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);

  TF1 *fit_r = new TF1("f_r","landau",0.10,1);
  TF1 *fit_l = new TF1("f_l","landau",0.10,1);
  
  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  // digiTree->SetBranchAddress("LED30",&LED30);
  //digiTree->SetBranchAddress("LED50",&LED50);
  
  for(k=0; k<digiTree->GetEntries(); k++){
    digiTree->GetEntry(k);
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    if(amp_max[3]>max) {max=amp_max[3];}
    if(amp_max[4]>max) {max=amp_max[4];}
    if(time[3]-time[4]>tmax && time[3]-time[4]<10) {tmax = time[3]-time[4];}

  }/*chiudo for */
 
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }/*chiudo for*/ 
   
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));

  
  cout<< max << endl;
  
  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");
  
  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,-0.4,6,nbiny,rxmin,rxmax); 
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,-0.4,6,nbiny,rxmin,rxmax);
  TH2F* h2_t= new TH2F("h2_t", "histo h2_t",nbinx,-0.4,1,nbinx,12,26); 
 
  for(k=0;k<digiTree->GetEntries();k++){
    
    digiTree->GetEntry(k);
    
    if (amp_max[3]/max>0.9*fit_r->GetParameter(1) &&  amp_max[3]/max<1.5*fit_r->GetParameter(1)) h2_r->Fill(time[3]-time[4],amp_max[3]/max);
    if (amp_max[4]/max>0.9*fit_l->GetParameter(1) &&  amp_max[4]/max<1.5*fit_l->GetParameter(1)){
      h2_l->Fill(time[3]-time[4],amp_max[4]/max);
   

	//h2_t->Fill((time[3]-time[4])/tmax,(time[3]+time[4])/2-time[0]);
      
      if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[4]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[3]-time[4] <<endl;
    }
    }
  
//chiudo for k

  h2_r->Draw("BOX");
  for(k=0;k<nbinx;k++){
    TH1D* histotemp_r;
    TH1D* histotemp_l;

   
    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_tprojY",k,k);
   
    y_r[k] = histotemp_r->GetMean();
    y_l[k] = histotemp_l->GetMean();

    x_r[k] = x_l[k]=-0.4 +(0.8-(-0.4))/nbinx*k;
    // y_r[k]=rxmin+(Float_t)(rxmax-rxmin)/nbiny*maxbin_r;
    // y_l[k]=rxmin+(Float_t)(rxmax-rxmin)/nbiny*maxbin_l;
    rmsy_r[k] = histotemp_r->GetRMS();
    rmsy_l[k] = histotemp_l->GetRMS();
    // cout << maxbin_l << endl;
    
    delete histotemp_r;
    delete histotemp_l;
  }
  
  
  /* for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;
  
    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
    histotemp_t=h2_t->ProjectionY("h2_tprojY",k,k);
    
    
    maxbin_l=histotemp_l->GetMaximumBin();
    maxbin_r=histotemp_r->GetMaximumBin();
    maxbin_t=histotemp_t->GetMaximumBin();
    
    xt[k]=-0.4+(Float_t)(0.8-(-0.4))/nbinx*k;
    yt[k]=12+(Float_t)(26-12)/nbinx*maxbin_t;
    rmsyt[k]=histotemp_t->GetRMS();
       
    x_l[k]=(rxmax-rxmin)/nbinx*k;
    y_l[k]=rymin+(rymax+3-rymin)/nbiny*maxbin_l;
    rmsy_l[k]=histotemp_l->GetRMS();

    x_r[k]=(rxmax-rxmin)/nbinx*k;
    y_r[k]=rymin+(rymax-rymin)/nbiny*maxbin_r;
    rmsy_r[k]=histotemp_r->GetRMS();

    
    
    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;
    
    if(k%20==0) cout << k << " / " << nbinx << endl;
    }//chiudo for k*/

  
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",600,550);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  // TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);
  //TF1* hyp_r = new TF1("hyp_r","[0] - [1]/(x**[3]-[2])",rxmin,rxmax);
  //TF1* hyp_l = new TF1("hyp_l","[0] - [1]/(x**[3]-[2])",rxmin,rxmax);
  //TF1* hyp_t = new TF1("hyp_t","[0] + [1]*x",-0.5,1);
  //  hyp_r->SetParameter(0,10);
  // hyp_l->SetParameter(0,10);

 
  // hyp_r->SetParLimits(0,1,8);
  gStyle->SetOptStat("");
  
  
  
 
  // graph_l->Fit("hyp_l","R");
  graph_l->SetMarkerSize(.8);
  graph_l->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.8);
  graph_r->SetMarkerStyle(8);
  graph_r->SetLineColor(kRed);
  graph_r->SetMarkerColor(kRed);

  graph_r->Draw("AP");
  graph_l->Draw("sameP");
   
}


