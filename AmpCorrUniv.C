//To run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3
//commento


void AmpCorrUniv(const char * filename){
  
  
  TFile*  file= TFile::Open(filename);
  TTree* digiTree = (TTree*)file->Get("digi");

  

  Float_t amp_max[225], time[225];
  Int_t k,j,maxbin_l,maxbin_r,maxbin_t,i;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  bool debug=false;
  bool blind=true;
  Double_t max=0;
  Int_t LED300,LED100,LED50,LED30,LED;
  Int_t LEDi,CFD,AMP1,AMP2,NINO1,NINO2;
  rxmin=0;
  rxmax=0.3;
  
  
  
  const Int_t  nbinx=350,nbiny=400;
  
  
  
  
    Double_t x_r[nbinx],y_r[nbiny], x_l[nbinx],y_l[nbiny],rmsy_l[nbiny],rmsy_r[nbiny];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];
  const Int_t nentries=digiTree->GetEntries();
  Int_t counter1=0, counter2=0, counter3 =0;
  Float_t Times1[nentries],Times2[nentries],Times3[nentries];

  
  
  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1F *mcp_amp =new TH1F("mcp_amp","histomcp_ampl",nbinx,0.0,1);
  
  
  TF1 *fit_r = new TF1("f_r","landau",0.01,0.3);
  TF1 *fit_l = new TF1("f_l","landau",0.055,1);
  
  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  digiTree->SetBranchAddress("LED",&LED300);
  digiTree->SetBranchAddress("CFD",&CFD);
  digiTree->SetBranchAddress("AMP1",&AMP1);
  digiTree->SetBranchAddress("AMP2",&AMP2);
  digiTree->SetBranchAddress("NINO1",&NINO1);
  digiTree->SetBranchAddress("NINO2",&NINO2);
  
  /* digiTree->SetBranchAddress("LED100",&LED100);
  digiTree->SetBranchAddress("LED50",&LED50);
  digiTree->SetBranchAddress("LED30",&LED30);
  */
   digiTree->GetEntry(3);
  LEDi=LED300;
  /*
  for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    //cout<<"HERE"<<endl;
    if(time[1+LEDi]-time[0]<150 && time[1+LEDi]-time[0]>50) {
      counter1++;
      Times1[counter1]=time[1+LEDi]-time[0];
      cout << Times1[k] << endl;}
    //  else{Times1[k]=0;}
    if(time[2+LEDi]-time[0]<150 && time[2+LEDi]-time[0]>50) {
      counter2++;
      Times2[counter2]=time[2+LEDi]-time[0];}
    // else{Times2[k]=0;}  
    if((time[1+LEDi]+time[2+LEDi])/2-time[0]<150 && (time[1+LEDi]+time[2+LEDi])/2-time[0]>50) {
      counter3++;
      Times3[counter3]=(time[1+LEDi]+time[2+LEDi])/2-time[0];}
    // else{Times3[k]=0;}
  }
  
  Double_t mean1=TMath::Mean(counter1,Times1);
  Double_t rms1=TMath::RMS(counter1,Times1);
  cout<<mean1<<"_________"<<rms1<<endl;
  Double_t mean2=TMath::Mean(counter2,Times2);
  Double_t rms2=TMath::RMS(counter2,Times2);
  cout<<mean2<<"_________"<<rms2<<endl;
  Double_t mean3=TMath::Mean(counter3,Times3);
  Double_t rms3=TMath::RMS(counter3,Times3);
  cout<<mean3<<"_________"<<rms3<<endl;

  rymin_l=mean1-1.1*rms1;
  rymax_l=mean1+5*rms1;
  rymin_r=mean2-1.1*rms2;
  rymax_r=mean2+5*rms2;
    
    

  tymin=mean3-1.1*rms3;
  tymax=mean3+5*rms3;
  */
  rymin_l=0;
  rymax_l=10;
  rymin_r=0;
  rymax_r=10;
    
  
  
  tymin=0;
  tymax= 5;
  
  


  txmin=-0.8;
  txmax=0.8;
  
  
  
  max=4096;

  for(k=0;k<digiTree->GetEntries();k++){
    if (k%3000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[AMP1]/max);
    hl_amp->Fill(amp_max[AMP2]/max);
  }//chiudo for k
   
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  cout<< max << endl;
 
  
  hr_amp->Fit("f_r","R0");
  hl_amp->Fit("f_l","R0");
  cout << 0.8*fit_l->GetParameter(1) << "____" << 3*fit_l\
	   ->GetParameter(1) << endl;  

 
  TCanvas* lcheck = new TCanvas("landaucheck","landaucheck plot", 1200, 550);
  lcheck->Divide(2,1);
  lcheck->cd(1)->SetLogy();
  hr_amp->Draw("HIST");
  fit_r->Draw("same");
  lcheck->cd(2)->SetLogy();
  hl_amp->Draw("HIST");
  fit_l->Draw("same");
  

  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2F* h2_t= new TH2F("h2_t", "histo h2_t",nbinx,txmin,txmax,nbiny,tymin,tymax);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1))/* && amp_max[0]/max > 0.1 && amp_max[0]/max < 0.75*/)
    
      {
	/*	cout << "##############################"<< endl;
	cout << time[1]-time[0] << endl;
	cout << time[2+LEDi]-time[0] << endl;
	cout << "##############################"<< endl; */
	h2_l->Fill(amp_max[AMP1]/max,time[NINO1+LEDi]-time[CFD]);
	h2_r->Fill(amp_max[AMP2]/max,time[NINO2+LEDi]-time[CFD]);
	h2_t->Fill((time[NINO1+LEDi]-time[NINO2+LEDi]),(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[CFD]);

      }//chiudo if      

  }//chiudo for k
 
   
  for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
    histotemp_t=h2_t->ProjectionY("h2_tprojY",k,k);
   
    xt[k]=txmin+(Float_t)(txmax-(txmin))/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    rmsyt[k]=histotemp_t->GetMeanError();
    RMS[2][k]= histotemp_t->GetRMS();

    x_l[k]=rxmin+(rxmax-rxmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    rmsy_l[k]=histotemp_l->GetMeanError();
    RMS[0][k]= histotemp_l->GetRMS();
    
    
    x_r[k]=rxmin+(rxmax-rxmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    rmsy_r[k]=histotemp_r->GetMeanError();
    RMS[1][k]= histotemp_r->GetRMS();


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k
  
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1800,1100);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);

  TF1* hyp_r = new TF1("hyp_r","[0]+[2]*log(x)",0.05,0.25);
  TF1* hyp_l = new TF1("hyp_l","[0]+[2]*log(x)",0.0,0.09);

  TF1* hyp_t = new TF1("hyp_t","[1]*x**2+[2]*x+[0]",-0.8,0.65);
  
  gStyle->SetOptStat("");


  /* SetParameters
   */    hyp_l->SetParameter(0, 3.51);
   hyp_l->SetParameter(1, 3);
   /* hyp_l->SetParameter(2, 1.2);
  //hyp_l->SetParameter(3, -2.43e-2);
   
  hyp_r->SetParameter(0, 3);
  hyp_r->SetParameter(1, -9e-2);
  hyp_r->SetParameter(2, -1e-1);
  
  */
  
  wf_c->Divide(3,2);

  wf_c->cd(1);

  h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");

   h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_l","0R");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("SAMEP");
  hyp_l->DrawF1(0,1,"same");


  wf_c->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","0R");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  graph_r->Draw("SAMEP");
  hyp_r->DrawF1(0,1,"same");
  
  wf_c->cd(3);
  h2_t->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
  h2_t->GetXaxis()->SetTitle("t_left-t_right [ns]");
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","RL0");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  graph_t->Draw("SAMEP");
  hyp_t->DrawF1(txmin,txmax,"same");


  rymin_lc=rymin_l-hyp_l->Eval(0.25)+hyp_l->GetParameter(0);
  rymax_lc=rymax_l-hyp_l->Eval(0.25)+hyp_l->GetParameter(0);
  rymin_rc=rymin_r-hyp_r->Eval(0.25)+hyp_r->GetParameter(0);
  rymax_rc=rymax_r-hyp_r->Eval(0.25)+hyp_r->GetParameter(0);
  tymin_c=tymin-hyp_t->Eval(0.2)+hyp_t->GetParameter(0);
  tymax_c=tymax-hyp_t->Eval(0.2)+hyp_t->GetParameter(0);



  TH2F* hc_l= new TH2F("hc_l", "histo hc_l",nbinx,rxmin,rxmax,nbiny,rymin_lc,rymax_lc);
  TH2F* hc_r= new TH2F("hc_r", "histo hc_r",nbinx,rxmin,rxmax,nbiny,rymin_rc,rymax_rc);
  TH2F* hc_t= new TH2F("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,-5,5);


  
   for(k=0;k<digiTree->GetEntries();k++){

     digiTree->GetEntry(k);
   
     if (0.8*(fit_l->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_l->GetParameter(1))/* && amp_max[0]/max > 0.1 && amp_max[0]/max < 0.75*/)  
       {
	hc_l->Fill(amp_max[AMP1]/max,time[NINO1+LEDi]-time[PTK1+CFD]-hyp_l->Eval(amp_max[AMP1]/max)+hyp_l->GetParameter(0));
        hc_r->Fill(amp_max[AMP2]/max,time[NINO2+LEDi]-time[PTK1+CFD]-hyp_r->Eval(amp_max[AMP2]/max)+hyp_r->GetParameter(0));
	hc_t->Fill((time[NINO1+LEDi]-time[NINO2+LEDi]-(hyp_r->Eval(amp_max[AMP1]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[AMP1]/max)+hyp_l->GetParameter(0))),(time[NINO1+LEDi]+time[2+LEDi])/2-time[32]-(hyp_r->Eval(amp_max[6]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[7]/max)-hyp_l->GetParameter(0))/2);
	
       }//chiudo if

  }//chiudo for k

   for(k=0;k<nbinx;k++){
     TH1D* histotemp_l;
     TH1D* histotemp_r;
     TH1D* histotemp_t;
     
     histotemp_l=hc_l->ProjectionY("hc_lprojY",k,k);
     histotemp_r=hc_r->ProjectionY("hc_rprojY",k,k);
     histotemp_t=hc_t->ProjectionY("hc_tprojY",k,k);
     
     
     yt[k]=histotemp_t->GetMean();
     RMS[2][k]= histotemp_t->GetMeanError();
     
     y_l[k]=histotemp_l->GetMean();
     RMS[0][k]= histotemp_l->GetMeanError();
     
     y_r[k]=histotemp_r->GetMean();
     RMS[1][k]= histotemp_r->GetMeanError();
     
     
     delete histotemp_l;
     delete histotemp_r;
     delete histotemp_t;
     
     
   }//chiudo for k
   

   TGraphErrors* graph_lc = new TGraphErrors(nbinx-1,x_l,y_l,0,RMS[0]);
   TGraphErrors* graph_rc = new TGraphErrors(nbinx-1,x_r,y_r,0,RMS[1]);
   TGraphErrors* graph_tc = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);
   TF1* corr_tc= new TF1("corr_tc","[0]*x+[1]",-0.8,0.8);

    wf_c->cd(4);
     hc_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
    hc_l->GetXaxis()->SetTitle("max.amplitude [mV]");
    hc_l->Draw("COLZ");
    graph_lc->SetMarkerStyle(8);
    graph_lc->SetMarkerSize(.5);
    graph_lc->Draw("P");



   wf_c->cd(5);

    hc_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
   hc_r->GetXaxis()->SetTitle("max.amplitude [mV]");
   hc_r->Draw("COLZ");
   graph_rc->SetMarkerStyle(8);
   graph_rc->SetMarkerSize(.5);
   graph_rc->Draw("P");


   wf_c->cd(6);

   hc_t->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
   hc_t->GetXaxis()->SetTitle("t_left-t_right [ns]");
   


   hc_t->Draw("COLZ");
   graph_tc->SetMarkerStyle(8);
   graph_tc->SetMarkerSize(.5);
   graph_tc->Draw("P");
   graph_tc->Fit("corr_tc","R");
   corr_tc->Draw("SAME");



   TH1D* histo_cl;
   TH1D* histo_cr;
   TH1D* histo_ct;
   TF1* gaus_cl = new TF1("gaus_cl","gaus",-2.5,-0.7);
   TF1* gaus_cr = new TF1("gaus_cr","gaus",-2.5,-0.7);
   TF1* gaus_ct = new TF1("gaus_ct","gaus",1,3);
   histo_cl = hc_l->ProjectionY("histo_cl",0,nbinx);
   histo_cr = hc_r->ProjectionY("histo_cr",0,nbinx);
   histo_ct = hc_t->ProjectionY("histo_ct",0,nbinx);


   histo_ct->SetLineColor(kBlack);
   histo_cl->SetLineColor(kBlue);
   histo_cr->SetLineColor(kRed);
   gaus_ct->SetLineColor(kBlack);
   gaus_cl->SetLineColor(kBlue);
   gaus_cr->SetLineColor(kRed);
   
   TCanvas* canv= new TCanvas();
   TLegend* l2 = new TLegend();

   histo_ct->GetXaxis()->SetTitle("t_{ave}-t_{MCP} (ns)");
   histo_ct->GetYaxis()->SetTitle("counts");
   histo_ct->SetAxisRange(0,histo_ct->GetBinContent(histo_ct->GetMaximumBin())+100,"Y");
   //   histo_cr->Draw();
   histo_ct->Draw("SAME");
   // histo_cl->Draw("SAME");
   histo_ct->Fit("gaus_ct", "R","same");
   //  histo_cr->Fit("gaus_cr", "","same");
   // histo_cl->Fit("gaus_cl", "","same");

   l2->AddEntry(histo_ct,((string) "#sigma_{t}^{ave}="+to_string(sqrt(gaus_ct->GetParameter(2)*gaus_ct->GetParameter(2)-0.015*0.015))+(string)"ns").c_str(), "l");
   // l2->AddEntry(histo_cl,((string) "#sigma_{t}^{left}="+to_string(sqrt(gaus_cl->GetParameter(2)*gaus_cl->GetParameter(2)-0.015*0.015))+(string)"ns").c_str(),"l");
   // l2->AddEntry(histo_cr,((string) "#sigma_{t}^{right}="+to_string(sqrt(gaus_cr->GetParameter(2)*gaus_cr->GetParameter(2)-0.015*0.015))+(string)"ns").c_str(),"l");
   l2->Draw();
   
}
