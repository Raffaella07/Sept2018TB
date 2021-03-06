TF1* fitGaus( TH1D* histo, float nSigma, bool addFunc ) {

  float mean_histo = histo->GetMean();
  float rms_histo  = histo->GetRMS();

  TF1* f1_gaus = new TF1( Form("gaus_%s", histo->GetName()), "gaus", mean_histo-rms_histo, mean_histo+rms_histo );
  f1_gaus->SetLineColor( 46 );

  histo->Fit( f1_gaus->GetName(), "RQ0" );

  float xMin_fit = f1_gaus->GetParameter(1) - nSigma*f1_gaus->GetParameter(2);
  float xMax_fit = f1_gaus->GetParameter(1) + nSigma*f1_gaus->GetParameter(2);

  f1_gaus->SetRange( xMin_fit, xMax_fit );


  Int_t n_iter = 5;

  for( Int_t i=0; i<n_iter; ++i ) { // iterative fit

    if( i==n_iter-1 && addFunc )
      histo->Fit( f1_gaus->GetName(), "RQ+" );
    else {
      histo->Fit( f1_gaus->GetName(), "RQ0" );
      xMin_fit = f1_gaus->GetParameter(1) - nSigma*f1_gaus->GetParameter(2);
      xMax_fit = f1_gaus->GetParameter(1) + nSigma*f1_gaus->GetParameter(2);
      f1_gaus->SetRange( xMin_fit, xMax_fit );
    }

  } // for iter


  return f1_gaus;

}


float getSigmaEff( TH1D* histo ) {

  float percentIntegral = 0.683;
  float Int_tegral = histo->Integral();
  // first fit to find mode:
  TF1* f1_gaus = fitGaus( histo, 1.5, false );
 
 
 
  float mode = f1_gaus->GetParameter(1);
  Int_t maxBin = histo->FindBin( mode );

  Int_t nBins = histo->GetNbinsX();
  float xMin = histo->GetXaxis()->GetXmin();
  float xMax = histo->GetXaxis()->GetXmax();

  TH1D* newHisto = new TH1D( Form("newHisto_%s", histo->GetName()), "", nBins, xMin, xMax);
  newHisto->SetBinContent( maxBin, histo->GetBinContent(maxBin) );
  newHisto->SetBinError( maxBin, histo->GetBinError(maxBin) );

  Int_t iBin = maxBin;
  Int_t delta_iBin =1;
  Int_t sign  = 1;

  float width = histo->GetBinWidth( maxBin );

  while( newHisto->Integral() < percentIntegral*Int_tegral ) {

    iBin += sign*delta_iBin;

    newHisto->SetBinContent( iBin, histo->GetBinContent(iBin) );
    newHisto->SetBinError( iBin, histo->GetBinError(iBin) );

    width += histo->GetBinWidth( iBin );

    delta_iBin += 1;
    sign *= -1;

  }

  return width/2.;

}



double subtractResoPTK( double reso ) {

  return sqrt( pow(reso,2) - pow(0.018,2) );

}




void AWtdiff(const char * filename){


  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");
  TTree* hodoTree = (TTree*)file->Get("hodo");

  unsigned int timetypes;

  digiTree->SetBranchAddress("n_timetypes",&timetypes);

  digiTree->GetEntry(3);
  Float_t amp_max[timetypes], time[timetypes],X[2],Y[2];
  Int_t k,j,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  bool debug=false;
  bool blind=true;
  Double_t max=0;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi,PTK1,AMP1,AMP2,NINO1,NINO2,CFD;
  rxmin=0;
  rxmax=1;


  const Int_t  nbinx=150,nbiny=300;

  Int_t i;
  Double_t sigma[50],erry[50],cut[50],errx[50];
  
  txmin=-20;
  txmax=20;
  
  
  Double_t x_r[nbinx],y_r[nbinx], x_l[nbinx],y_l[nbinx],rmsy_l[nbinx],rmsy_r[nbinx];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];
  
  
  Int_t nentries=digiTree->GetEntries(), counter1=0,counter2=0, counter3=0;
  Float_t Times1[nentries],Times2[nentries],Times3[nentries];
  
  for(k=0;k<digiTree->GetEntries();k++){
    Times1[k]=0;
    Times2[k]=0;
    Times3[k]=0;
  }  
  
  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1F *mcp_amp =new TH1F("mcp_amp","histomcp_amp",nbinx,0.0,1);


  TF1 *fit_r = new TF1("f_r","landau",0.15,1);
  TF1 *fit_l = new TF1("f_l","landau",0.01,0.5);

  

  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  digiTree->SetBranchAddress("CFD",&CFD);
  digiTree->SetBranchAddress("LED",&LED300);
  digiTree->SetBranchAddress("LED100",&LED100);
  digiTree->SetBranchAddress("LED50",&LED50);
  digiTree->SetBranchAddress("LED30",&LED30);
  digiTree->SetBranchAddress("PTK1",&PTK1);
  
  digiTree->SetBranchAddress("AMP1",&AMP1);
  digiTree->SetBranchAddress("AMP2",&AMP2);
  digiTree->SetBranchAddress("NINO1",&NINO1);
  digiTree->SetBranchAddress("NINO2",&NINO2);

  hodoTree->SetBranchAddress("Y",&Y);
  hodoTree->SetBranchAddress("X",&X);

  digiTree->GetEntry(3);
  
  cout << "timetypes   " << timetypes << endl;
  cout << "CFD   " << CFD << endl;
  cout << "LED300   " << LED300 << endl;
  cout << "PTK1   " << PTK1 << endl;
  cout << "NINO1   " << NINO1 << endl;
  cout << "NINO2   " << NINO2 << endl;
  cout << "AMP1   " << AMP1 << endl;
  cout << "AMP2   " << AMP2 << endl;
  
 
 
  LEDi=LED300;
  for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    
    if(time[1+LEDi]-time[0]<15 && time[1+LEDi]-time[0]>0) {
      counter1++;
      Times1[counter1]=time[1+LEDi]-time[0];
    }
    
    if(time[2+LEDi]-time[0]<15 && time[2+LEDi]-time[0]>0) {
      counter2++;
      Times2[counter2]=time[2+LEDi]-time[0];
    }
    
    if((time[1+LEDi]+time[2+LEDi])/2-time[0]<15 && (time[1+LEDi]+time[2+LEDi])/2-time[0]>-1) {
      counter3++;
      Times3[counter3]=(time[1+LEDi]+time[2+LEDi])/2-time[0];
    }
    
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

  rymin_l=mean1-1.2*rms1;
  rymax_l=mean1+0.8*rms1;
  rymin_r=mean2-1.2*rms2;
  rymax_r=mean2+0.8*rms2;
    
  

  tymin=mean3-1.2*rms3;
  tymax=mean3+0.8*rms3;

  rymin_l=0;
  rymax_l=10;
  rymin_r=0;
  rymax_r=10;
    
  
  
  tymin=0;
  tymax= 5;
  
  
  max=4096;


  for(k=0;k<digiTree->GetEntries();k++){
    if (k%100000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);

    hr_amp->Fill(amp_max[AMP2]/max);
    hl_amp->Fill(amp_max[AMP1]/max);
    mcp_amp->Fill(amp_max[PTK1]/max);

  }/*chiudo for */

  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  //cout << tmax <<endl;
  cout<< max << endl;

  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");

   TCanvas* lcheck = new TCanvas("landaucheck","landaucheck plot", 1800, 550);
  lcheck->Divide(3,1);
  lcheck->cd(1)->SetLogy();
  hr_amp->GetXaxis()->SetTitle("max.amplitude(mV)");
  hr_amp->Draw("HIST");
  fit_r->Draw("same");
  lcheck->cd(2)->SetLogy();
  hl_amp->GetXaxis()->SetTitle("max.amplitude(mV)");
  hl_amp->Draw("HIST");
  fit_l->Draw("same");
  lcheck->cd(3)->SetLogy();
  mcp_amp->GetXaxis()->SetTitle("max.amplitude(mV)");
  mcp_amp->Draw("HIST");


  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2F* h2_t= new TH2F("h2_t", "histo h2_t",nbinx/1.5,txmin,txmax,nbiny,tymin,tymax);
  TH2F* h2_time= new TH2F("h2_t", "histo h2_t",nbinx/1.5,-0.8,0.8,nbiny,tymin,tymax);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);
    hodoTree->GetEntry(k);

    if  (amp_max[PTK1]/max > 0.1 && amp_max[PTK1]/max < 0.55 && GetHodoPosition(nFibresOnX,X)>-11 && GetHodoPosition(nFibresOnY,Y)>-2 && GetHodoPosition(nFibresOnY,Y)<5){
      
	 if ((0.8*(fit_l->GetParameter(1)) < (amp_max[5]/max) && (amp_max[5]/max) < (3*fit_l->GetParameter(1))) ) {
	   	h2_l->Fill(amp_max[AMP1]/max,time[NINO1+LEDi]-time[PTK1+CFD]);
	 }
       
	 if ((0.8*(fit_r->GetParameter(1))) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1))) {
	      h2_r->Fill(amp_max[AMP2]/max,time[NINO2+LEDi]-time[0+CFD]);
	   }
       
	   	   
	       if ((0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1))) || ((0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1)))) ){
	     
	     
	     
	     h2_time->Fill(time[NINO1+LEDi]-time[NINO2+LEDi],(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[0+CFD]);
	     h2_t->Fill(GetHodoPosition(nFibresOnX,X),(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[0+CFD]);
	     //	cout << "__________________" << GetHodoPosition(nFibresOnX,X) << endl;
	    
	      
	      
	       }
       }
  }
	   
	   
	   
  TH1D* largeh_r = h2_t->ProjectionY("largeh_r",0,nbinx);
  largeh_r->GetXaxis()->SetTitle("t_{ave}-t_{MCP}");
  TCanvas* largeg = new TCanvas("large_g","plot large_g", 550, 600);
  largeh_r->Fit("gaus");
  largeh_r->Draw("same");
  
  for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
    histotemp_t=h2_time->ProjectionY("h2_tprojY",k,k);

    TF1* gaus_l = new TF1("gaus_l","gaus",histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())-3*histotemp_l->GetRMS(), histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())+3*histotemp_l->GetRMS());
    TF1* gaus_r = new TF1("gaus_r","gaus",histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())-3*histotemp_r->GetRMS(), histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())+3*histotemp_r->GetRMS());
    TF1* gaus_t = new TF1("gaus_t","gaus",histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())-3*histotemp_t->GetRMS(), histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())+3*histotemp_t->GetRMS());
      
    histotemp_l->Fit("gaus_l", "0R");
    histotemp_r->Fit("gaus_r", "0R");
    histotemp_t->Fit("gaus_t","0R");

    //xt[k]=-0.8+(Float_t)(0.8-(-0.8))/(nbinx/1.5)*k;
    xt[k]=txmin+(Float_t)(txmax-(txmin))/(nbinx/1.5)*k;
    yt[k]=histotemp_t->GetMean();
    rmsyt[k]=histotemp_t->GetMeanError();
    RMS[2][k]= histotemp_t->GetRMS();

    x_l[k]=(rxmax-rxmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    rmsy_l[k]=histotemp_l->GetMeanError();
    RMS[0][k]= histotemp_l->GetRMS();
    
    
    x_r[k]=(rxmax-rxmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    rmsy_r[k]=histotemp_r->GetMeanError();
    RMS[1][k]= histotemp_r->GetRMS();


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;
    delete gaus_l;
    delete gaus_r;
    delete gaus_t;

    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k
  
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1800,1100);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx/1.5,xt,yt,0,rmsyt);

  TF1* hyp_r = new TF1("hyp_r","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*1/(x+[5])",0.15,0.80);
  TF1* hyp_l = new TF1("hyp_l","[0]+[1]*x+[2]*x**2+[3]*x**3",0.03,0.15);

  TF1* hyp_t = new TF1("hyp_t","[1]*x**2+[2]*x+[0]",-0.1,0.65);
  
  gStyle->SetOptStat("");


    /* SetParameters*/
  hyp_l->SetParameter(0, 3.51);
  hyp_l->SetParameter(1, 5);
  hyp_l->SetParameter(2, 1.2);
  hyp_r->SetParameter(5, 0.2);
  hyp_r->SetParLimits(5, 0.1,0.25);
  /* hyp_l->SetParameter(3, -2.43e-2);
  
  hyp_r->SetParameter(0, 3.51);
 
  hyp_r->SetParameter(2, 1.2);
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
  hyp_l->DrawF1(0,0.5,"same");


  wf_c->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","0R");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  graph_r->Draw("SAMEP");
  hyp_r->DrawF1(0.15,0.85,"same");
  
  wf_c->cd(3);
  h2_t->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
  h2_t->GetXaxis()->SetTitle("t_{left}-t_{right}[ns]");
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","0");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  graph_t->Draw("SAMEP");
  //  hyp_t->DrawF1(txmin,txmax,"same");

  rymin_lc=rymin_l-hyp_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))+hyp_l->GetParameter(0);
  rymax_lc=rymax_l-hyp_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))+hyp_l->GetParameter(0);
  rymin_rc=rymin_r-hyp_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))+hyp_r->GetParameter(0);
  rymax_rc=rymax_r-hyp_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))+hyp_r->GetParameter(0);
  tymin_c=tymin-(hyp_l->Eval(0.25)-hyp_l->GetParameter(0)+hyp_r->Eval(0.25)-hyp_r->GetParameter(0))/2;
  tymax_c=tymax-(hyp_l->Eval(0.25)-hyp_l->GetParameter(0)+hyp_r->Eval(0.25)-hyp_r->GetParameter(0))/2;

  
  TH2D* hc_l= new TH2D("hc_l", "histo hc_l",nbinx,rxmin,rxmax,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_r= new TH2D("hc_r", "histo hc_r",nbinx,rxmin,rxmax,nbiny,rymin_rc,rymax_rc);
  TH2D* hc_t= new TH2D("hc_t", "histo hc_t",nbinx/1.5,txmin,txmax,nbiny,tymin_c,tymax_c);
  TH2D* hc_tdiff= new TH2D("hc_tdiff", "histo hc_tdiff",nbinx/1.5,txmin,txmax,nbiny,-10,10);
  

  
   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);
    hodoTree->GetEntry(k);
   
    if  (amp_max[PTK1]/max > 0.1 && amp_max[PTK1]/max < 0.55 && GetHodoPosition(nFibresOnX,X)>-11 && GetHodoPosition(nFibresOnY,Y)>-2 && GetHodoPosition(nFibresOnY,Y)<5){
     if ((0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1))) ) 	hc_l->Fill(amp_max[AMP1]/max,time[NINO1+LEDi]-time[0+CFD]-hyp_l->Eval(amp_max[AMP1]/max)+hyp_l->GetParameter(0));
      
      if (((0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1)))) ) hc_r->Fill(amp_max[AMP2]/max,time[NINO2+LEDi]-time[0+CFD]-hyp_r->Eval(amp_max[AMP2]/max)+hyp_r->GetParameter(0)); 
      
      if ((0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1))) || ((0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1)))) )
	 {
	   
	   hc_t->Fill(GetHodoPosition(nFibresOnX,X),(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[0+CFD]-(hyp_r->Eval(amp_max[AMP2]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[AMP1]/max)-hyp_l->GetParameter(0))/2);	
	   
	   //	cout << "__________________" << GetHodoPosition(nFibresOnX,X) << endl;
	 }//chiudo if
    }
   }//chiudo for k

    TH1D* corrh_r = hc_t->ProjectionY("corrh_r",0,nbinx);
    corrh_r->GetXaxis()->SetTitle("t_{ave}-t_{MCP}");
    TCanvas* corrg = new TCanvas("corr_g","plot corr_g", 550, 600);
    corrh_r->Fit("gaus");
    corrh_r->Draw("same");

for(k=0;k<nbinx/1.5;k++){
  TH1D* histotemp_l;
  TH1D* histotemp_r;
  TH1D* histotemp_t;
  
  histotemp_l=hc_l->ProjectionY("hc_lprojY",k,k);
  histotemp_r=hc_r->ProjectionY("hc_rprojY",k,k);
  histotemp_t=hc_t->ProjectionY("hc_tprojY",k,k);

      TF1* gaus_l = new TF1("gaus_l","gaus",histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())-3*histotemp_l->GetRMS(), histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())+3*histotemp_l->GetRMS());
    TF1* gaus_r = new TF1("gaus_r","gaus",histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())-3*histotemp_r->GetRMS(), histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())+3*histotemp_r->GetRMS());
    TF1* gaus_t = new TF1("gaus_t","gaus",histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())-3*histotemp_t->GetRMS(), histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())+3*histotemp_t->GetRMS());
      
    histotemp_l->Fit("gaus_l", "0R");
    histotemp_r->Fit("gaus_r", "0R");
    histotemp_t->Fit("gaus_t","0R");
    
    xt[k]=txmin+(Float_t)(txmax-(txmin))/(nbinx/1.5)*k;
    yt[k]=gaus_t->GetParameter(1);
    rmsyt[k]=gaus_t->GetParError(1);
    RMS[2][k]= histotemp_t->GetRMS();

    x_l[k]=(rxmax-rxmin)/nbinx*k;
    y_l[k]=gaus_l->GetParameter(1);
    rmsy_l[k]=gaus_l->GetParError(1);
    RMS[0][k]= histotemp_l->GetRMS();
    
    
    x_r[k]=(rxmax-rxmin)/nbinx*k;
    y_r[k]=gaus_r->GetParameter(1);
    rmsy_r[k]=gaus_r->GetParError(1);
    RMS[1][k]= histotemp_r->GetRMS();


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;
    delete gaus_l;
    delete gaus_r;
    delete gaus_t;
  
 

  
 }//chiudo for k


TGraphErrors* graph_lc = new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
TGraphErrors* graph_rc = new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
TGraphErrors* graph_tc = new TGraphErrors(nbinx/1.5,xt,yt,0,rmsyt);
TF1* fit_tdiff = new TF1("fit_tdiff","[0]+[1]*x+[2]*x**2",-20,20);

graph_tc->Fit("fit_tdiff","R0");

for(k=0;k<digiTree->GetEntries();k++){
     
  
  digiTree->GetEntry(k);
  hodoTree->GetEntry(k);
   if  (amp_max[PTK1]/max > 0.1 && amp_max[PTK1]/max < 0.55 && GetHodoPosition(nFibresOnX,X)>-11 && GetHodoPosition(nFibresOnY,Y)>-2 && GetHodoPosition(nFibresOnY,Y)<5){
  if ((0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1))) || ((0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1)))) ) 
    
    {
      //      cout << (time[1+LEDi]+time[2+LEDi])/2-time[0+CFD]-(hyp_r->Eval(amp_max[6]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[5]/max)-hyp_l->GetParameter(0))/2-fit_tdiff->Eval(GetHodoPosition(nFibresOnX,X))+fit_tdiff->GetParameter(0) << endl;	
	hc_tdiff->Fill(GetHodoPosition(nFibresOnX,X),(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[0+CFD]-(hyp_r->Eval(amp_max[AMP2]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[AMP1]/max)-hyp_l->GetParameter(0))/2-fit_tdiff->Eval(GetHodoPosition(nFibresOnY,Y))+fit_tdiff->GetParameter(0));	
       }
   }
 }//chiudo for k



 

 TF1* retta = new TF1("retta","[0]+[1]*x",-15,15);
 retta->SetParameter(0,0);
 retta->SetParameter(1,0);
 retta->SetParLimits(0,0,1);
 retta->SetParLimits(0,0,0.5);

 
for(k=0;k<nbinx;k++){
     
  
  TH1D* histotemp_t;
  histotemp_t=hc_tdiff->ProjectionY("hc_tprojY",k,k);
  
  TF1* gaus_t = new TF1("gaus_t","gaus",histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())-2*histotemp_t->GetRMS(), histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())+2*histotemp_t->GetRMS());
    

  histotemp_t->Fit("gaus_t","0R");
  yt[k]=gaus_t->GetParameter(1);
  rmsyt[k]=gaus_t->GetParError(1);
  // TCanvas* contol = new TCanvas("cntol","contol_plt",600,550);
  gStyle->SetOptFit(1111);
  histotemp_t->GetXaxis()->SetRangeUser(1.5,3);
  // histotemp_t->Draw();
  // gaus_t->Draw("same");

  
 
  delete gaus_t;
  delete histotemp_t;
     
     
 }//chiudo for k

TGraphErrors* graph_tcdiff = new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);
graph_tcdiff->Fit("retta","0R");

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
hc_t->GetXaxis()->SetTitle("X[mm]");
hc_t->Draw("COLZ");
graph_tc->SetMarkerStyle(8);
graph_tc->SetMarkerSize(.5);
graph_tc->Draw("P");

/*#################################################################################################################################*/ 

 TH1D* histo_cl;
 TH1D* histo_cr;
 TH1D* histo_ct;
 TH1D* histo_ctdiff;
 
histo_cl = hc_l->ProjectionY("histo_cl",0,nbinx);
histo_cr = hc_r->ProjectionY("histo_cr",0,nbinx);

histo_ct = hc_t->ProjectionY("histo_ct",0,nbinx);
histo_ctdiff = hc_tdiff->ProjectionY("histo_ctdiff",0,nbinx);

//TF1* gaus_cl = new TF1("gaus_cl","gaus");
//TF1* gaus_cr = new TF1("gaus_cr","gaus");
//TF1* gaus_ct = new TF1("gaus_ct","gaus",2,3);
TF1* gaus_ctdiff = new TF1("gaus_ctdiff","gaus", 2, 2.5);
 

 
histo_ct->SetLineColor(kBlack);
histo_cl->SetLineColor(kBlue);
histo_cr->SetLineColor(kRed);


TCanvas* tdiff = new TCanvas("tdiff","plot_tdiff",600,550);
 TLegend* l2=new TLegend(0.1,0.7,0.48,0.9);
 TF1* gaus_cl = fitGaus(histo_cl,2,false);
 TF1* gaus_cr = fitGaus(histo_cr,2,false); 
 TF1* gaus_ct = fitGaus(histo_ct,2,false); 
 gaus_ct->SetLineColor(kBlack);
 // histo_ct->Fit("gaus_ct","0R","SAME");
 //gaus_ctdiff->SetParameter(0,gaus_ct->GetParameter(0));
 // gaus_ctdiff->SetParameter(1,gaus_ct->GetParameter(1));
 //gaus_ctdiff->SetParameter(2,gaus_ct->GetParameter(2));
 gaus_ctdiff->SetLineColor(kGreen);
 histo_ctdiff->SetLineColor(kGreen);
 histo_ctdiff->Fit("gaus_ctdiff", "0R");


 histo_ctdiff->Draw();
 gaus_ctdiff->Draw("same");
 gaus_ct->Draw("same");
 histo_ct->Draw("same");


l2->SetHeader("t_{ave}-t_{MCP} distrib");
l2->AddEntry(histo_ct,"t_ave-t_MCP");
l2->AddEntry(gaus_ct,("#sigma="+to_string(gaus_ct->GetParameter(2))).c_str());
l2->AddEntry(histo_ctdiff,"t_ave-t_MCP(tdiff corr)");
l2->AddEntry(gaus_ctdiff,("#sigma="+to_string(gaus_ctdiff->GetParameter(2))).c_str());
l2->Draw();


cout << "########################### "<< gaus_ct->GetParameter(2)/gaus_ct->GetParameter(1) << "_________" << gaus_ctdiff->GetParameter(2)/gaus_ctdiff->GetParameter(1) << endl;

   
Int_t npt=13;
bool control=false;

TH1D* histotemp_t[(Int_t)nbinx/npt];
TF1* fit[(Int_t)nbinx/npt];
double SigmaEff[(Int_t)nbinx/npt];
gSystem->Exec("mkdir slicetdiff");
gSystem->Exec("cd slicetdiff");
gROOT->SetBatch(kTRUE);
TCanvas* rest_gaussine[nbinx/npt];
//   if(control)rest_gaussine->Divide(nbinx/(npt*4)+1,4);

for (i=0;i<=nbinx/npt;i++){
  cut[i] =txmin+(Float_t)(txmax-txmin)*((Float_t)i*npt)/nbinx;
 }
cut[nbinx/npt+1]=txmax-0.001;

for (i=0;i<=nbinx/npt;i++){
  // if(control)rest_gaussine->cd(i+1);
  rest_gaussine[i]  = new TCanvas(((string)"rest_gaussine"+ to_string(i)).c_str(),((string)"rest_plotgaus"+to_string(i)).c_str(),600,550);
 
  
  histotemp_t[i]=hc_t->ProjectionY(((string)("histoY"+to_string(i))).c_str(), hc_t->GetXaxis()->FindBin(cut[i]), hc_t->GetXaxis()->FindBin(cut[i+1]));
  fit[i] = fitGaus(histotemp_t[i],2,false);

  cout <<i<< "____________" << cut[i]<<"__________"<<cut[i+1]<<"__________" <<hc_t->GetXaxis()->FindBin(cut[i])<<"_____________"<<hc_t->GetXaxis()->FindBin(cut[i+1])<< endl;
  gStyle->SetOptFit(00010);
  histotemp_t[i]->SetTitle((to_string(cut[i])+"->"+to_string(cut[i+1])).c_str());
  
  histotemp_t[i]->Draw();
  fit[i]->Draw("same");
  rest_gaussine[i]->SaveAs(((string)"slicetdiff/gaus_"+to_string(i)+(string)".eps").c_str());
  sigma[i]=fit[i]->GetParameter(2);
  gPad->SetBatch(kFALSE);
  SigmaEff[i]=getSigmaEff(histotemp_t[i]);
  gPad->SetBatch(kTRUE);
  erry[i]=fit[i]->GetParError(2);
  errx[i]= (txmax-txmin)*npt/(2*nbinx);
  
 }
   //   if(!control) delete rest_gaussine;

 gROOT->SetBatch(kFALSE);
 TCanvas* rest_plot = new TCanvas("rest","rest_plot",600,550);
 TGraphErrors* rest = new TGraphErrors(nbinx/npt,cut,sigma,errx,erry);
 TGraphErrors* rest_eff = new TGraphErrors(nbinx/npt,cut,SigmaEff,errx,erry);
 TLegend* l4=new TLegend(0.5,0.8,0.5,0.8);
 l4->SetHeader("Resolutions","C");
 l4->AddEntry(rest,"#sigma (2RMS range fit)","P");
 l4->AddEntry(rest_eff,"#sigma_{eff}: 0.683*gaussArea ","P");
 rest->GetXaxis()->SetTitle("X(mm)");
 rest->GetYaxis()->SetTitle("#sigma_{t_{ave}}(ns)");
 rest->SetMarkerStyle(8);
 rest->SetMarkerSize(.8);
 rest->SetMarkerColor(kRed);
 rest_eff->SetMarkerStyle(20);
 rest_eff->SetMarkerSize(.8);
 rest->GetYaxis()->SetRangeUser(0.05,0.03);
 rest->Draw("AP");
 rest_eff->Draw("SAMEP");

 l4->Draw();
 rest_plot->SaveAs("FinalPlots/ResVsX.eps");
 
 bool ConfrontoTdiff=true;
 
   if(ConfrontoTdiff){
     TCanvas* ConfCanv = new TCanvas("confronto","",1200,800);
     gStyle->SetOptFit();
     ConfCanv->Divide(3,1);
     
     ConfCanv->cd(1);
     h2_t->Draw("COLZ");
     graph_t->Draw("P");
     
     ConfCanv->cd(2);
     hc_t->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
     hc_t->GetXaxis()->SetTitle("t_left-t_right [ns]");
     hc_t->Draw("COLZ");
     graph_tc->Draw("P");
     fit_tdiff->Draw("SAME");
     
     ConfCanv->cd(3);
     hc_tdiff->Draw("COLZ");
     graph_tcdiff->Draw("P");
     graph_tcdiff->Fit("retta");
   }


}
