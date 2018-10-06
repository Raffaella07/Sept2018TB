//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3
float GetHodoPosition( Int_t nFibres[2], Float_t var[2] ) {

  Int_t nFibresMax = 10;

  bool showering0 = nFibres[0]>=nFibresMax;
  bool showering1 = nFibres[1]>=nFibresMax;

  bool bad0 = nFibres[0]==0 || var[0]<-900.;
  bool bad1 = nFibres[1]==0 || var[1]<-900.;

  float pos = -999.;

  if( showering0 || showering1 ) pos = -999.;  // showering (at least one plane)

  else if( bad0 && bad1 ) pos = -999.;  // empty (both planes have 0 fibres)

  else {

    if( bad0 ) pos = var[1];

    else if( bad1 ) pos = var[0];

    else pos = 0.5*(var[0]+var[1]);

  } // else

  return pos;

}




void plotWF_time(const char * filename){


  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");
 TTree* hodoTree = (TTree*)file->Get("hodo");

  unsigned int timetypes;

  digiTree->SetBranchAddress("n_timetypes",&timetypes);

  digiTree->GetEntry(3);
  Float_t amp_max[timetypes], time[timetypes],Y[2],X[2];


 
  Int_t k,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  bool debug=false;
  Double_t max=0;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi,PTK1,AMP1,AMP2,NINO1,NINO2,CFD;
 
  Int_t nentries=digiTree->GetEntries();
  Float_t Times1[nentries],Times2[nentries],Times3[nentries];
  rxmin=0;
  rxmax=0.5;


  const Int_t  nbinx=200,nbiny=500;

  

  txmin=-20;
  txmax=20;


  Float_t x_r[nbinx],y_r[nbiny], x_l[nbinx],y_l[nbiny],rmsy_l[nbiny],rmsy_r[nbiny];
  Float_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Float_t x_rx[nbinx],y_rx[nbiny], x_lx[nbinx],y_lx[nbiny],rmsy_lx[nbiny],rmsy_rx[nbiny];
  Float_t xtx[nbinx],ytx[nbinx],rmsytx[nbinx];


  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);

  TF1 *fit_r = new TF1("f_r","landau",0.15,1);
  TF1 *fit_l = new TF1("f_l","landau",0.01,0.5);


  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  digiTree->SetBranchAddress("LED",&LED300);
  digiTree->SetBranchAddress("LED100",&LED100);
  digiTree->SetBranchAddress("LED50",&LED50);
  digiTree->SetBranchAddress("LED30",&LED30);
  digiTree->SetBranchAddress("PTK1",&PTK1);
  digiTree->SetBranchAddress("CFD",&CFD);
  
  digiTree->SetBranchAddress("AMP1",&AMP1);
  digiTree->SetBranchAddress("AMP2",&AMP2);
  digiTree->SetBranchAddress("NINO1",&NINO1);
  digiTree->SetBranchAddress("NINO2",&NINO2);

  hodoTree->SetBranchAddress("X",&X);
  hodoTree->SetBranchAddress("Y",&Y);
  hodoTree->SetBranchAddress( "nFibresOnX", nFibresOnX );
  
  hodoTree->SetBranchAddress( "nFibresOnY", nFibresOnY );


  digiTree->GetEntry(3);
  LEDi=LED300;

  for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    //cout<<"HERE"<<endl;
    if(time[1+LEDi]-time[0]<50 && time[1+LEDi]-time[0]>0) {Times1[k]=time[1+LEDi]-time[0];}
    else{Times1[k]=Times2[k-1];}
    if(time[2+LEDi]-time[0]<50 && time[2+LEDi]-time[0]>0) {Times2[k]=time[2+LEDi]-time[0];}
    else{Times2[k]=Times2[k-1];}  
    if((time[1+LEDi]+time[2+LEDi])/2-time[0]<50 && (time[1+LEDi]+time[2+LEDi])/2-time[0]>-10) {Times3[k]=(time[1+LEDi]+time[2+LEDi])/2-time[0];}    else{Times3[k]=Times3[k-1];}
    
    
  }
  Double_t mean1=TMath::Mean(nentries,Times1)-1.2;
  Double_t rms1=TMath::RMS(nentries,Times1);
  cout<<mean1<<"_________"<<rms1<<endl;
  Double_t mean2=TMath::Mean(nentries,Times2)-1.2;
  Double_t rms2=TMath::RMS(nentries,Times2);
  cout<<mean2<<"_________"<<rms2<<endl;
  Double_t mean3=TMath::Mean(nentries,Times3)-1.2;
  Double_t rms3=TMath::RMS(nentries,Times3);
  cout<<mean3<<"_________"<<rms3<<endl;

  rymin_l=mean1-0.5*rms1;
  rymax_l=mean1+0.5*rms1;
  rymin_r=mean2-0.5*rms2;
  rymax_r=mean2+0.5*rms2;
    
  

  tymin=mean3-0.5*rms3;
  tymax=mean3+0.5*rms3;

  rymin_l=0;
  rymax_l=10;
  rymin_r=0;
  rymax_r=10;
    
  
  
  tymin=-5;
  tymax= 5;
  
  max=4096;

  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);

    hr_amp->Fill(amp_max[AMP2]/max);
    hl_amp->Fill(amp_max[AMP1]/max);
  }/*chiudo for */

  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));

  //cout << tmax <<endl;
  cout<< max << endl;

  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");


  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2F* h2_t= new TH2F("h2_t", "histo h2_t",nbinx,txmin,txmax,nbiny,tymin,tymax);

   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);
    hodoTree->GetEntry(k);
       if (amp_max[0]/max>0.08 && amp_max[0]/max < 0.55){
	 if ((0.8*(fit_l->GetParameter(1)) < (amp_max[5]/max) && (amp_max[5]/max) < (3*fit_l->GetParameter(1))) ) 	h2_l->Fill(amp_max[5]/max,time[1+LEDi]-time[0+CFD]);
	 
	 if (((0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1)))) ) h2_r->Fill(amp_max[AMP2]/max,time[NINO2+LEDi]-time[0+CFD]);
	 
	 if ((0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1))) || ((0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1)))) )
	   {
	     
	     
	     h2_t->Fill(GetHodoPosition(nFibresOnY,Y),(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[0+CFD]);
	     //	cout << "__________________" << GetHodoPosition(nFibresOnX,X) << endl;
	   }//chiudo if
       }
    
   }//chiudo for k

    
 
 for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
    histotemp_t=h2_t->ProjectionY("h2_tprojY",k,k);

    TF1* gaus_l = new TF1("gaus_l","gaus",histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())-3*histotemp_l->GetRMS(), histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())+3*histotemp_l->GetRMS());
    TF1* gaus_r = new TF1("gaus_r","gaus",histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())-3*histotemp_r->GetRMS(), histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())+3*histotemp_r->GetRMS());
    TF1* gaus_t = new TF1("gaus_t","gaus",histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())-3*histotemp_t->GetRMS(), histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())+3*histotemp_t->GetRMS());
      
    histotemp_l->Fit("gaus_l", "0R");
    histotemp_r->Fit("gaus_r", "0R");
    histotemp_t->Fit("gaus_t","0R");
    
    xt[k]=txmin+(Float_t)(txmax-(txmin))/nbinx*k;
    yt[k]=gaus_t->GetParameter(1);
    rmsyt[k]=gaus_t->GetParError(1);
 

    x_l[k]=(rxmax-rxmin)/nbinx*k;
    y_l[k]=gaus_l->GetParameter(1);
    rmsy_l[k]=gaus_l->GetParError(1);
 
    
    
    x_r[k]=(rxmax-rxmin)/nbinx*k;
    y_r[k]=gaus_r->GetParameter(1);
    rmsy_r[k]=gaus_r->GetParError(1);
 


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;
    delete gaus_l;
    delete gaus_r;
    delete gaus_t;

    if(k%20==0) cout << k << " / " << nbinx << endl;
 }//chiudo for k



  TCanvas* wf_c =new TCanvas("wf","Plot wf",1800,550);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);
 
  
  // TF1* hyp_r = new TF1("hyp_r","[0]-[1]*log(x+[2])",0.135,0.35);
  //  TF1* hyp_l = new TF1("hyp_l","[0]-[1]*log(x+[2])",0.11,0.35);

  TF1* hyp_r = new TF1("hyp_r","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*1/(x+[5])",0.15,0.80);
  TF1* hyp_l = new TF1("hyp_l","[0]+[1]*x+[2]*x**2+[3]*x**3",0.03,0.15);
  TF1* hyp_t = new TF1("hyp_t","[0]*x+[1]",0.7,txmax);


    /* SetParameters*/
  hyp_l->SetParameter(0, 8.51);
  hyp_l->SetParameter(1, 5);
  hyp_l->SetParameter(2, 1.2);
  // hyp_l->SetParameter(3, -2.43e-2);
  
  hyp_r->SetParameter(0, 7);
  hyp_r->SetParameter(1, 5);
  hyp_r->SetParameter(2, 1.2);
  


  wf_c->Divide(3,2);
  wf_c->cd(1);
  h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
  h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_l","0RL");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("P");
  hyp_l->Draw("same");
  
  
  
  wf_c->cd(2);
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","0");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  graph_r->Draw("P");
  hyp_r->Draw("same");

  wf_c->cd(3);
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","R");
  
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  graph_t->Draw("P");
  

   
  TH2D* h3_l= new TH2D("h3_l", "histo h3_l",nbinx,-0.8,0.8,nbiny,tymin,tymax);
  TH2D* h3_r= new TH2D("h3_r", "histo h3_r",nbinx,-0.8,0.8,nbiny,tymin,tymax);
  TH2D* h3_t= new TH2D("h3_t", "histo h3_t",nbinx,-0.8,0.8,nbiny,tymin,tymax);
  TH2D* h3Y_l= new TH2D("h3Yl", "histo h3Y_l",nbinx,txmin,txmax,nbiny,tymin,tymax);
  TH2D* h3Y_r= new TH2D("h3Y_r", "histo h3Y_r",nbinx,txmin,txmax,nbiny,tymin,tymax);
  TH2D* h3Y_t= new TH2D("h3Y_t", "histo h3Y_t",nbinx,txmin,txmax,nbiny,tymin,tymax);
  TH1D* h = new TH1D("h","h plot",nbinx,txmin,txmax);
  TH1D*distrib_r = new TH1D("h","h plot",nbinx,txmin,100);
  TH1D*distrib_l = new TH1D("h","h plot",nbinx,txmin,100);
  TH2D* scat= new TH2D("h_scat", "histo h_scat",nbinx,txmin,txmax,nbinx,txmin,txmax);
  TH2D* e= new TH2D("aa","aaa", nbinx,txmin,txmax,nbinx,txmin,txmax);
  for(k=0;k<digiTree->GetEntries();k++){
        digiTree->GetEntry(k);
       hodoTree->GetEntry(k);
  
    h->Fill(time[NINO1+LEDi]-time[NINO2+LEDi]);
    distrib_l->Fill(time[NINO1+LEDi]-time[PTK1+CFD]);
    distrib_r->Fill(time[NINO2+LEDi]-time[PTK1+CFD]);
    scat->Fill(time[NINO1+LEDi]-time[PTK1+CFD],time[NINO2+LEDi]-time[PTK1+CFD]);
    e->Fill(GetHodoPosition(nFibresOnX,X),GetHodoPosition(nFibresOnY,Y),amp_max[AMP1]/max);
   if  (amp_max[PTK1]/max > 0.1 && amp_max[PTK1]/max < 0.55 && GetHodoPosition(nFibresOnX,X)>-11 && GetHodoPosition(nFibresOnY,Y)>-2 && GetHodoPosition(nFibresOnY,Y)<5){
      if (0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1)))
      {
	h3_l->Fill(-(time[NINO1+LEDi]-time[NINO2+LEDi]),time[NINO1+LEDi]-time[PTK1+CFD]);
      }
    if (0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1)) ){
      h3_r->Fill(-(time[NINO1+LEDi]-time[NINO2+LEDi]),time[NINO2+LEDi]-time[PTK1+CFD]);
    }
    if ((0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1))) || ( 0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_l->GetParameter(1)) ) ) h3_t->Fill(-(time[NINO1+LEDi]-time[NINO2+LEDi]),(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[PTK1+CFD]);//chiudo for k
    // }
  }
  }
   for(k=0;k<digiTree->GetEntries();k++){

       digiTree->GetEntry(k);
       hodoTree->GetEntry(k);
     
       if  (amp_max[PTK1]/max > 0.1 && amp_max[PTK1]/max < 0.55 && GetHodoPosition(nFibresOnX,X)>-11 && GetHodoPosition(nFibresOnY,Y)>-2 && GetHodoPosition(nFibresOnY,Y)<5){
      if (0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1)))
      {
	 h3Y_l->Fill(GetHodoPosition(nFibresOnX,X),time[NINO1+LEDi]-time[PTK1+CFD]);
      }
    if (0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1)) ){
	h3Y_r->Fill(GetHodoPosition(nFibresOnX,X),time[NINO2+LEDi]-time[PTK1+CFD]);
    }
    if ((0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1))) || ( 0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_l->GetParameter(1)) ) ) h3Y_t->Fill(GetHodoPosition(nFibresOnX,X),(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[PTK1+CFD]);//chiudo for k
    }
   }


  
  TCanvas* hp = new TCanvas("HP","HPplot",1200,550);
  
  hp->Divide(2,1);
  hp->cd(1);
  h3_l->Draw();
  hp->cd(2)->SetLogz();
  e->GetXaxis()->SetTitle("X[mm]");
  e->GetYaxis()->SetTitle("Y[mm]");
  e->Draw("COLZ");
    
  TCanvas* hscatt = new TCanvas("Hscat","Hscatplot",600,550);
  scat->Draw();
  
  for(k=0;k<nbinx;k++){
  TH1D* histotemp_l;
  TH1D* histotemp_r;
  TH1D* histotemp_t;
  
  histotemp_l=h3Y_l->ProjectionY("hc_lprojY",k,k);
  histotemp_r=h3Y_r->ProjectionY("hc_rprojY",k,k);
  histotemp_t=h3Y_t->ProjectionY("hc_tprojY",k,k);

  TF1* gaus_l = new TF1("gaus_l","gaus",histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())-3*histotemp_l->GetRMS(), histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())+3*histotemp_l->GetRMS());
  TF1* gaus_r = new TF1("gaus_r","gaus",histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())-3*histotemp_r->GetRMS(), histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())+3*histotemp_r->GetRMS());
  TF1* gaus_t = new TF1("gaus_t","gaus",histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())-3*histotemp_t->GetRMS(), histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())+3*histotemp_t->GetRMS());


        
    // histotemp_l->Fit("gaus_l", "0R");
    // histotemp_r->Fit("gaus_r", "0R");
    // histotemp_t->Fit("gaus_t","0R");
    
    xt[k]=txmin+(Float_t)(txmax-(txmin))/(nbinx)*k;
    yt[k]=gaus_t->GetParameter(1);
    rmsyt[k]=gaus_t->GetParError(1);
 

    x_l[k]=txmin+(txmax-txmin)/(nbinx)*k;
    y_l[k]=gaus_l->GetParameter(1);
    rmsy_l[k]=gaus_l->GetParError(1);
    
    
    x_r[k]=txmin+(txmax-txmin)/(nbinx)*k;
    y_r[k]=gaus_r->GetParameter(1);
    rmsy_r[k]=gaus_r->GetParError(1);
 


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;
    delete gaus_l;
    delete gaus_r;
    delete gaus_t;
  
 

  
 }//chiudo for k


 
  
  for(k=0;k<nbinx;k++){
     
  TH1D* histotemp_l;
  TH1D* histotemp_r;
  TH1D* histotemp_t;
  
  histotemp_l=h3_l->ProjectionY("hc_lprojY",k,k);
  histotemp_r=h3_r->ProjectionY("hc_rprojY",k,k);
  histotemp_t=h3_t->ProjectionY("hc_tprojY",k,k);

  // TF1* gaus_l = new TF1("gaus_l","gaus",histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())-3*histotemp_l->GetRMS(), histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())+3*histotemp_l->GetRMS());
  // TF1* gaus_r = new TF1("gaus_r","gaus",histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())-3*histotemp_r->GetRMS(), histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())+3*histotemp_r->GetRMS());
  // TF1* gaus_t = new TF1("gaus_t","gaus",histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())-3*histotemp_t->GetRMS(), histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())+3*histotemp_t->GetRMS());


        
    // histotemp_l->Fit("gaus_l", "0R");
    // histotemp_r->Fit("gaus_r", "0R");
    // histotemp_t->Fit("gaus_t","0R");
    
    xt[k]=-0.8+(Float_t)(0.8-(-0.8))/(nbinx)*k;
    yt[k]=histotemp_t->GetMean();
    rmsyt[k]=histotemp_t->GetMeanError();
 

    x_l[k]=-0.8+(0.8-(-0.8))/(nbinx)*k;
    y_l[k]=histotemp_l->GetMean();
    rmsy_l[k]=histotemp_l->GetMeanError();
    
    
    x_r[k]=-0.8+(0.8-(-0.8))/(nbinx)*k;
    y_r[k]=histotemp_r->GetMean();
    rmsy_r[k]=histotemp_r->GetMeanError();
 


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;
    //  delete gaus_l;
    // delete gaus_r;
    //delete gaus_t;
  
 

  
 }//chiudo for k
   
 
  

  TCanvas* time_p =new TCanvas("time","Plot time",600,550);
  TLegend* l1= new TLegend(0.1,0.7,0.48,0.9);
  TGraphErrors* gr2_r=new TGraphErrors(nbinx,x_r,y_r,0,rmsy_r);
  TGraphErrors* gr2_l=new TGraphErrors(nbinx,x_l,y_l,0,rmsy_l);
  TGraphErrors* gr2_t=new TGraphErrors(nbinx,xt,yt,0,rmsyt);
 
  
  
  
  TF1* fitgr2_l = new TF1("fitgr2_l","[0]+x*[1]",tymin,tymax);
  TF1* fitgr2_r = new TF1("fitgr2_r","[0]+x*[1]",tymin,tymax);
  TF1* fitgr2_t = new TF1("fitgr2_t","[0]+x*[1]",tymin,tymax);
  


  fitgr2_l->SetLineStyle(2);
  fitgr2_r->SetLineStyle(2);
  fitgr2_t->SetLineStyle(2);
  fitgr2_l->SetLineColor(kRed);
  fitgr2_r->SetLineColor(kBlue);
  fitgr2_t->SetLineColor(kBlack);
  fitgrx_l->SetLineColor(kRed);
  fitgrx_r->SetLineColor(kBlue);
  fitgrx_t->SetLineColor(kBlack);
   
   

  l1->AddEntry(fitgr2_l,"t_{left}-t_{MCP}","L");
  l1->AddEntry(fitgr2_r,"t_{right}-t_{MCP}","L");
  l1->AddEntry(fitgr2_t,"t_{ave}-t_{MCP}","L"); 

  gr2_l->GetYaxis()->SetRangeUser(rymin_l,rymax_l);
  gr2_l->GetXaxis()->SetRangeUser(-0.8,0.8);
  gr2_l->GetYaxis()->SetTitle("t-t_{MCP}(ns)");
  gr2_l->GetXaxis()->SetTitle("t_{left}-t_{right}(ns)");
  gr2_l->SetMarkerColor(kRed);
  gr2_l->SetMarkerStyle(8);
  gr2_l->SetMarkerSize(.5);
  gr2_l->Fit("fitgr2_l","R");
  gr2_l->Draw("AP");
  
  gr2_r->SetMarkerColor(kBlue);
  gr2_r->SetMarkerStyle(8);
  gr2_r->SetMarkerSize(.5);
  gr2_r->Fit("fitgr2_r","R");
  gr2_r->Draw("SAMEP");
    
    
  gStyle->SetOptFit(0001);  
  gr2_t->SetMarkerColor(kBlack);
  gr2_t->SetMarkerStyle(8);
  gr2_t->SetMarkerSize(.5);
  gr2_t->Fit("fitgr2_t");
  gr2_t->Draw("SAMEP");

   l1->Draw();
  
  TCanvas* time_px =new TCanvas("timex","Plot timex",600,550);
  TLegend* l2= new TLegend(0.1,0.7,0.48,0.9);
  TGraphErrors* grx_r=new TGraphErrors(nbinx,x_rx,y_rx,0,rmsy_rx);
  TGraphErrors* grx_l=new TGraphErrors(nbinx,x_lx,y_lx,0,rmsy_lx);
  TGraphErrors* grx_t=new TGraphErrors(nbinx,xtx,ytx,0,rmsytx);


  TF1* fitgrx_l = new TF1("fitgrx_l","[0]+x*[1]",txmin,txmax);
  TF1* fitgrx_r = new TF1("fitgrx_r","[0]+x*[1]",txmin,txmax);
  TF1* fitgrx_t = new TF1("fitgrx_t","[0]+x*[1]",txmin,txmax);

  l2->AddEntry(fitgrx_l,"t_{left}-t_{MCP}","P");
  l2->AddEntry(fitgrx_r,"t_{right}-t_{MCP}","P");
  l2->AddEntry(fitgrx_t,"t_{ave}-t_{MCP}","P");
  grx_t->GetYaxis()->SetTitle("t-t_{MCP}(ns)");
  grx_t->GetXaxis()->SetTitle("X(mm)");
  grx_t->SetMarkerColor(kBlack);
  grx_t->SetMarkerStyle(8);
  grx_t->SetMarkerSize(.5);
  grx_t->Fit("fitgrx_t","R");
  grx_t->Draw("SAMEP");

  grx_r->GetYaxis()->SetTitle("t-t_{MCP}(ns)");
  grx_r->GetXaxis()->SetTitle("X(mm)");
  grx_r->SetMarkerColor(kBlue);
  grx_r->SetMarkerStyle(8);
  grx_r->SetMarkerSize(.5);
  grx_r->Fit("fitgrx_r","R");
  grx_r->Draw("SAMEP");

  
  grx_l->GetYaxis()->SetTitle("t-t_{MCP}(ns)");
  grx_l->GetXaxis()->SetTitle("X(mm)");
  grx_l->SetMarkerColor(kRed);
  grx_l->SetMarkerStyle(8);
  grx_l->SetMarkerSize(.5);
  grx_l->Fit("fitgrx_l","R");
  grx_l->Draw("SAMEP");
  
  l2->Draw();

}
  /*
   TH2D* ht_l= new TH2D("hc_l", "histo hc_l",nbinx,txmin,txmax,nbiny,tymin,tymax);
   TH2D* ht_r= new TH2D("hc_r", "histo hc_r",nbinx,txmin,txmax,nbiny,tymin,tymax);
   TH2D* ht= new TH2D("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin,tymax);


  
   for(k=0;k<digiTree->GetEntries();k++){

     digiTree->GetEntry(k);
    
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	ht_l->Fill(time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0),time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0));
	
        ht_r->Fill(time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0),time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0));

	ht->Fill(time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[3]/max)-hyp_l->GetParameter(0))/2);
		
      }

   }//chiudo for k

 
 
 
   
    for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=ht_l->ProjectionY("ht_lprojY",k,k);
    histotemp_r=ht_r->ProjectionY("ht_rprojY",k,k);
    histotemp_t=ht->ProjectionY("ht_projY",k,k);


    xt[k]=txmin+(Float_t)(txmax-(txmin))/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    rmsyt[k]=histotemp_t->GetMeanError();

    x_l[k]=txmin +(Float_t)(txmax-txmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    rmsy_l[k]=histotemp_l->GetMeanError();

    x_r[k]=txmin+(Float_t)(txmax-txmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    rmsy_r[k]=histotemp_r->GetMeanError();



    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k



  TCanvas* time_pc =new TCanvas("timecorr","Plot timecorr",600,550);
  TLegend* l2= new TLegend(0.1,0.7,0.48,0.9);
  TGraphErrors* g_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* g_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* g_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);
 
  TF1* fitg_l = new TF1("fitg_l","[0]+x*[1]",tymin,tymax);
  TF1* fitg_r = new TF1("fitg_r","[0]+x*[1]",tymin,tymax);
  TF1* fitg_t = new TF1("fitg_t","[0]+x*[1]",tymin,tymax);

  fitg_l->SetLineColor(kRed);
  fitg_r->SetLineColor(kBlue); 
  fitg_t->SetLineColor(kBlack);
 
  
   
  l2->AddEntry(g_l,"t_{left}-t_{MCP}","P");
  l2->AddEntry(g_r,"t_{right}-t_{MCP}","P");
  l2->AddEntry(g_t,"t_{ave}-t_{MCP}","P");
  l2->AddEntry(fitg_l,"t_{left}-t_{MCP}","L");
  l2->AddEntry(fitg_r,"t_{right}-t_{MCP}","L");
  l2->AddEntry(fitg_t,"t_{ave}-t_{MCP}","L"); 
  g_l->GetYaxis()->SetRangeUser(rymin_l,rymax_l);
  g_l->GetXaxis()->SetRangeUser(-0.08,0.6);
  g_l->GetYaxis()->SetTitle("t-t_{MCP}(amp.walk corr)(ns)");
  g_l->GetXaxis()->SetTitle("t_{left}-t_{right}(amp.walk corr)(ns)");
  g_l->SetMarkerColor(kRed);
  g_l->SetMarkerStyle(8);
  g_l->SetMarkerSize(.5);
  g_l->Fit("fitg_l");
  g_l->Draw("AP");
 
  
  g_r->SetMarkerColor(kBlue);
  g_r->SetMarkerStyle(8);
  g_r->SetMarkerSize(.5);
  g_r->Fit("fitg_r");
  g_r->Draw("SAMEP");
    
  gStyle->SetOptFit(0001);
  g_t->SetMarkerColor(kBlack);
  g_t->SetMarkerStyle(8);
  g_t->SetMarkerSize(.5);
  cout << "*********************************" <<endl;
  cout << "*********************************" <<endl;
  g_t->Draw("SAMEP");
  l2->Draw();


  
  
  TH2D* htdiff= new TH2D("hc_tdiff", "histo hc_tdiff",nbinx,txmin,txmax,nbiny,tymin,tymax);
  g_t->Fit("fitg_t");

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);
    
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {

	htdiff->Fill(time[1+LEDi]-time[2+LEDi]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[3]/max)-hyp_l->GetParameter(0))/2-fitg_t->Eval(time[1+LEDi]-time[2+LEDi]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)))+fitg_t->GetParameter(0)); //GIUSTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      }

  }//chiudo for k 

  for(k=0;k<nbinx;k++){
    
    TH1D* histotemp_t;
    histotemp_t=htdiff->ProjectionY("ht_projY",k,k);
    
    
    xt[k]=txmin+(Float_t)(txmax-(txmin))/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    rmsyt[k]=histotemp_t->GetMeanError();
    
    delete histotemp_t;
    
    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k
  
  TGraphErrors* g_tdiff=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);
  TF1* fitc_tdiff = new TF1("fitc_tdiff","[0]+x*[1]",tymin,tymax);
  
  fitc_tdiff->SetLineColor(kGreen);
  
  
  g_tdiff->SetMarkerColor(kGreen);
  g_tdiff->SetMarkerStyle(8);
  g_tdiff->SetMarkerSize(.5);
  cout << "__________________________" << endl;
  g_tdiff->Fit("fitc_tdiff");
  cout << "__________________________" << endl;
  g_tdiff->Draw("SAMEP");
  l2->Draw();
  
  wf_c->cd(4);
  ht->Draw("COLZ");
  g_t->Fit("fitg_t","0L");
  g_t->SetMarkerStyle(8);
  g_t->SetMarkerSize(.5);
  g_t->Draw("P");
  fitg_t->Draw("same");

  

  }
*/
