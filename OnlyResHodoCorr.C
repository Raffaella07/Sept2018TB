#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "stdio.h"
#include "string.h"


Float_t GetHodoPosition( Int_t* nFibres, Float_t* var ) {

  Int_t nFibresMax = 10;

  bool showering0 = nFibres[0]>=nFibresMax;
  bool showering1 = nFibres[1]>=nFibresMax;

  bool bad0 = nFibres[0]==0 || var[0]<-900.;
  bool bad1 = nFibres[1]==0 || var[1]<-900.;

  Float_t pos = -999.;

  if( showering0 || showering1 ) pos = -999.;  // showering (at least one plane)

  else if( bad0 && bad1 ) pos = -999.;  // empty (both planes have 0 fibres)

  else {

    if( bad0 ) pos = var[1];

    else if( bad1 ) pos = var[0];

    else pos = 0.5*(var[0]+var[1]);

  } // else

  return pos;

}

TF1* fitGaus( TH1D* histo, float nSigma, bool addFunc, string type ) {

  float mean_histo = histo->GetMean();
  float rms_histo  = histo->GetRMS();
  
  TCanvas* iterCanvas = new TCanvas("iterFitCanvas","",1200,800);
  
  TF1* f1_gaus = new TF1( Form("gaus_%s", histo->GetName()), "gaus", mean_histo-rms_histo, mean_histo+rms_histo );
  f1_gaus->SetLineColor( 46 );
  histo->Draw("HISTO");
  histo->Fit( f1_gaus->GetName(), "RQ" );
  f1_gaus->Draw("SAME");

  cout <<"Fit numero "<< "0"<< "#sigma=" <<f1_gaus->GetParameter(2)<<"->"<<sqrt(f1_gaus->GetParameter(2)*f1_gaus->GetParameter(2)-0.018*0.018) <<endl; 

  iterCanvas->SaveAs(("IterfitControl/"+type+"/FitZero.pdf").c_str());
  
  delete iterCanvas;
  
  float xMin_fit = f1_gaus->GetParameter(1) - nSigma*f1_gaus->GetParameter(2);
  float xMax_fit = f1_gaus->GetParameter(1) + nSigma*f1_gaus->GetParameter(2);
  
  f1_gaus->SetRange( xMin_fit, xMax_fit );
  
  
  int n_iter = 5;
  TLegend* legenda = new TLegend();  
  
  for( int i=0; i<n_iter; ++i ) { // iterative fit
    TCanvas* iterCanvas = new TCanvas("iterFitCanvas","",1200,800);
    if( i==n_iter-1 && addFunc )
      histo->Fit( f1_gaus->GetName(), "RQ+" );
    else {
      histo->Draw("HISTO");
      histo->Fit( f1_gaus->GetName(), "RQ" );
      f1_gaus->Draw("SAME");
      xMin_fit = f1_gaus->GetParameter(1) - nSigma*f1_gaus->GetParameter(2);
      xMax_fit = f1_gaus->GetParameter(1) + nSigma*f1_gaus->GetParameter(2);
      f1_gaus->SetRange( xMin_fit, xMax_fit );
    }
    if(i==n_iter-1){
      legenda->AddEntry(f1_gaus,("#sigma_t="+to_string(f1_gaus->GetParameter(2))).c_str());
      legenda->Draw();
    }
    
    iterCanvas->SaveAs(("IterfitControl/"+type+"/Fit"+to_string(i)+type+".pdf").c_str());
    cout <<"Fit numero "<< i<< "#sigma=" <<f1_gaus->GetParameter(2)<<"->"<<sqrt(f1_gaus->GetParameter(2)*f1_gaus->GetParameter(2)-0.018*0.018) <<endl; 
    delete iterCanvas;
    
  } // for iter
  

  return f1_gaus;
  
}


float GetSigmaEff( TH1D* histo, string dir) {

  float percentIntegral = 0.683;
  float integral = histo->Integral();

  // first fit to find mode:
  TF1* f1_gaus = fitGaus( histo, 1.5, false, "a");

  float mode = f1_gaus->GetParameter(1);
  int maxBin = histo->FindBin( mode );

  int nBins = histo->GetNbinsX();
  float xMin = histo->GetXaxis()->GetXmin();
  float xMax = histo->GetXaxis()->GetXmax();

  TH1D* newHisto = new TH1D( Form("newHisto_%s", histo->GetName()), "", nBins, xMin, xMax);
  newHisto->SetBinContent( maxBin, histo->GetBinContent(maxBin) );
  newHisto->SetBinError( maxBin, histo->GetBinError(maxBin) );

  Int_t iBin = maxBin;
  Int_t delta_iBin = 1;
  Int_t sign  = 1;

  float width = histo->GetBinWidth( maxBin );
  
  TCanvas* ControlloSigmaEff;
  int i=0;
  
  while( newHisto->Integral() < percentIntegral*integral ) {
    
    ControlloSigmaEff= new TCanvas("ControlloSigmaEff","",1200,800);
    
    iBin += sign*delta_iBin; 
    
    newHisto->SetBinContent( iBin, histo->GetBinContent(iBin) );
    newHisto->SetBinError( iBin, histo->GetBinError(iBin) );
    
    width += histo->GetBinWidth( iBin );
    
    delta_iBin += 1;
    sign *= -1;
    newHisto->Draw("HISTO");
    ControlloSigmaEff->SaveAs(("IterfitControl/"+dir+"/"+to_string(i)+".pdf").c_str());
    i++;
    delete ControlloSigmaEff;
    
  }

  return width/2.;

}


void AWtdiff(const char * filename){

  gSystem->Exec("rm -r IterfitControl");
  gSystem->Exec("mkdir IterfitControl");
  gSystem->Exec("mkdir IterfitControl/t");
  gSystem->Exec("mkdir IterfitControl/tdiff");
  gSystem->Exec("mkdir IterfitControl/Efft");
  gSystem->Exec("mkdir IterfitControl/Efftdiff");
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptFit(1111);

  
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

  const Int_t  nbinx=150,nbiny=800;

  Int_t i;
  Double_t sigma[50],erry[50],cut[50],errx[50];
  
  txmin=-0.8;
  txmax=0.8;
  
  
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

  Int_t nFibresOnX[2],nFibresOnY[2];

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
  hodoTree->SetBranchAddress( "nFibresOnX", nFibresOnX );
  hodoTree->SetBranchAddress( "nFibresOnY", nFibresOnY );


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
    
    if(time[NINO1+LEDi]-time[0+CFD]<15 && time[NINO1+LEDi]-time[0+CFD]>0) {
      counter1++;
      Times1[counter1]=time[1+LEDi]-time[0+CFD];
    }
    
    if(time[NINO2+LEDi]-time[0+CFD]<15 && time[NINO2+LEDi]-time[0+CFD]>0) {
      counter2++;
      Times2[counter2]=time[NINO2+LEDi]-time[0+CFD];
    }
    
    if((time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[0+CFD]<15 && (time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[0+CFD]>-1) {
      counter3++;
      Times3[counter3]=(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[0+CFD];
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

  rymin_l=mean1-0.8*rms1;
  rymax_l=mean1+0.5*rms1;
  rymin_r=mean2-0.8*rms2;
  rymax_r=mean2+0.5*rms2;
    
  

  tymin=mean3-0.7*rms3;
  tymax=mean3+0.3*rms3;

  /*rymin_l=0;
  rymax_l=10;
  rymin_r=0;
  rymax_r=10;
    
  
  
  tymin=0;
  tymax= 5;*/
  
  
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

  TCanvas* lcheck = new TCanvas("landaucheck","landaucheck plot", 1200, 550);
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

  lcheck->SaveAs("IterfitControl/landau.pdf");

  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2F* h2_time= new TH2F("h2_time", "histo h2_t",nbinx,-0.8,0.8,nbiny,tymin,tymax);

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
      }
    }
  }
  
  
	   
  TH1D* DistribUncorr = h2_time->ProjectionY("RMSUncorr",0,nbinx);
  DistribUncorr->GetXaxis()->SetTitle("t_{ave}-t_{MCP}");
  TCanvas* largeg = new TCanvas("large_g","plot large_g", 550, 600);
  DistribUncorr->Fit("gaus");
  DistribUncorr->Draw("same");
  largeg->SaveAs("IterfitControl/ResUncorr.pdf");
  
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
    
    histotemp_l->Fit("gaus_l", "Q0");
    histotemp_r->Fit("gaus_r", "Q0");
    histotemp_t->Fit("gaus_t","Q0");

    //xt[k]=-0.8+(Float_t)(0.8-(-0.8))/(nbinx/1.5)*k;
    xt[k]=txmin+(Float_t)(txmax-(txmin))/(nbinx)*k;
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
    //delete gaus_l;
    //delete gaus_r;
    //delete gaus_t;
    
    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k
  
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1200,600);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx,xt,yt,0,rmsyt);
  
  //TF1* hyp_r = new TF1("hyp_r","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*1/(x+[5])",0.15,0.80);
  //TF1* hyp_l = new TF1("hyp_l","[0]+[1]*x+[2]*x**2+[3]*x**3",0.03,0.15);
  
  TF1* hyp_r = new TF1("hyp_r","pol5",0.15,0.80);
  TF1* hyp_l = new TF1("hyp_l","pol5",0.03,0.15);
  
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
  h2_time->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
  h2_time->GetXaxis()->SetTitle("t_{left}-t_{right}[ns]");
  h2_time->Draw("COLZ");
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
  TH2D* hc_time= new TH2D("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin,tymax);
  TH2D* hc_tdiff= new TH2D("hc_tdiff", "histo hc_tdiff",nbinx,txmin,txmax,nbiny,tymin,tymax);
  

  
  for(k=0;k<digiTree->GetEntries();k++){
    
    digiTree->GetEntry(k);
    hodoTree->GetEntry(k);
    
    if  (amp_max[PTK1]/max > 0.1 && amp_max[PTK1]/max < 0.55 && GetHodoPosition(nFibresOnX,X)>-11 && GetHodoPosition(nFibresOnY,Y)>-2 && GetHodoPosition(nFibresOnY,Y)<5){
      if ((0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1))) ) 	hc_l->Fill(amp_max[AMP1]/max,time[NINO1+LEDi]-time[0+CFD]-hyp_l->Eval(amp_max[AMP1]/max)+hyp_l->GetParameter(0));
      
      if (((0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1)))) ) hc_r->Fill(amp_max[AMP2]/max,time[NINO2+LEDi]-time[0+CFD]-hyp_r->Eval(amp_max[AMP2]/max)+hyp_r->GetParameter(0)); 
      
      if ((0.8*(fit_l->GetParameter(1)) < (amp_max[AMP1]/max) && (amp_max[AMP1]/max) < (3*fit_l->GetParameter(1))) || ((0.8*(fit_r->GetParameter(1)) < (amp_max[AMP2]/max) && (amp_max[AMP2]/max) < (3*fit_r->GetParameter(1)))) )
	{
	  hc_time->Fill(time[NINO1+LEDi]-time[NINO2+LEDi],(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[0+CFD]-(hyp_r->Eval(amp_max[AMP2]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[AMP1]/max)-hyp_l->GetParameter(0))/2);	
	}//chiudo if
    }
  }//chiudo for k
  

  TH1D* DistribAWCorr = hc_time->ProjectionY("DistribAWCorr",0,nbinx);
  DistribAWCorr->GetXaxis()->SetTitle("t_{ave}-t_{MCP}");
  TCanvas* corrg = new TCanvas("corr_g","plot corr_g", 550, 600);
  DistribAWCorr->Fit("gaus");
  DistribAWCorr->Draw("same");
  
  corrg->SaveAs("IterfitControl/AWCorr.pdf");

  for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;
     
    histotemp_l=hc_l->ProjectionY("hc_lprojY",k,k);
    histotemp_r=hc_r->ProjectionY("hc_rprojY",k,k);
    histotemp_t=hc_time->ProjectionY("hc_tprojY",k,k);
    
    TF1* gaus_l = new TF1("gaus_l","gaus",histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())-3*histotemp_l->GetRMS(), histotemp_l->GetBinCenter(histotemp_l->GetMaximumBin())+3*histotemp_l->GetRMS());
    TF1* gaus_r = new TF1("gaus_r","gaus",histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())-3*histotemp_r->GetRMS(), histotemp_r->GetBinCenter(histotemp_r->GetMaximumBin())+3*histotemp_r->GetRMS());
    TF1* gaus_t = new TF1("gaus_t","gaus",histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())-3*histotemp_t->GetRMS(), histotemp_t->GetBinCenter(histotemp_t->GetMaximumBin())+3*histotemp_t->GetRMS());
    
    histotemp_l->Fit("gaus_l", "Q0");
    histotemp_r->Fit("gaus_r", "Q0");
    histotemp_t->Fit("gaus_t","Q0");
    
    xt[k]=txmin+(Float_t)(txmax-(txmin))/(nbinx)*k;
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
  TGraphErrors* graph_tc = new TGraphErrors(nbinx,xt,yt,0,rmsyt);
  
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
  
  hc_time->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
  hc_time->GetXaxis()->SetTitle("t_left-t_right[mm]");
  hc_time->Draw("COLZ");
  graph_tc->SetMarkerStyle(8);
  graph_tc->SetMarkerSize(.5);
  graph_tc->Draw("P");
  
  wf_c->SaveAs("IterfitControl/6plot.pdf");
  
  float sigmaEffT,sigmaEffTdiff;

    
  TH1D* histoAWcorr = hc_time->ProjectionY("histo projection",0,nbinx);
  hc_time->GetXaxis()->SetRangeUser(hc_time->GetMean()-3*hc_time->GetRMS(),hc_time->GetMean()+3*hc_time->GetRMS());
  TF1* iterFitCt = fitGaus(histoAWcorr,1.7,false,"t");
  sigmaEffT=GetSigmaEff(histoAWcorr ,"Efft");
  
  //Correzione posizione con Hodoscopio

  TH2D* TimeAveVsX= new TH2D("TimeAveVsX", "TimeAve vs HodoPosition",nbinx,-20,20,nbiny*2,tymin,tymax);
  TH2D* AmpLVsX= new TH2D("AmpLVsX", "AmpL vs HodoPosition",nbinx,-20,20,nbiny*2,0,1);
  TH2D* AmpRVsX= new TH2D("AmpRVsX", "AmpR vs HodoPosition",nbinx,-20,20,nbiny*2,0,1);
  
  Int_t counter=0;
  Float_t Pos;
  cout<<"HERE"<<endl;
  for(k=0;k<digiTree->GetEntries();k++){
    hodoTree->GetEntry(k);
    Pos=GetHodoPosition(nFibresOnX,X);
    digiTree->GetEntry(k);
    
    if (  0.8*fit_l->GetParameter(1) < amp_max[AMP1]/max && amp_max[AMP1]/max < 3*fit_l->GetParameter(1) ) {
      AmpLVsX->Fill( Pos, amp_max[AMP1]/max);
    }
    
    if (  0.8*fit_r->GetParameter(1) < amp_max[AMP2]/max && amp_max[AMP2]/max < 3*fit_r->GetParameter(1) ) {
      AmpRVsX->Fill( Pos,amp_max[AMP2]/max);
    }
    
    if (  0.8*fit_r->GetParameter(1) < amp_max[AMP2]/max && amp_max[AMP2]/max < 3*fit_r->GetParameter(1) ){ 
      TimeAveVsX->Fill(Pos,(time[NINO1+LEDi]+time[NINO2+LEDi])/2-time[0+CFD]-(hyp_r->Eval(amp_max[AMP2]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[AMP1]/max)-hyp_l->GetParameter(0))/2);
      counter++;
    }//chiudo if
    
  }//chiudo for
  
  TCanvas* CanvTimeAveX = new TCanvas("Canvas T_aveCorrVsX","",1200,800);
  TimeAveVsX->Draw("COLZ");
  TCanvas* CanvAmpRVsX = new TCanvas("Canvas T_aveCorrVsX","",1200,800);
  AmpRVsX->Draw("COLZ");
  TCanvas* CanvAmpLVsX = new TCanvas("Canvas T_aveCorrVsX","",1200,800);
  AmpLVsX->Draw("COLZ");
  
  CanvTimeAveX->SaveAs("IterfitControl/t_aveCorrVsX.pdf");
  CanvAmpRVsX->SaveAs("IterfitControl/AmpRVsX.pdf");
  CanvAmpLVsX->SaveAs("IterfitControl/AmpLVsX.pdf");
  
  /*
    TF1* iterFitCtdiff = fitGaus(histo_ctdiff,1.7,false,"tdiff");
    sigmaEffTdiff=GetSigmaEff(histo_ctdiff ,"Efftdiff");
  */

  //è una prova, non è tdiff
  TH1D* histoProjFromX = TimeAveVsX->ProjectionY("histo projection",0,nbinx);
  histoProjFromX->GetXaxis()->SetRangeUser(histoProjFromX->GetMean()-3*histoProjFromX->GetRMS(),histoProjFromX->GetMean()+3*histoProjFromX->GetRMS());
  TF1* iterFittdiff = fitGaus(histoProjFromX,1.7,false,"tdiff");
  sigmaEffT=GetSigmaEff(histoProjFromX ,"Efftdiff");
  //è una prova, non è tdiff


  cout<< "sigmaEffT=" << sigmaEffT << "     " << "sigmaEffTdiff=" << sigmaEffT<< endl;
  cout<<counter<<endl;
}
