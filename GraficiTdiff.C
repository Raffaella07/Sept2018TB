//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h> 
#include <TLine.h>
#include <TAxis.h>


void GetH2MaximumBins(TH2D* histo, Int_t* XMaxBin, Int_t* YMaxBin){
  
  TH1D* proj;
  Int_t k,nxbin,nybin;
  Double_t maximum=0;
  
  nxbin=histo->GetXaxis()->GetNbins();
  nybin=histo->GetYaxis()->GetNbins();
   
  for(k=0;k<nxbin;k++){
    
    proj=histo->ProjectionY("Yprojection",k,k);
    proj->FindBin(proj->GetMaximum());
    
    if(maximum<proj->GetMaximum()){
      cout<<k<<endl;
      maximum=proj->GetMaximum();
      *XMaxBin = k;
      *YMaxBin = proj->GetMaximumBin();
      
    }
  
  }
  cout<<"infunction_______________________xbin="<<*XMaxBin<<"/"<<nxbin<<"          ybin"<<*YMaxBin<<"/"<<nybin<<endl;
}

void GraficiTdiff(const char * filename){

  
  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");



  Float_t amp_max[54], time[54];
  int k,j,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  bool debug=false;
  bool blind=true;
  Double_t max=0;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi;
  rxmin=0;
  rxmax=0.5;

  Int_t nentries=digiTree->GetEntries();
  Float_t Times1[nentries],Times2[nentries],Times3[nentries];

  const Int_t  nbinx=500,nbiny=450;
  
  int i;
  Double_t sigma[50],erry[50],cut[50],errx[50];
  
  txmin=-0.4;
  txmax=0.6;
  
  
  
  Double_t x_r[nbinx],y_r[nbinx], x_l[nbinx],y_l[nbinx],rmsy_l[nbinx],rmsy_r[nbinx];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];


  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1F *mcp_amp =new TH1F("mcp_amp","histomcp_ampl",nbinx,0.0,1);


  TF1 *fit_r = new TF1("f_r","landau",0.14,1);
  TF1 *fit_l = new TF1("f_l","landau",0.14,1);


  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
  digiTree->SetBranchAddress("LED300",&LED300);
  digiTree->SetBranchAddress("LED100",&LED100);
  digiTree->SetBranchAddress("LED50",&LED50);
  digiTree->SetBranchAddress("LED30",&LED30);

  digiTree->GetEntry(3);
  LEDi=LED300;
  
for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    //cout<<"HERE"<<endl;
    if(time[1+LEDi]-time[0]<50 && time[1+LEDi]-time[0]>0) {Times1[k]=time[1+LEDi]-time[0];}
    else{Times1[k]=Times2[k-1];}
    if(time[2+LEDi]-time[0]<50 && time[2+LEDi]-time[0]>0) {Times2[k]=time[2+LEDi]-time[0];}
    else{Times2[k]=Times2[k-1];}  
    if((time[1+LEDi]+time[2+LEDi])/2-time[0]<50 && (time[1+LEDi]+time[2+LEDi])/2-time[0]>-10) {Times3[k]=(time[1+LEDi]+time[2+LEDi])/2-time[0];}
    else{Times3[k]=Times3[k-1];}
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
  
  



  max=4096;
  
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);

    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
    mcp_amp->Fill(amp_max[0]/max);

  }/*chiudo for */

  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  //cout << tmax <<endl;
  cout<< max << endl;

  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");


  TH2D* h2_l= new TH2D("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2D* h2_r= new TH2D("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2D* h2_t= new TH2D("h2_t", "histo h2_t",nbinx,txmin,txmax,nbiny,tymin,tymax);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);


    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	h2_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]);
	if (amp_max[4]/max < 0.35)h2_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]);
	h2_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]);

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

    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1800,1100);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);

  TF1* hyp_r = new TF1("hyp_r","[0]+[2]*log(x+[1])",0.135,0.35);
  TF1* hyp_l = new TF1("hyp_l","[0]+[2]*log(x+[1])",0.11,0.35);

  TF1* hyp_t = new TF1("hyp_t","[1]*x**2+[2]*x+[0]",-0.1,0.65);
  
  gStyle->SetOptStat("");


 /* SetParameters*/
  hyp_l->SetParameter(0, 8.51);
  hyp_l->SetParameter(1, 5);
  hyp_l->SetParameter(2, 1.2);
  /* hyp_l->SetParameter(3, -2.43e-2);
  */
  hyp_r->SetParameter(0, 7);
  hyp_r->SetParameter(1, -9e-2);
  hyp_r->SetParameter(2, -1e-1);

 
 wf_c->Divide(3,2);

  wf_c->cd(1);

  h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");

   h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_l","0RL");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("SAMEP");
  hyp_l->Draw("same");


  wf_c->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","R0L");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  graph_r->Draw("SAMEP");
  hyp_r->Draw("same");
  
  wf_c->cd(3);
  h2_t->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
  h2_t->GetXaxis()->SetTitle("t_left-t_right [ns]");
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","RL0M");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  graph_t->Draw("SAMEP");
  hyp_t->Draw("same");

  rymin_lc=rymin_l-hyp_l->Eval(0.25)+hyp_l->GetParameter(0);
  rymax_lc=rymax_l-hyp_l->Eval(0.25)+hyp_l->GetParameter(0);
  rymin_rc=rymin_r-hyp_r->Eval(0.25)+hyp_r->GetParameter(0);
  rymax_rc=rymax_r-hyp_r->Eval(0.25)+hyp_r->GetParameter(0);
  tymin_c=tymin;
  tymax_c=tymax;

  
  TH2D* hc_l= new TH2D("hc_l", "histo hc_l",nbinx,txmin,txmax,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_r= new TH2D("hc_r", "histo hc_r",nbinx,txmin,txmax,nbiny,rymin_rc,rymax_rc);
  TH2D* hc_t= new TH2D("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  TH2D* hc_tdiff= new TH2D("hc_tdiff", "histo hc_tdiff",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  
  
  
   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	hc_l->Fill((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))),time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0));
        hc_r->Fill((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))),time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0));
	hc_t->Fill((time[1+LEDi]-(time[2+LEDi])),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0))/2);

	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
      }

  }//chiudo for k

    for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=hc_l->ProjectionY("hc_lprojY",k,k);
    histotemp_r=hc_r->ProjectionY("hc_rprojY",k,k);
    histotemp_t=hc_t->ProjectionY("hc_tprojY",k,k);
   
    xt[k]=txmin+(txmax-txmin)/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    RMS[2][k]= histotemp_t->GetMeanError();
    
    x_l[k]=txmin+(txmax-txmin)/nbinx*k;
    y_l[k]=histotemp_l->GetMean();
    RMS[0][k]= histotemp_l->GetMeanError();
    
    x_r[k]=txmin+(txmax-txmin)/nbinx*k;
    y_r[k]=histotemp_r->GetMean();
    RMS[1][k]= histotemp_r->GetMeanError();


    delete histotemp_l;
    delete histotemp_r;
    delete histotemp_t;

    
  }//chiudo for k


   TGraphErrors* graph_lc = new TGraphErrors(nbinx-1,x_l,y_l,0,RMS[0]);
   TGraphErrors* graph_rc = new TGraphErrors(nbinx-1,x_r,y_r,0,RMS[1]);
   TGraphErrors* graph_tc = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);
   TF1* fit_tdiff = new TF1("fit_tdiff","[0]+[1]*x",-0.1,0.65);
   TF1* fit_lc = new TF1("fit_lc","[0]+[1]*x",-0.1,0.65);
   TF1* fit_rc = new TF1("fit_rc","[0]+[1]*x",-0.1,0.65);
  
   graph_tc->Fit("fit_tdiff","R");
   graph_rc->Fit("fit_rc","R");
   graph_lc->Fit("fit_lc","R");


   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
		 hc_tdiff->Fill(time[1+LEDi]-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[3]/max)-hyp_l->GetParameter(0))/2-fit_tdiff->Eval(time[1+LEDi]-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0))+fit_tdiff->GetParameter(0));
      }

  }//chiudo for k
   

   for(k=0;k<nbinx;k++){ 
     TH1D* histotemp_t;
     
     histotemp_t=hc_tdiff->ProjectionY("hc_tprojY",k,k);
     
     yt[k]=histotemp_t->GetMean();
     RMS[2][k]= histotemp_t->GetMeanError();
     
     delete histotemp_t;
   }//chiudo for k
   
    TGraphErrors* graph_tcdiff = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);

    
    wf_c->cd(4);
    hc_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
    hc_l->GetXaxis()->SetTitle("t_left-t_right [mV]");
    hc_l->Draw("COLZ");
    graph_lc->SetMarkerStyle(8);
    graph_lc->SetMarkerSize(.5);
    graph_lc->Draw("P");

  

    wf_c->cd(5);
    
    hc_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
    hc_r->GetXaxis()->SetTitle("t_left-t_right [mV]");
    hc_r->Draw("COLZ");
    graph_rc->SetMarkerStyle(8);
    graph_rc->SetMarkerSize(.5);
    graph_rc->Draw("P");
    
  

   wf_c->cd(6);

   hc_tdiff->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
   hc_tdiff->GetXaxis()->SetTitle("t_left-t_right [ns]");
   hc_tdiff->Draw("COLZ");
   graph_tcdiff->SetMarkerStyle(8);
   graph_tcdiff->SetMarkerSize(.5);
   graph_tcdiff->Draw("P");
   
   Int_t MaxBinX,MaxBinY;
   
   GetH2MaximumBins(hc_tdiff, &MaxBinX, &MaxBinY);
  
   TLine* lineax= new TLine(hc_tdiff->GetXaxis()->GetBinCenter(MaxBinX),0,hc_tdiff->GetXaxis()->GetBinCenter(MaxBinX),10);
   TLine* lineay= new TLine(-1,hc_tdiff->GetYaxis()->GetBinCenter(MaxBinY),10,hc_tdiff->GetYaxis()->GetBinCenter(MaxBinY));
   
   cout<<"____________________________________________________________________________maxX="<<hc_tdiff->GetXaxis()->GetBinCenter(MaxBinX)<<"    maxY="<<hc_tdiff->GetYaxis()->GetBinCenter(MaxBinY)<<endl; 
   lineax->Draw("SAME");
   lineay->Draw("SAME");

  

   
   hc_r->Reset();
   hc_l->Reset();
   
   for(k=0;k<digiTree->GetEntries();k++){
     digiTree->GetEntry(k);

     if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
       {
	 hc_l->Fill((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))),time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-fit_lc->Eval((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))))+fit_lc->GetParameter(0));
	 
	 hc_r->Fill((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))),time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0)-fit_rc->Eval((time[1+LEDi]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)-(time[2+LEDi]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0))))+fit_rc->GetParameter(0)); 	 
       }// chiudo if
     
   }//chiudo for k
   
   /*
   plottini->cd(1);
   hc_l->Draw("COLZ");
   plottini->cd(2);
   hc_r->Draw("COLZ");
   */
   
   TH1D* histo_cl;
   TH1D* histo_cr;
   TH1D* histo_ct;
   TH1D* histo_ctdiff;
   
   TF1* gaus_cl = new TF1("gaus_cl","gaus",-2.5,-0.7);
   TF1* gaus_cr = new TF1("gaus_cr","gaus",-2.5,-0.7);
   TF1* gaus_ct = new TF1("gaus_ct","gaus",-2.5,-0.7);
   TF1* gaus_ctdiff = new TF1("gaus_ctdiff","gaus",6,8);

   histo_cl = hc_l->ProjectionY("histo_cl",0,nbinx);
   histo_cr = hc_r->ProjectionY("histo_cr",0,nbinx);
   histo_ct = hc_t->ProjectionY("histo_ct",0,nbinx);
   histo_ctdiff = hc_tdiff->ProjectionY("histo_ctdiff",0,nbinx);


   histo_ct->SetLineColor(kBlack);
   histo_cl->SetLineColor(kBlue);
   histo_cr->SetLineColor(kRed);
   gaus_ct->SetLineColor(kBlack);


   

   TCanvas* tdiff = new TCanvas("tdiff","plot_tdiff",600,550);
   TLegend* l2=new TLegend(0.1,0.7,0.48,0.9);
   l2->SetHeader("time stamps");
   l2->AddEntry(histo_ct,"t_ave-t_MCP");
   l2->AddEntry(histo_ctdiff,"t_ave-t_MCP(tdiff corr)");
   gaus_ctdiff->SetLineColor(kGreen);
   histo_ct->Fit("gaus_ct");
   histo_ctdiff->Fit("gaus_ctdiff");
   histo_ctdiff->SetLineColor(kGreen);
   histo_ct->Draw();
   histo_ctdiff->Draw("same");
   l2->Draw();
 
   
   cut[0]=0.4;
   cut[1]=0.2;
   cut[2]=0.1;
   cut[3]=0.05;
   cut[4]=0.02;
   cut[5]=0.01;
   
   TCanvas* canvino = new TCanvas("canvino","",1200,600);
   TF1* fittino[3][6];
   canvino->Clear();
   canvino->Divide(3,2);
   
   TH1D* istogrammi[3][6];

 
   Int_t DeltaBin;
   for (i=0;i<6;i++){
     
     
     fittino[0][i] = new TF1("argaerg"+i ,"gaus",6,8);
     fittino[1][i] = new TF1("xzcbxcn"+i ,"gaus",6,8);
     fittino[2][i] = new TF1("jklnjbm"+i ,"gaus",6,8);
     
     Double_t XCenter=hc_tdiff->GetXaxis()->GetBinCenter(MaxBinX);
     DeltaBin= (int)(hc_tdiff->GetXaxis()->FindBin(XCenter+cut[i])-hc_tdiff->GetXaxis()->FindBin(XCenter-cut[i]))/2;
     
  
     istogrammi[0][i]=hc_tdiff->ProjectionY("ghijklxz"+i,MaxBinX-DeltaBin,MaxBinX+DeltaBin);
     istogrammi[1][i]=hc_l->ProjectionY("abcdefuv"+i, hc_l->GetXaxis()->FindBin(TMath::Max(-cut[i],(Double_t)txmin)), hc_l->GetXaxis()->FindBin(cut[i]));
     istogrammi[2][i]=hc_r->ProjectionY("mnopqrst"+i, hc_r->GetXaxis()->FindBin(TMath::Max(-cut[i],(Double_t)txmin)), hc_r->GetXaxis()->FindBin(cut[i]));

     
     istogrammi[0][i]->SetLineColor(kBlack);
     istogrammi[1][i]->SetLineColor(kBlue);
     istogrammi[2][i]->SetLineColor(kRed);
     istogrammi[0][i]->GetXaxis()->SetTitle("t_{ave}-t_{MCP}");
     istogrammi[0][i]->GetYaxis()->SetTitle("counts");
     istogrammi[1][i]->GetXaxis()->SetTitle("t_{ave}-t_{MCP}");
     istogrammi[1][i]->GetYaxis()->SetTitle("counts");
     istogrammi[2][i]->GetXaxis()->SetTitle("t_{ave}-t_{MCP}");
     istogrammi[2][i]->GetYaxis()->SetTitle("counts");
     
     cout<<"__________________"<<MaxBinX-DeltaBin <<"     "<<MaxBinX+DeltaBin<<"__________"<<hc_tdiff->GetXaxis()->GetBinCenter(MaxBinX-DeltaBin) <<"        "<<hc_tdiff->GetXaxis()->GetBinCenter(MaxBinX+DeltaBin) <<endl;

     //canvino->cd(i+1);
     istogrammi[0][i]->Fit("argaerg"+i,"0");
     cout<< "Gaussiana centrale" <<endl;
     istogrammi[1][i]->Fit("xzcbxcn"+i,"0");
     istogrammi[2][i]->Fit("jklnjbm"+i,"0");
     
       
     
   }

   TLegend* legenda[6];
   TString LegSigma;
   TString TTempo;

   for(i=0;i<6;i++){ 
     LegSigma="";
     LegSigma.Append("sigma=");
     TTempo="";
     TTempo.Append("|dt|<");
     TTempo.Append(to_string(cut[i]));
     TTempo.Resize(TTempo.Sizeof()-4);
     TTempo.Append("ns");
     canvino->cd(i+1);
     
     //istogrammi[1][i]->Draw();
     //     istogrammi[2][i]->Draw("SAME");
     istogrammi[0][i]->Draw("SAME");
     
     fittino[0][i]->Draw("SAME");
     //fittino[1][i]->Draw("SAME");
     //fittino[2][i]->Draw("SAME");
     LegSigma.Append(to_string(sqrt(fittino[0][i]->GetParameter(2)*fittino[0][i]->GetParameter(2)-0.015*0.015)));
     legenda[i] = new TLegend();
     legenda[i]->AddEntry(fittino[0][i],LegSigma);
     legenda[i]->AddEntry(istogrammi[0][i],TTempo);
     legenda[i]->Draw();
   }
   
}
