//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <TGraphErrors.h>
#include <TStyle.h>

void AmplitudeStudy(const char * filename){


  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");



  Float_t amp_max[54], time[54];
  int k,maxbin_l,maxbin_r,maxbin_t;
  Float_t max;
  const Int_t  nbinx=170,nbiny=240;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  rxmin=0;
  rxmax=0.7;
  txmin=-1.2;
  txmax=1.2;

  Double_t x_r[nbinx],y_r[nbinx], x_l[nbinx],y_l[nbinx],rmsy_l[nbinx],rmsy_r[nbinx];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];
  Int_t i;


  TH1D *hr_amp =new TH1D("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1D *hl_amp =new TH1D("hl_amp","histos_ampl",nbinx,0.0,1);

  TF1 *fit_r = new TF1("f_r","landau",0.10,1);
  TF1 *fit_l = new TF1("f_l","landau",0.10,1);



  Int_t nentries=digiTree->GetEntries(), counter1=0,counter2=0, counter3=0;
  Float_t Times1[nentries],Times2[nentries],Times3[nentries];
  Int_t LEDi=24; 
  for(k=0;k<digiTree->GetEntries();k++){
    Times1[k]=0;
    Times2[k]=0;
    Times3[k]=0;
  }  
  
  
  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
 
  digiTree->GetEntry(3);
 
  for(k=0;k<digiTree->GetEntries();k++){
    digiTree->GetEntry(k);
    
    if(time[1+LEDi]-time[0]<40 && time[1+LEDi]-time[0]>0) {
      counter1++;
      Times1[counter1]=time[1+LEDi]-time[0];
    }
    
    if(time[2+LEDi]-time[0]<40 && time[2+LEDi]-time[0]>0) {
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

  rymin_l=-20;//mean1-1.2*rms1;
  rymax_l=12;//mean1+0.8*rms1;
  rymin_r=-20;//mean2-1.2*rms2;
  rymax_r=12;//mean2+0.8*rms2;
  
  tymin=mean3-1.2*rms3;
  tymax=mean3+0.8*rms3;


  max=4096;
  
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);
    
    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }/*chiudo for */

  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));

  
  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");


  TH2D* h2_l= new TH2D("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2D* h2_r= new TH2D("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2D* h2_m= new TH2D("h2_m", "histo h2_m",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	h2_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]);
	h2_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]);
	h2_m->Fill((amp_max[3]+amp_max[4])/(2*max),(time[1+LEDi]+time[2+LEDi])/2-time[0]);
      }

  }//chiudo for k



  for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_t;

    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
    histotemp_t=h2_m->ProjectionY("h2_mprojY",k,k);
   
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
  
  
  TCanvas* wf_c =new TCanvas("wf","Plot wf",1800,1000);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);


  TF1* hyp_r = new TF1("hyp_r","[0]+[2]*x+[3]*x**2+[4]*x**3+[5]*x**4",0.8*fit_r->GetParameter(1),1.5*fit_r->GetParameter(1));

  TF1* hyp_l = new TF1("hyp_l","[0]+[2]*x+[3]*x**2+[4]*x**3+[5]*x**4",0.5*fit_l->GetParameter(1),1.5*fit_l->GetParameter(1));
  TF1* hyp_t = new TF1("hyp_t","[1]*x**2+[2]*x+[0]",-0.1,0.65);

  gStyle->SetOptStat("");
   
  wf_c->Divide(3,1);

  wf_c->cd(1);

  h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");

  h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_l","0");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  //graph_l->Draw("SAMEP");
  //hyp_l->DrawF1(0,0.5,"same");


  wf_c->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","0");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  //graph_r->Draw("SAMEP");
  //hyp_r->DrawF1(0,0.5,"same");

 cout <<"________________________________________________"  << endl;
  
  wf_c->cd(3);
  h2_m->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
  h2_m->GetXaxis()->SetTitle("t_left-t_right [ns]");
  h2_m->Draw("COLZ");
  graph_t->Fit("hyp_t","0");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  //graph_t->Draw("SAMEP");
  //hyp_t->DrawF1(txmin,txmax,"same");  


  TH1D* histotemp_lup;
  TH1D* histotemp_ldown;
  TH1D* histotemp_rup;
  TH1D* histotemp_rdown;
  
  Double_t downr=0;
  Double_t upr=10.0;
  
  Double_t downl=0;
  Double_t upl=10.0;
   
  histotemp_rup=h2_r->ProjectionX("h2_rprojYup",h2_r->GetYaxis()->FindBin(downr),h2_r->GetYaxis()->FindBin(upr));
  histotemp_rdown=h2_r->ProjectionX("h2_rprojYdown",h2_r->GetYaxis()->FindBin(rymin_r),h2_r->GetYaxis()->FindBin(downr));
  histotemp_lup=h2_l->ProjectionX("h2_lprojYup",h2_l->GetYaxis()->FindBin(downl),h2_l->GetYaxis()->FindBin(upl));
  histotemp_ldown=h2_l->ProjectionX("h2_lprojYdown",h2_l->GetYaxis()->FindBin(rymin_l),h2_l->GetYaxis()->FindBin(downr));
   
 

  TCanvas* cnv =new TCanvas("cnv","Plot cnv",600,550);
  cnv->Divide(2,2);
  
  cnv->cd(1);
  histotemp_rup->GetXaxis()->SetTitle("amp_max [mV]");
  histotemp_rup->GetYaxis()->SetTitle("counts");
  histotemp_rup->Draw("HIST");
  
  
  cnv->cd(2);
  histotemp_lup->GetXaxis()->SetTitle("amp_max [mV]");
  histotemp_lup->GetYaxis()->SetTitle("counts");
  histotemp_lup->Draw("HIST");
  
  cnv->cd(3);
  histotemp_rdown->GetXaxis()->SetTitle("amp_max [mV]");
  histotemp_rdown->GetYaxis()->SetTitle("counts");
  histotemp_rdown->Draw("HIST");
  
  cnv->cd(4);
  histotemp_ldown->GetXaxis()->SetTitle("amp_max [mV]");
  histotemp_ldown->GetYaxis()->SetTitle("counts");
  histotemp_ldown->Draw("HIST");
  
  
  TCanvas* normalizz =new TCanvas("normalizz","Plot normalized",600,550);
  normalizz->Divide(2,3);
  
  /* histotemp_rup->Scale(1/(histotemp_rup->Integral()));
  histotemp_lup->Scale(1/(histotemp_lup->Integral()));
  histotemp_rdown->Scale(1/(histotemp_rdown->Integral()));
  histotemp_ldown->Scale(1/(histotemp_ldown->Integral()));
  */
  normalizz->cd(1);
  histotemp_rup->Draw("HISTO");
  normalizz->cd(2);
  histotemp_lup->Draw("HISTO");
  normalizz->cd(3);
  histotemp_rdown->Draw("HISTO");
  normalizz->cd(4);
  histotemp_ldown->Draw("HISTO");
  
  cout<<"________________" << histotemp_rup->GetMaximumBin()<<"______"<< histotemp_rdown->GetMaximumBin()<<endl;
  
  TH1D* HistoDiffr = new TH1D("histo_diffr","differenzeR",nbinx,rxmin,rxmax);
  TH1D* HistoDiffl = new TH1D("histo_diffl","differenzeL",nbinx,rxmin,rxmax);
  
  for(i=0;i<nbinx-1;i++){
    HistoDiffr->SetBinContent(i,histotemp_rup->GetBinContent(i)-histotemp_rdown->GetBinContent(i));
    HistoDiffl->SetBinContent(i,histotemp_lup->GetBinContent(i)-histotemp_ldown->GetBinContent(i));
  }
  
  normalizz->cd(5);
  HistoDiffr->Draw("HISTO");
  normalizz->cd(6);
  HistoDiffl->Draw("HISTO");
 
}

