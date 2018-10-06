#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
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

void ProjAWC(TH2D* histo,string name,Int_t nbinx,Float_t index,string time){
  bool EditCanvas=false;

  cout<<"STO QUA_________________________________________________________________________________________________________________________"<<endl;

  TH1D* projection;
  projection=histo->ProjectionY("Proiezione",0,nbinx);
  if(EditCanvas) gROOT->SetBatch(kFALSE);
  TCanvas* projAW = new TCanvas("proiezioni AW","",1200,800);
  TF1* gaussFit = new TF1("projectionFit","gaus");
 
  gaussFit->SetRange(projection->GetBinCenter(projection->GetMaximumBin())-1, projection->GetBinCenter(projection->GetMaximumBin())+1);
  
  gStyle->SetOptFit(1111);
  gStyle->SetStatX(0.5);

  projection->Fit("projectionFit","R");
  gStyle->SetOptFit(1111);
  projection->Draw("");
  
  if(EditCanvas) gPad->WaitPrimitive();
  
  projAW->SaveAs(("HDCRPlot/"+time+"/ProiezioniR"+to_string((Int_t)index)+name+".pdf").c_str());
  projAW->Close();
  
  if(EditCanvas) gROOT->SetBatch(kTRUE);
  
  delete projAW;
  delete gaussFit;
}


//####################################################################################################################################################################################################################### 
void GaussianAmpWalk(TFile* file,Float_t index, Double_t rymin_l,Double_t rymax_l,Double_t rymin_r,Double_t rymax_r,Double_t tymin,Double_t tymax,Float_t* NAWRes,Float_t* errNAWRes){

  TTree* newtree2 = (TTree*)file->Get("digi");

  Float_t amp_max[54], time[54];
  Int_t k,j,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  Double_t max;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi;
  rxmin=0;
  rxmax=0.5;



  const Int_t  nbinx=290,nbiny=150;

  Int_t i;
  Double_t sigma[50],erry[50],cut[50],errx[50];
  
  txmin=-1.1;
  txmax=1.4;
  
  
  Double_t x_r[nbinx],y_r[nbinx], x_l[nbinx],y_l[nbinx],rmsy_l[nbinx],rmsy_r[nbinx];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];
  
  
  Int_t nentries=newtree2->GetEntries(), counter1=0,counter2=0, counter3=0;
  Float_t Times1[nentries],Times2[nentries],Times3[nentries];
  
  for(k=0;k<newtree2->GetEntries();k++){
    Times1[k]=0;
    Times2[k]=0;
    Times3[k]=0;
  }  
  
  TH1D *hr_amp =new TH1D("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1D *hl_amp =new TH1D("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1D *mcp_amp =new TH1D("mcp_amp","histomcp_ampl",nbinx,0.0,1);


  TF1 *fit_r = new TF1("f_r","landau",0.04,1);
  TF1 *fit_l = new TF1("f_l","landau",0.04,1);


  newtree2->SetBranchAddress("amp_max",&amp_max);
  newtree2->SetBranchAddress("time",&time);
  newtree2->SetBranchAddress("LED300",&LED300);
  newtree2->SetBranchAddress("LED100",&LED100);
  newtree2->SetBranchAddress("LED50",&LED50);
  newtree2->SetBranchAddress("LED30",&LED30);

     
  max=4096;


  for(k=0;k<newtree2->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    newtree2->GetEntry(k);

    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
    mcp_amp->Fill(amp_max[0]/max);

  }/*chiudo for */
  
  LEDi=LED300;
  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  //cout << tmax <<endl;
  cout<< max << endl;

  hr_amp->Fit("f_r","R0");
  hl_amp->Fit("f_l","R0");

  bool DebugLand=true;
  
  if(DebugLand){
    TCanvas* LandCanv = new TCanvas("mycanvas","",1200,700);
    LandCanv->Divide(2,1);
    
    LandCanv->cd(1)->SetLogy();
    hl_amp->Draw();
    fit_l->Draw("SAME");
    
    LandCanv->cd(2)->SetLogy();
    hr_amp->Draw();
    fit_r->Draw("SAME");
    
    LandCanv->SaveAs(("HDCRPlot/Landaubis"+to_string((Int_t)index)+".pdf").c_str());
    LandCanv->Close();
  }

  TH2D* h2_l= new TH2D("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2D* h2_r= new TH2D("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2D* h2_t= new TH2D("h2_t", "histo h2_t",nbinx,txmin,txmax,nbiny,tymin,tymax);

  for(k=0;k<newtree2->GetEntries();k++){

    newtree2->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	h2_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]);
	if (amp_max[4]/max < 0.35) h2_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]);
	h2_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]);

      }//chiudo if

  }//chiudo for k

  TH1D* histohampl;
  TH1D* histohampr;
  TF1* gaush_r;
  TF1* gaush_l;

  TCanvas* ProiezioniTemp;
  Int_t npoInt_t=50;
  
  gSystem->Exec(("mkdir HDCRPlot/Gauss"+to_string((Int_t)index)).c_str());
  
  for(i=0;i<npoInt_t;i++){
    x_r[i]=0;
    x_l[i]=0;
    y_l[i]=0;
    y_r[i]=0;
    rmsy_r[i]=0;
    rmsy_l[i]=0;
  }
    
  for(i=0;i<npoInt_t;i++){
    
    histohampl=h2_l->ProjectionY("FettaAmpL",h2_l->FindBin(rxmin+(rxmax-rxmin)/npoInt_t*i),h2_l->FindBin(rxmin+(rxmax-rxmin)/npoInt_t*(i+1)));
    histohampr=h2_r->ProjectionY("FettaAmpR",h2_r->FindBin(rxmin+(rxmax-rxmin)/npoInt_t*i),h2_r->FindBin(rxmin+(rxmax-rxmin)/npoInt_t*(i+1)));
  
    x_r[i]=rxmin+(rxmax-rxmin)/npoInt_t*i+((rxmax-rxmin)/npoInt_t)/2;
    x_l[i]=rxmin+(rxmax-rxmin)/npoInt_t*i+((rxmax-rxmin)/npoInt_t)/2;
    
    gaush_r = new TF1("gaushr","gaus",7.7,9);
    gaush_l = new TF1("gaushl","gaus",7.7,9);
    ProiezioniTemp = new TCanvas("mycanv","",1200,800);
    ProiezioniTemp->Divide(2,1);
    
    ProiezioniTemp->cd(1);
    
    histohampl->Draw();
    histohampl->Fit("gaushl","R");
    y_l[i]=gaush_l->GetParameter(1);
    rmsy_l[i]=gaush_l->GetParError(1);

    ProiezioniTemp->cd(2);
    
    histohampr->Draw();
    histohampr->Fit("gaushr","R");
    y_r[i]=gaush_r->GetParameter(1);
    rmsy_r[i]=gaush_r->GetParError(1);
   
    
    ProiezioniTemp->SaveAs(("HDCRPlot/Gauss"+to_string((Int_t)index)+"/GaussianeAmp"+to_string(i)+".pdf").c_str());
    ProiezioniTemp->Close();
    delete ProiezioniTemp;
    delete gaush_l;
    delete gaush_r;
    //delete mcp_amp;
  
  }//chiudo for
  
  TGraphErrors* graphGAmp_r=new TGraphErrors(npoInt_t,x_r,y_r,0,rmsy_r);
  TGraphErrors* graphGAmp_l=new TGraphErrors(npoInt_t,x_l,y_l,0,rmsy_l);

  TF1* FitLog_r = new TF1("FitLog_r","[0]+[2]*log(x+[1])",0.8*fit_r->GetParameter(1),1.5*fit_r->GetParameter(1));
  TF1* FitLog_l = new TF1("FitLog_l","[0]+[2]*log(x+[1])",0.5*fit_l->GetParameter(1),1.2*fit_l->GetParameter(1));
  
  gStyle->SetOptFit();
  
     /* SetParameters*/
  FitLog_l->SetParameter(0, 8.51);
  FitLog_l->SetParameter(1, 5);
  FitLog_l->SetParameter(2, 1.2);
  
  FitLog_r->SetParameter(0, 8.51);
  FitLog_r->SetParameter(1, 5);
  FitLog_r->SetParameter(2, 1.2);
  

  bool controlFit=false;
  
  if(controlFit)gROOT->SetBatch(kFALSE);
 
  TCanvas* CanvAmpGauss = new TCanvas("canv amplitude gauss","",1200,800);
  CanvAmpGauss->Divide(3,2);
  
  CanvAmpGauss->cd(1);
  h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
  h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_l->Draw("COLZ");
  graphGAmp_l->SetMarkerStyle(8);
  graphGAmp_l->SetMarkerSize(.8);
  graphGAmp_l->Draw("SAMEP");
  graphGAmp_l->Fit("FitLog_l");

  CanvAmpGauss->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graphGAmp_r->SetMarkerStyle(8);
  graphGAmp_r->SetMarkerSize(.8);
  graphGAmp_r->Draw("SAMEP");
  graphGAmp_r->Fit("FitLog_r");
  
  CanvAmpGauss->cd(3);
  h2_t->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_t->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_t->Draw("COLZ");
  


  rymin_lc=rymin_l-FitLog_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))+FitLog_l->GetParameter(0);
  rymax_lc=rymax_l-FitLog_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))+FitLog_l->GetParameter(0);
  rymin_rc=rymin_r-FitLog_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))+FitLog_r->GetParameter(0);
  rymax_rc=rymax_r-FitLog_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))+FitLog_r->GetParameter(0);
  tymin_c=tymin-(FitLog_l->Eval(0.25)-FitLog_l->GetParameter(0)+FitLog_r->Eval(0.25)-FitLog_r->GetParameter(0))/2;
  tymax_c=tymax-(FitLog_l->Eval(0.25)-FitLog_l->GetParameter(0)+FitLog_r->Eval(0.25)-FitLog_r->GetParameter(0))/2;
  


  TH2D* hc_l= new TH2D("hc_l", "histo hc_l",nbinx,rxmin,rxmax,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_r= new TH2D("hc_r", "histo hc_r",nbinx,rxmin,rxmax,nbiny,rymin_rc,rymax_rc);
  TH2D* hc_t= new TH2D("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);

  for(k=0;k<newtree2->GetEntries();k++){
    
    newtree2->GetEntry(k);
    
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	hc_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]-FitLog_l->Eval(amp_max[3]/max)+FitLog_l->GetParameter(0));
        hc_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]-FitLog_r->Eval(amp_max[4]/max)+FitLog_r->GetParameter(0));
	hc_t->Fill(time[1+LEDi]-time[2+LEDi],(time[1+LEDi]+time[2+LEDi])/2-time[0]-(FitLog_r->Eval(amp_max[4]/max)-FitLog_r->GetParameter(0)+FitLog_l->Eval(amp_max[3]/max)-FitLog_l->GetParameter(0))/2);
      }
    
  }//chiudo for k
  gSystem->Exec("mkdir HDCRPlot/NAW");
  
  
  
  CanvAmpGauss->cd(4);
  hc_l->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  hc_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  hc_l->Draw("COLZ");

  CanvAmpGauss->cd(5);
  hc_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  hc_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  hc_r->Draw("COLZ");

  CanvAmpGauss->cd(6);
  hc_t->GetYaxis()->SetTitle("t_{ave}-t_MCP [ns]");
  hc_t->GetXaxis()->SetTitle("t_{left}-t_{right} [ns]");
  hc_t->Draw("COLZ");
  
  if(controlFit) gPad->WaitPrimitive();
  CanvAmpGauss->SaveAs(("HDCRPlot/GaussianeAmp"+to_string((Int_t)index)+".pdf").c_str());
  CanvAmpGauss->Close();
  if(controlFit) gROOT->SetBatch(kTRUE);

  
  TH1D* tProjection = hc_t->ProjectionY("ProjectionTave",0,nbinx);
  TCanvas* RisDef = new TCanvas("RisDef","",1200,700);

  TF1* fitgaussiano = new TF1("fitgaussiano","gaus");
  gStyle->SetOptFit(0111);
  fitgaussiano->SetRange(tProjection->GetBinCenter(tProjection->GetMaximumBin())-1, tProjection->GetBinCenter(tProjection->GetMaximumBin())+1);
  
  tProjection->Fit("fitgaussiano","R");
  tProjection->Draw("SAME");
  
  RisDef->SaveAs(("HDCRPlot/NewAWTimeRes"+to_string((Int_t)index)+".pdf").c_str());
  RisDef->Close();
  
  *NAWRes=fitgaussiano->GetParameter(2);
  *errNAWRes=fitgaussiano->GetParError(2);

  ProjAWC(h2_r,"uncorrNAW",nbinx,index,"NAW");
  ProjAWC(hc_r,"corrNAW",nbinx,index,"NAW");

  delete hr_amp;
  delete hl_amp;
  delete h2_l;
  delete h2_r;
  delete h2_t;
  delete hc_r;
  delete hc_l;
  delete hc_t;
  delete RisDef;
  delete fitgaussiano;

}//CHIUDO FUNZIONE

//#######################################################################################################################################################################################################################

void plotWF_lsig(TFile* file,Float_t* Resol,Float_t* errResol,Float_t index,Float_t* Resolc,Float_t* errResolc){

  TTree* newtree = (TTree*)file->Get("digi");
  

  Float_t amp_max[54], time[54];
  Int_t k,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi;
  bool debug=false;
  Double_t max=0,tmax=0;

  rxmin=-0.3;
  rxmax=0.6;


  
  const Int_t  nbinx=100,nbiny=50;

  rymin_l=6;
  rymax_l=10;
  rymin_r=6;
  rymax_r=10;
  tymin=7;
  tymax=10;
  txmin =-0.5;
  txmax=0.8;
 
  Double_t RMS[3][nbinx];
  Double_t x_r[nbinx],y_r[nbiny], x_l[nbinx],y_l[nbiny],rmsy_l[nbiny],rmsy_r[nbiny];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];


  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);

  TF1 *fit_r = new TF1("f_r","landau",0.10,1);
  TF1 *fit_l = new TF1("f_l","landau",0.10,1);


  newtree->SetBranchAddress("amp_max",&amp_max);
  newtree->SetBranchAddress("time",&time);

  newtree->SetBranchAddress("LED300",&LED300);
  newtree->SetBranchAddress("LED100",&LED100);
  newtree->SetBranchAddress("LED50",&LED50);
  newtree->SetBranchAddress("LED30",&LED30);
  
  newtree->GetEntry(1);
  
  LEDi=LED300;

  max=4096;

  for(k=0;k<newtree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    newtree->GetEntry(k);

    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }/*chiudo for */

  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));

  
  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");


  TH2D* h2_l= new TH2D("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2D* h2_r= new TH2D("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2D* histotempi= new TH2D("histotempi", "histo histotempi",nbinx,txmin,txmax,nbiny,tymin,tymax);

  histotempi->Reset();
  
  newtree->GetEntry(9);
  cout<<time[1+LEDi]-time[2+LEDi]<<endl;
  
  for(k=0;k<newtree->GetEntries();k++){
    
    newtree->GetEntry(k);
   
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75){
     
      histotempi->Fill(time[1+LEDi]-time[2+LEDi],(time[1+LEDi]+time[2+LEDi])/2-time[0]);
    }
    
  }//chiudo for k
  
  cout <<    histotempi->GetEntries() << endl;
  
  // for(k=0;k<nbinx;k++){
  TH1D* histotemp_l;
  TH1D* histotemp_r;
  TH1D* histotemp_m;
  
  histotemp_l=h2_l->ProjectionY(((string)"h2_lprojY"+to_string((Int_t)index)).c_str(),0,nbinx);
  histotemp_r=h2_r->ProjectionY(((string)"h2_rprojY"+to_string((Int_t)index)).c_str(),0,nbinx);
  histotemp_m=histotempi->ProjectionY(((string)"h2_tprojY"+to_string((Int_t)index)).c_str(),0,nbinx);
  

  TCanvas* wf_c =new TCanvas("wfbis","Plot wfbis",600,550);
  
  TF1* g_r = new TF1(((string)"g_r"+to_string((Int_t)index)).c_str(),"gaus",0,21);
  TF1* g_l = new TF1(((string)"g_l"+to_string((Int_t)index)).c_str(),"gaus",0,21);

  TF1* g_m = new TF1(((string)"g_m"+to_string((Int_t)index)).c_str(),"gaus",7,10);
  //  hyp_r->SetParameter(0,10);
  // hyp_l->SetParameter(0,10);
  histotemp_m->Fit(((string)"g_m"+to_string((Int_t)index)).c_str(),"0");
 
  
  gStyle->SetOptStat("");
  gStyle->SetOptFit();
  histotemp_r->SetLineColor(kRed);
  histotemp_l->SetLineColor(kBlue);
  

  *Resol = g_m->GetParameter(2);
  *errResol = g_m->GetParError(2);

  histotemp_m->GetXaxis()->SetTitle("t_ave-t_MCP");
  histotemp_m->GetYaxis()->SetTitle("counts");

  
  
  histotemp_m->Draw();

  g_m->Draw("same");


  wf_c->SaveAs(("HDCRPlot/LargeGaus"+to_string((Int_t)index)+".pdf").c_str());
  wf_c->Close();

  
 for(k=0;k<nbinx;k++){
   
    TH1D* histotemp_t;
   

    histotemp_t=histotempi->ProjectionY("hc_tprojY",k,k);
   
    xt[k]=rxmin+(Float_t)(rxmax-(rxmin))/nbinx*k;
    yt[k]=histotemp_t->GetMean();
    RMS[2][k]= histotemp_t->GetMeanError();

            
    delete histotemp_t;
  }//chiudo for k
  
  TH2D* hc_tdiff= new TH2D("hc_tdiff", "histo hc_tdiff",nbinx,txmin,txmax,nbiny,tymin,tymax);
  TF1* fit_tdiff = new TF1("fit_tdiff","[0]+[1]*x+[2]*x**2",-0.3,0.6);
  TGraphErrors* graph_tc = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);

  graph_tc->Fit("fit_tdiff","R0");
  TCanvas * prova = new TCanvas("","",600,500);
  histotempi->Draw("COLZ");
  graph_tc->Draw("SAME");
  fit_tdiff->Draw("same");
  //prova->SaveAs(((string)"HDCRPlot/prova"+ to_string(index)+(string)".png").c_str());
  prova->Close();
       

  

  for(k=0;k<newtree->GetEntries();k++){
    
    
    newtree->GetEntry(k);
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	hc_tdiff->Fill(time[1+LEDi]-time[2+LEDi],(time[1+LEDi]+time[2+LEDi])/2-time[0]-fit_tdiff->Eval(time[1+LEDi]-time[2+LEDi])+fit_tdiff->GetParameter(0));
      }//chiudo for k
    
    
  }
  
  
  TF1* retta = new TF1("retta","[0]+[1]*x",txmin+0.3,txmax-0.3);
  for(k=0;k<nbinx;k++){
    TH1D* histotemp_t;
    histotemp_t=hc_tdiff->ProjectionY("hc_tprojY",k,k);
    
    yt[k]=histotemp_t->GetMean();
    RMS[2][k]= histotemp_t->GetMeanError();
    
    delete histotemp_t;
  }//chiudo for k
  
  TH1D* histo_ctdiff;
  TGraphErrors* graph_tcdiff = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);
  graph_tcdiff->Fit("retta","0R");
  histo_ctdiff = hc_tdiff->ProjectionY("histo_ctdiff",0,nbinx);
  
  TF1* gaus_ctdiff = new TF1(((string)"gaus_ctdiff"+to_string((Int_t)index)).c_str(),"gaus", 7, 10);
  TCanvas* tdiff = new TCanvas("tdiff","plot_tdiff",600,550);
  TLegend* l2=new TLegend(0.1,0.7,0.48,0.9);
  g_m->SetLineColor(kBlack);
  
  // g_m->Draw("same");
  gaus_ctdiff->SetLineColor(kGreen);
  gaus_ctdiff->SetParameter(0,histo_ctdiff->GetBinContent(histo_ctdiff->GetMaximumBin()));
  gaus_ctdiff->SetParLimits(0,histo_ctdiff->GetBinContent(histo_ctdiff->GetMaximumBin())-5,TMath::Infinity());
  
  histo_ctdiff->Fit(((string)"gaus_ctdiff"+to_string((Int_t)index)).c_str(),"0MR");
  histo_ctdiff->SetLineColor(kGreen);
  histo_ctdiff->Draw();
  gaus_ctdiff->Draw("same");
  histotemp_m->Draw("same");
  
  *Resolc=gaus_ctdiff->GetParameter(2);
  *errResolc=gaus_ctdiff->GetParError(2);
  
  l2->SetHeader("t_{ave}-t_{MCP} distrib");
  l2->AddEntry(histotemp_m,"t_ave-t_MCP");
  l2->AddEntry(g_m,("#sigma="+to_string(g_m->GetParameter(2))).c_str());
  l2->AddEntry(histo_ctdiff,"t_ave-t_MCP(tdiff corr)");
  l2->AddEntry(gaus_ctdiff,("#sigma="+to_string(gaus_ctdiff->GetParameter(2))).c_str());
  //  l2->Draw();
  tdiff->SaveAs(("HDCRPlot/Gaussianenowalk"+to_string((Int_t)index)+".pdf").c_str());
  tdiff->Close();
  
  delete histo_ctdiff;
  delete gaus_ctdiff;
  delete histotempi;
  delete h2_r;
  delete h2_l;
  delete fit_tdiff;
  delete hr_amp;
  delete hl_amp;
  delete hc_tdiff;
  
  
}
//####################################################################################################################################################################################################################### 
void Proiezione(TH2D* hl,TH2D* hr,TH2D* ht,Int_t nxbin,Float_t DCR){
    TH1D* proj[3];
    proj[0]=hl->ProjectionY("uncorrR",0,nxbin);
    proj[1]=hr->ProjectionY("uncorrL",0,nxbin);
    proj[2]=ht->ProjectionY("uncorrT",0,nxbin);
    
    TCanvas* pro = new TCanvas("leproj","",500,400);
    pro->Divide(3,1);
    pro->cd(1);
    proj[0]->GetXaxis()->SetTitle("t_{left}-t_{MCP}");
    proj[0]->GetYaxis()->SetTitle("counts");
    proj[0]->SetTitle("SiPM_left");
    proj[0]->Draw();
    
    pro->cd(2);
    proj[1]->GetXaxis()->SetTitle("t_{right}-t_{MCP}");
    proj[1]->GetYaxis()->SetTitle("counts");
    proj[1]->SetTitle("SiPM_right");
    proj[1]->Draw("SAME");
    
    pro->cd(3);
    proj[2]->GetXaxis()->SetTitle("t_{ave}-t_{MCP}");
    proj[2]->GetYaxis()->SetTitle("counts");
    proj[2]->SetTitle("t_{ave}");
    proj[2]->Draw("SAME");
    
    pro->SaveAs(("HDCRPlot/ProiezioniLR"+to_string((Int_t)DCR)+".pdf").c_str());
    
    
    pro->Close();
    delete pro;
    
    
  }
//####################################################################################################################################################################################################################### 


void Ris(TFile* file,Float_t* Resol,Float_t* errResol,Float_t* ResolHisto,Float_t index,Float_t* NAWRes,Float_t* errNAWRes){

  TTree* newtree3 = (TTree*)file->Get("digi");

  Float_t amp_max[54], time[54];
  Int_t k,j,maxbin_l,maxbin_r,maxbin_t;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  bool debug=false;
  bool blind=true;
  Double_t max=0;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi;
  rxmin=0;
  rxmax=0.5;



  const Int_t  nbinx=290,nbiny=200;

  Int_t i;
  Double_t sigma[50],erry[50],cut[50],errx[50];
  
  txmin=-1.1;
  txmax=1.4;
  
  
  Double_t x_r[nbinx],y_r[nbinx], x_l[nbinx],y_l[nbinx],rmsy_l[nbinx],rmsy_r[nbinx];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];
  
  
  Int_t nentries=newtree3->GetEntries(), counter1=0,counter2=0, counter3=0;
  Float_t Times1[nentries],Times2[nentries],Times3[nentries];
  
  for(k=0;k<newtree3->GetEntries();k++){
    Times1[k]=0;
    Times2[k]=0;
    Times3[k]=0;
  }  
  
  TH1D *hr_amp =new TH1D("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1D *hl_amp =new TH1D("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1D *mcp_amp =new TH1D("mcp_amp","histomcp_ampl",nbinx,0.0,1);


  TF1 *fit_r = new TF1("f_r","landau",0.04,1);
  TF1 *fit_l = new TF1("f_l","landau",0.04,1);


  newtree3->SetBranchAddress("amp_max",&amp_max);
  newtree3->SetBranchAddress("time",&time);
  newtree3->SetBranchAddress("LED300",&LED300);
  newtree3->SetBranchAddress("LED100",&LED100);
  newtree3->SetBranchAddress("LED50",&LED50);
  newtree3->SetBranchAddress("LED30",&LED30);

  newtree3->GetEntry(3);
  LEDi=LED300;
  for(k=0;k<newtree3->GetEntries();k++){
    newtree3->GetEntry(k);
    
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
    

  tymin=mean3-1.3*rms3;
  tymax=mean3+1.3*rms3;
 
  max=4096;

  
  cout<<"Exit"<<endl;
  for(k=0;k<newtree3->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    newtree3->GetEntry(k);

    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
    mcp_amp->Fill(amp_max[0]/max);

  }/*chiudo for */

  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));


  //cout << tmax <<endl;
  cout<< max << endl;

  hr_amp->Fit("f_r","R0");
  hl_amp->Fit("f_l","R0");

  bool DebugLand=true;
  
  if(DebugLand){
    TCanvas* LandCanv = new TCanvas("mycanvas","",1200,700);
    LandCanv->Divide(2,1);
    
    LandCanv->cd(1)->SetLogy();
    hl_amp->Draw();
    fit_l->Draw("SAME");
    
    LandCanv->cd(2)->SetLogy();
    hr_amp->Draw();
    fit_r->Draw("SAME");
    
    LandCanv->SaveAs(("HDCRPlot/Landau"+to_string((Int_t)index)+".pdf").c_str());
    LandCanv->Close();
  }

  TH2D* h2_l= new TH2D("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2D* h2_r= new TH2D("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2D* h2_t= new TH2D("h2_t", "histo h2_t",nbinx,txmin,txmax,nbiny,tymin,tymax);

  for(k=0;k<newtree3->GetEntries();k++){

    newtree3->GetEntry(k);

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


  TF1* hyp_r = new TF1("hyp_r","[0]+[2]*x + [3]*x**2+[4]*x**3+[5]*x**4+[6]*x**5+ [7]/x",0.8*fit_r->GetParameter(1),1.5*fit_r->GetParameter(1));
  TF1* hyp_l = new TF1("hyp_l","[0]+[2]*x + [3]*x**2+[4]*x**3+[5]*x**4+[6]*x**5+ [7]/x",0.5*fit_l->GetParameter(1),1.2*fit_l->GetParameter(1));


  TF1* hyp_t = new TF1("hyp_t","[1]*x**2+[2]*x+[0]",-0.1,0.65);

  Proiezione(h2_l,h2_r,h2_t,nbinx,index);
  
  gStyle->SetOptStat("");
   
  wf_c->Divide(3,2);

  wf_c->cd(1);

  h2_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");

  h2_l->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_l->Draw("COLZ");
  graph_l->Fit("hyp_l","0");
  graph_l->SetMarkerStyle(8);
  graph_l->SetMarkerSize(.5);
  graph_l->Draw("SAMEP");
  hyp_l->DrawF1(0,0.5,"same");


  wf_c->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","0");
  graph_r->SetMarkerStyle(8);
  graph_r->SetMarkerSize(.5);
  graph_r->Draw("SAMEP");
  hyp_r->DrawF1(0,0.5,"same");
  
  wf_c->cd(3);
  h2_t->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
  h2_t->GetXaxis()->SetTitle("t_left-t_right [ns]");
  h2_t->Draw("COLZ");
  graph_t->Fit("hyp_t","0");
  graph_t->SetMarkerStyle(8);
  graph_t->SetMarkerSize(.5);
  graph_t->Draw("SAMEP");
  hyp_t->DrawF1(txmin,txmax,"same");


  rymin_lc=rymin_l-hyp_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))/*+hyp_l->GetParameter(0)*/;
  rymax_lc=rymax_l-hyp_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))/*+hyp_l->GetParameter(0)*/;
  rymin_rc=rymin_r-hyp_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))/*+hyp_r->GetParameter(0)*/;
  rymax_rc=rymax_r-hyp_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))/*+hyp_r->GetParameter(0)*/;
  tymin_c=tymin-(hyp_l->Eval((fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2)+fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2)/2))/*-hyp_l->GetParameter(0)*/+hyp_r->Eval((fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2)+fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))/2)/*-hyp_r->GetParameter(0)*/)/2;
  tymax_c=tymax-(hyp_l->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))/*-hyp_l->GetParameter(0)*/+hyp_r->Eval((fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2)+fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))/2)/*-hyp_r->GetParameter(0)*/)/2;


  
  TH2D* hc_l= new TH2D("hc_l", "histo hc_l",nbinx,rxmin,rxmax,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_r= new TH2D("hc_r", "histo hc_r",nbinx,rxmin,rxmax,nbiny,rymin_rc,rymax_rc);
  TH2D* hc_t= new TH2D("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  TH2D* hc_tdiff= new TH2D("hc_tdiff", "histo hc_tdiff",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  

  
   for(k=0;k<newtree3->GetEntries();k++){

    newtree3->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {

	hc_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/max));
        hc_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/max));
	hc_t->Fill(time[1+LEDi]-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_l->Eval(amp_max[3]/max),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)+hyp_l->Eval(amp_max[3]/max))/2);

      }
    
  }//chiudo for k

   gSystem->Exec("mkdir HDCRPlot/OAW");
   ProjAWC(h2_r,"uncorrOAW",nbinx,index,"OAW");
   ProjAWC(hc_r,"corrOAW",nbinx,index,"OAW");
   
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
   TF1* fit_tdiff = new TF1("fit_tdiff","[0]+[1]*x",-0.2,0.4);
  
   graph_tc->Fit("fit_tdiff","R0L");
   
   for(k=0;k<newtree3->GetEntries();k++){
     newtree3->GetEntry(k);
     
     if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
       {

	 hc_tdiff->Fill(time[1+LEDi]-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_l->Eval(amp_max[3]/max),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)+hyp_l->Eval(amp_max[3]/max))/2-fit_tdiff->Eval(time[1+LEDi]-time[2+LEDi]+hyp_r->Eval(amp_max[4]/max)-hyp_l->Eval(amp_max[3]/max)));

       }
     
   }//chiudo for k
   
  

 
   
   TF1* retta = new TF1("retta","[0]+[1]*x",-0.2,0.4);
   for(k=0;k<nbinx;k++){
     TH1D* histotemp_t;
     histotemp_t=hc_tdiff->ProjectionY("hc_tprojY",k,k);
     
   
     yt[k]=histotemp_t->GetMean();
     RMS[2][k]= histotemp_t->GetMeanError();
     
     delete histotemp_t;
     
     
   }//chiudo for k
    
   TGraphErrors* graph_tcdiff = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);
   graph_tcdiff->Fit("retta","0");
   
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

   hc_tdiff->GetYaxis()->SetTitle("t_ave-t_MCP [ns]");
   hc_tdiff->GetXaxis()->SetTitle("t_left-t_right [ns]");
   hc_tdiff->Draw("COLZ");
   graph_tcdiff->SetMarkerStyle(8);
   graph_tcdiff->SetMarkerSize(.5);
   graph_tcdiff->Draw("P");
   
   wf_c->SaveAs(("HDCRPlot/Controllo"+to_string((Int_t)index)+".pdf").c_str());
   wf_c->Close();

   
   TH1D* histo_cl;
   TH1D* histo_cr;
   TH1D* histo_ct;
   TH1D* histo_ctdiff;
   
   histo_cl = hc_l->ProjectionY("histo_cl",0,nbinx);
   histo_cr = hc_r->ProjectionY("histo_cr",0,nbinx);
   
   histo_ct = hc_t->ProjectionY("histo_ct",0,nbinx);
   histo_ctdiff = hc_tdiff->ProjectionY("histo_ctdiff",0,nbinx);

   TF1* gaus_cl = new TF1("gaus_cl","gaus");
   TF1* gaus_cr = new TF1("gaus_cr","gaus");
   TF1* gaus_ct = new TF1("gaus_ct","gaus", -0.5, 2);
   TF1* gaus_ctdiff = new TF1("gaus_ctdiff","gaus", -0.5, 2);

   histo_ct->SetLineColor(kBlack);
   histo_cl->SetLineColor(kBlue);
   histo_cr->SetLineColor(kRed);
   gaus_ct->SetLineColor(kBlack);
      
   TCanvas* tdiff = new TCanvas("tdiff","plot_tdiff",600,550);
   TLegend* l2=new TLegend(0.1,0.7,0.48,0.9);
  
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   

   gaus_ct->SetParameter(0,80);
   gaus_ct->SetParameter(1, 0.0);
   gaus_ct->SetParameter(2,0.03);
   histo_ct->Fit("gaus_ct","R");
   
   gaus_ctdiff->SetLineColor(kGreen);
   gaus_ctdiff->SetParameter(0,80);
   gaus_ctdiff->SetParameter(1,0.0);
   gaus_ctdiff->SetParameter(2,0.04);
   
   histo_ctdiff->Fit("gaus_ctdiff","R");
   histo_ctdiff->SetLineColor(kGreen);
   histo_ctdiff->Draw();   
   histo_ct->Draw("SAME");
   
   *Resol=gaus_ctdiff->GetParameter(2);
   *errResol=gaus_ctdiff->GetParError(2);
   *ResolHisto=histo_ctdiff->GetRMS();
   
   l2->SetHeader("t_{ave}-t_{MCP} distrib");
   l2->AddEntry(histo_ct,"t_ave-t_MCP");
   l2->AddEntry(gaus_ct,("#sigma="+to_string(gaus_ct->GetParameter(2))).c_str());
   l2->AddEntry(histo_ctdiff,"t_ave-t_MCP(tdiff corr)");
   l2->AddEntry(gaus_ctdiff,("#sigma="+to_string(gaus_ctdiff->GetParameter(2))).c_str());
   l2->AddEntry(histo_ctdiff,("#sigma_{histo}="+to_string(histo_ctdiff->GetRMS())).c_str());
   l2->Draw();
   tdiff->SaveAs(("HDCRPlot/Gaussiane"+to_string((Int_t)index)+".pdf").c_str());
   tdiff->Close();

   

   bool ConfrontoTdiff=true;
   
   if(ConfrontoTdiff){
     TCanvas* ConfCanv = new TCanvas("confronto","",1200,600);
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
     graph_tcdiff->Fit("retta","0");
     retta->DrawF1(-1,1,"SAME");
     
     ConfCanv->SaveAs(("HDCRPlot/Confronto"+to_string((Int_t)index)+".pdf").c_str());
     ConfCanv->Close();
   }

   delete h2_r;
   delete h2_l;
   delete h2_t;
   delete hc_r;
   delete hc_l;
   delete hc_t;
   delete hc_tdiff;
   delete hr_amp;
   delete hl_amp;
   delete mcp_amp;
  
   GaussianAmpWalk(file,index, rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,NAWRes,errNAWRes);
}

//####################################################################################################################################################################################################################### 
void RisVsDcr(string version){
  //Int_t nfiles=19;
  Int_t nfiles=5;
  
  TFile* myfile[nfiles];

  Float_t sigma[nfiles],errsigma[nfiles],sigmaHisto[nfiles], lsigma[nfiles],errlsigma[nfiles],lcsigma[nfiles],errlcsigma[nfiles];
  Float_t NAWRes[nfiles],errNAWRes[nfiles];


  Int_t i;
  Float_t bias[nfiles],NINOthr[nfiles],DCR[nfiles];
  Float_t biasPl[nfiles],NINOthrPl[nfiles],DCRPl[nfiles];
  
  gSystem->Exec("rm -r -f HDCRPlot");
  gSystem->Exec("mkdir HDCRPlot");

  for(i=0;i<nfiles;i++){
    if(version=="old")myfile[i]=TFile::Open(("Pd/DCR10072/"+to_string(i)+".root").c_str());    
    if(version=="new")myfile[i]=TFile::Open(("Pd/NewRecoDcr/"+to_string(i)+".root").c_str());    
    
  }
  
  TTree* info[nfiles];
  
  for(i=0;i<nfiles;i++){
    info[i]=(TTree*)myfile[i]->Get("info");
    info[i]->SetBranchAddress("SiPMCurrent_bar",&DCR[i]);
    info[i]->SetBranchAddress("NINOthr_bar",&NINOthr[i]);
    info[i]->SetBranchAddress("Vbias_bar",&bias[i]);
  }
  
  gROOT->SetBatch(kTRUE);
  
  for(i=0;i<nfiles;i++){
  
    info[i]->GetEntry(1);
    
    biasPl[i]=bias[i];
    NINOthrPl[i]= NINOthr[i];
    DCRPl[i]=DCR[i];

    Ris(myfile[i],&sigma[i],&errsigma[i],&sigmaHisto[i],DCRPl[i],&NAWRes[i],&errNAWRes[i]);
    
    plotWF_lsig(myfile[i],&lsigma[i],&errlsigma[i],DCRPl[i],&lcsigma[i],&errlcsigma[i]);
  }
  
  gROOT->SetBatch(kFALSE);
  

  TCanvas* defcanv = new TCanvas("RisvsDCR","",600,400);
  //defcanv->SetLogy();
  defcanv->SetLogy();
  TGraphErrors* graph = new TGraphErrors(nfiles,DCRPl,sigma,0,errsigma);
  graph->SetMarkerStyle(8);
  graph->SetMarkerSize(.8);
  graph->SetTitle("NINOthr=100 d.u. bias SiPM=72 ");
  
  
  TGraphErrors* graphHisto = new TGraphErrors(nfiles,DCR,sigmaHisto,0,0);

  graphHisto->GetXaxis()->SetTitle("DCR [#muA]");
  graphHisto->GetYaxis()->SetTitle("#sigma_{t_{ave}} [ns]");
  graphHisto->SetMarkerStyle(8);
  graphHisto->SetMarkerSize(.8);
  graphHisto->SetMarkerColor(kBlue);

  graph->SetLineColor(kBlack);
  graphHisto->GetXaxis()->SetLimits(-40.0,2100.0);
  
  TLegend* legenda= new TLegend();
  legenda->SetHeader("#sigma(DCR)");
  legenda->AddEntry(graph,"GaussRMS");
  legenda->AddEntry(graphHisto,"HistoRMS");
  
  graphHisto->Draw("AP");
  graph->Draw("SAMEP");
  legenda->Draw();
  defcanv->SaveAs("HDCRPlot/plot.pdf");

  //defcanv->Close();
  /*  TLatex* tex[nfiles];
  
   for(i=0;i<nfiles;i++){
    
    tex[i]= new TLatex(DCRPl[i],sigma[i],( (to_string((Int_t)NINOthrPl[i]))+"   "+to_string((Int_t)biasPl[i])).c_str());
    tex[i]->SetTextSize(0.035);
    tex[i]->SetLineWidth(2);
    tex[i]->SetLineColor(i);
    tex[i]->Draw();
     
    }*/
  

   TCanvas* nocorr_sigma = new TCanvas("nocorr_sigma","nocorr_#sigma plot",600,550);
   TGraphErrors* ncsigma = new TGraphErrors(nfiles,DCR,lsigma,0,errlsigma);
   TGraphErrors* csigma = new TGraphErrors(nfiles,DCR,lcsigma,0,errlcsigma);
   TGraphErrors* NAWGraph = new TGraphErrors(nfiles,DCR,NAWRes,0,errNAWRes);
   
   TLegend* l3 = new TLegend();
   
   //nocorr_sigma->cd();
    ncsigma->GetXaxis()->SetTitle("DCR [#muA]");

    ncsigma->GetXaxis()->Set(8,-70,2070);
    ncsigma->GetYaxis()->SetRangeUser(0.03,0.25);
    ncsigma->GetYaxis()->SetTitle("#sigma_{t_{ave}}(ns)");
    ncsigma->SetMarkerStyle(8);
    graph->SetMarkerColor(kRed);
    
    ncsigma->SetMarkerSize(.8);
    csigma->SetMarkerStyle(8);
    csigma->SetMarkerSize(.8);
    graph->SetMarkerSize(.8);
    csigma->SetMarkerColor(kBlue);
    NAWGraph->SetMarkerStyle(24);
    NAWGraph->SetMarkerSize(.8);
    


    ncsigma->Draw("AP");
    csigma->Draw("SAMEP");
    graph->Draw("SAMEP");
    NAWGraph->Draw("SAMEP");
    l3->AddEntry(ncsigma,"uncorr. #sigma ","P");
    l3->AddEntry(csigma,"tdiff_corr. #sigma ","P");
    l3->AddEntry(graph,"amp walk+ tdiff corr. #sigma","P");
    l3->AddEntry(NAWGraph,"new amp walk #sigma","P");
    l3->Draw();
 


    for(i=0;i<nfiles;i++){
      cout << DCRPl[i] << "   " <<sigmaHisto[i] << "  " << sigma[i] << "   " << errsigma[i]<< "      ____      "<<NAWRes[i]<<"    " << errNAWRes[i] << endl;
    }
}
//####################################################################################################################################################################################################################### 
