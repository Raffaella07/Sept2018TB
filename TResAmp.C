//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void TResAmp(const char * filename){


  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");



  Float_t amp_max[54], time[54];
  int k,j,maxbin_l,maxbin_r,maxbin_t;
  Double_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax,txmin,txmax,tymin_c,tymax_c,rymin_lc,rymax_lc,rymin_rc,rymax_rc;
  bool debug=false;
  bool blind=true;
  Double_t max=0;
  Int_t LED300,LED100,LED50,LED30;
  Int_t LEDi;
  rxmin=0.8;
  rxmax=1.5;



  const Int_t  nbinx=200,nbiny=250;

 int i;
 Double_t sigma[50],erry[50],cut[50],errx[50];
 Int_t nentries=digiTree->GetEntries();
 Float_t Times1[nentries],Times2[nentries],Times3[nentries];

  
  txmin=-0.3;
  txmax=0.8;


  Double_t x_r[nbinx],y_r[nbinx], x_l[nbinx],y_l[nbinx],rmsy_l[nbinx],rmsy_r[nbinx];
  Double_t xt[nbinx],yt[nbinx],rmsyt[nbinx];
  Double_t RMS[3][nbinx];


  TH1D *hr_amp =new TH1D("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1D *hl_amp =new TH1D("hl_amp","histos_ampl",nbinx,0.0,1);
  TH1D *mcp_amp =new TH1D("mcp_amp","histomcp_ampl",nbinx,0.0,1);


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
	h2_l->Fill(amp_max[3]/(max*fit_l->GetParameter(1)),time[1+LEDi]-time[0]);
	h2_r->Fill(amp_max[4]/(max*fit_r->GetParameter(1)),time[2+LEDi]-time[0]);
	h2_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]);

      }//chiudo if

  }//chiudo for k

  
  
  TH1D* histotemp_l;
  TH1D* histotemp_r;
  TH1D* histotemp_t;
  
  for(k=0;k<nbinx;k++){
    
    histotemp_l=h2_l->ProjectionY("h2_lprojY",k,k);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",k,k);
  
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
    
    
    if(k%20==0) cout << k << " / " << nbinx << endl;
  }//chiudo for k
    
  
  TCanvas* wf_c =new TCanvas("wf_altro","Plot wf_altro",1800,1100);
  TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);
    
  TF1* hyp_r = new TF1("hyp_r","[0]+[2]*log(x+[1])",0.8,3);
  TF1* hyp_l = new TF1("hyp_l","[0]+[2]*log(x+[1])",0.8,3);

  TF1* hyp_t = new TF1("hyp_t","[1]*x**2+[2]*x+[0]",-0.1,0.65);
  
  gStyle->SetOptStat("");

  
  /* SetParameters*/
  hyp_l->SetParameter(0, 8.51);
  hyp_l->SetParameter(1, 5);
  hyp_l->SetParameter(2, 1.2);
  
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

  rymin_lc=rymin_l-hyp_l->Eval(1);
  rymax_lc=rymax_l-hyp_l->Eval(1);
  rymin_rc=rymin_r-hyp_r->Eval(1);
  rymax_rc=rymax_r-hyp_r->Eval(1);
  tymin_c=tymin;
  tymax_c=tymax;

  cout<<"_____________________________________--"<<rymin_l<<"  "<<hyp_l->Eval(1.0)<<"   "<< hyp_l->GetParameter(0)<<endl;
  
  TH2D* hc_l= new TH2D("hc_l", "histo hc_l",nbinx,rxmin,rxmax,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_r= new TH2D("hc_r", "histo hc_r",nbinx,rxmin,rxmax,nbiny,rymin_rc,rymax_rc);
  TH2D* hc_t= new TH2D("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  TH2D* hc_tdiff= new TH2D("hc_tdiff", "histo hc_tdiff",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  TH2D* hc_tl= new TH2D("hc_tl", "histo hc_tl",nbinx,0.8,3,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_tr= new TH2D("hc_tr", "histo hc_tr",nbinx,0.8,3,nbiny,rymin_lc,rymax_lc);
  TH2D* hc_tot= new TH2D("hc_tot", "histo hc_tot",nbinx,0.8,3,nbiny,rymin_lc,rymax_lc);
    
  
   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	if(time[1+LEDi]-time[2+LEDi]<99 && time[1+LEDi]-time[2+LEDi]>-99)hc_l->Fill(amp_max[3]/(max*fit_l->GetParameter(1)),time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1))));
	if(time[1+LEDi]-time[2+LEDi]<99 && time[1+LEDi]-time[2+LEDi]>-99)hc_r->Fill(amp_max[4]/(max*fit_r->GetParameter(1)),time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1))));
	if(time[1+LEDi]-time[2+LEDi]<99 && time[1+LEDi]-time[2+LEDi]>-99)hc_t->Fill(time[1+LEDi]-hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1)))-time[2+LEDi]+hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1))),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1)))+hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1))))/2);

      }

  }//chiudo for k
  
    
   for(k=0;k<nbinx;k++){
     histotemp_l=hc_l->ProjectionY("hc_lprojY",k,k);
     histotemp_r=hc_r->ProjectionY("hc_rprojY",k,k);
     histotemp_t=hc_t->ProjectionY("hc_tprojY",k,k);
     
     xt[k]=txmin +(Double_t)(txmax-txmin)*k/nbinx;
     yt[k]=histotemp_t->GetMean();
     RMS[2][k]= histotemp_t->GetMeanError();
     
     x_l[k]=rxmin+(Double_t)(rxmax-rxmin)*k/nbinx;
     y_l[k]=histotemp_l->GetMean();
     RMS[0][k]= histotemp_l->GetMeanError();
     
     x_r[k]=rxmin+(Double_t)(rxmax-rxmin)*k/nbinx;
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
   
   graph_tc->Fit("fit_tdiff");

   wf_c->cd(4);
    hc_l->GetYaxis()->SetTitle("t_left-t_MCP [ns]");
    hc_l->GetXaxis()->SetTitle("max.amplitude [mV]");
    hc_l->Draw("COLZ");
    graph_lc->SetMarkerStyle(8);
    graph_lc->SetMarkerSize(.5);
    graph_lc->Draw("P");
 
    for(k=0;k<nbinx;k++){
   
   
    TH1D* histotemp_t;

   
   
    histotemp_t=hc_tdiff->ProjectionY("hc_tprojY",k,k);
   
   
    yt[k]=histotemp_t->GetMean();
    RMS[2][k]= histotemp_t->GetMeanError();

    delete histotemp_t;

    
  }//chiudo for k

    TGraphErrors* graph_tcdiff = new TGraphErrors(nbinx-1,xt,yt,0,RMS[2]);
   


 
  
  wf_c->cd(5);

   hc_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
   hc_r->GetXaxis()->SetTitle("max.amplitude [mV]");
   hc_r->Draw("COLZ");
   graph_rc->SetMarkerStyle(8);
   graph_rc->SetMarkerSize(.5);
   graph_rc->Draw("P");
   
   for(k=0;k<digiTree->GetEntries();k++){
     
     digiTree->GetEntry(k);
     
     if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
       {
	 
	 if(time[1+LEDi]-time[2+LEDi]<99 && time[1+LEDi]-time[2+LEDi]>-99) hc_tdiff->Fill(time[1+LEDi]-hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1)))-time[2+LEDi]+hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1))),(time[1+LEDi]+time[2+LEDi])/2-fit_tdiff->Eval(time[1+LEDi]-time[2+LEDi])+fit_tdiff->GetParameter(0)-time[0]-(hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1)))+hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1))))/2);
	 
	 if(time[1+LEDi]-time[2+LEDi]<99 && time[1+LEDi]-time[2+LEDi]>-99) hc_tl->Fill(amp_max[3]/(max*fit_l->GetParameter(1)),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1)))+hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1))))/2);
	 
	 if(time[1+LEDi]-time[2+LEDi]<99 && time[1+LEDi]-time[2+LEDi]>-99) hc_tl->Fill(amp_max[4]/(max*fit_r->GetParameter(1)),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/(max*fit_r->GetParameter(1)))+hyp_l->Eval(amp_max[3]/(max*fit_l->GetParameter(1))))/2);
	 
	 
      }
     
   }//chiudo for k

   int npoints = 10;
   rxmin+=0.1;
   rxmax-=0.2;
   for (i=1;i<nbinx/npoints+1;i++){
     cut[i] =rxmin+(Double_t)(rxmax-rxmin)*(i*npoints)/nbinx;  
   }
   
   TCanvas *proiettati;
   TString histoname;
   TString title;
   TF1* fit;
   bool SaveProjection=false;
   Double_t XMinPlot=-0.3, XMaxPlot=0.3;

   if(SaveProjection){
     gSystem->Exec("rm -r -f HistoTResAmp");
     gSystem->Exec("mkdir HistoTResAmp");
   }
   
   for (i=1;i<nbinx/npoints+1;i++){
     fit = new TF1("fit","gaus",XMinPlot,XMaxPlot);
     fit->SetParameter(0,200);
     fit->SetParameter(1,0);
     fit->SetParameter(2,30e-2);
     proiettati = new TCanvas("canvasProjection","",1200,800);

     
     histotemp_t=hc_tl->ProjectionY("hc_totprojY",hc_tot->GetXaxis()->FindBin(cut[i]),hc_tot->GetXaxis()->FindBin(cut[i+1]));
     cout << "____________" << hc_tot->GetXaxis()->FindBin(cut[i])<<"_____" << hc_tot->GetXaxis()->FindBin(cut[i+1]) << endl;
     histotemp_t->Fit("fit","R0");
     
     if(SaveProjection){
       histoname="";
       histoname.Append("HistoTResAmp/histo");
       histoname.Append(to_string(i));
       histoname.Append(".png");
       
       title="";
       title.Append("Amp/MIP -> ");
       title.Append(to_string(cut[i]));
       title.Resize(title.Sizeof()-4);
       title.Append("-");
       title.Append(to_string(cut[i+1]));
       title.Resize(title.Sizeof()-4);
       
       //histotemp_t->GetXaxis()->SetLimits(XMinPlot,XMaxPlot);
       //fit->GetXaxis()->SetLimits(XMinPlot,XMaxPlot);
       gStyle->SetOptFit(0110);
       fit->GetXaxis()->SetTitle("t_{ave}-t_{MCP} [ns]");
       fit->GetYaxis()->SetTitle("counts");
       fit->SetTitle(title);
       fit->Draw();
       histotemp_t->Draw("SAME");
       //fit->DrawF1(XMinPlot,XMaxPlot,"SAME");
       
       //gPad->WaitPrimitive();
       proiettati->SaveAs(histoname);
     }
     
     sigma[i]=sqrt((fit->GetParameter(2))*(fit->GetParameter(2))-0.015*0.015);
     erry[i]=fit->GetParError(2);
     errx[i]= (rxmax-rxmin)*npoints/(2*nbinx);
     

     delete histotemp_t;
     delete fit;
     delete proiettati;
   }

   
				     
   TCanvas* c_TresAmp = new TCanvas("c_TresAmp","c_rest_plot",600,550);
   TGraphErrors* TResAmp = new TGraphErrors(nbinx/npoints-3,cut,sigma,errx,erry);
   TF1* fitramp = new TF1("fitramp","[0]+[1]/sqrt(x)",rxmin,rxmax);
   TLine* intercept = new TLine(0.82,0.03,1.8,0.03);
   intercept->SetLineColor(kBlue);
   gStyle->SetOptFit(1111);

   fitramp->SetParameter(0,0.0);
   //fitramp->SetParLimits(0,0.0,10);
   fitramp->SetParameter(1,100);
   TResAmp->Fit("fitramp","0RL");
   TResAmp ->GetXaxis()->SetTitle("amp/mip peak");
   TResAmp ->GetYaxis()->SetTitle("sigma_{t_{ave}}(ns)");
   TResAmp ->SetMarkerStyle(8);
   TResAmp ->SetMarkerSize(.8);

   TResAmp->GetYaxis()->SetRangeUser(0.025,0.044);
   TResAmp->GetXaxis()->SetLimits(0.82,1.8);
   TResAmp ->Draw("AP");
   fitramp->DrawF1(0,3,"same");
   intercept->Draw("same");

      
}
