//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_tdiff(const char * filename, Double_t *sig, Double_t * ssig){


  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");



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



  const Int_t  nbinx=200,nbiny=150;

  Int_t i;
  Double_t sigma[50],erry[50],cut[50],errx[50];
  
  txmin=-0.3;
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
  TH1F *mcp_amp =new TH1F("mcp_amp","histomcp_ampl",nbinx,0.0,1);


  TF1 *fit_r = new TF1("f_r","landau",0.06,1);
  TF1 *fit_l = new TF1("f_l","landau",0.1,1);


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
    if(time[1+LEDi]-time[0]<15 && time[1+LEDi]-time[0]>0) {
      counter1++;
      Times1[counter1]=time[1+LEDi]-time[0];
      cout << Times1[k] << endl;}
    //  else{Times1[k]=0;}
    if(time[2+LEDi]-time[0]<15 && time[2+LEDi]-time[0]>0) {
      counter2++;
      Times2[counter2]=time[2+LEDi]-time[0];}
    // else{Times2[k]=0;}  
    if((time[1+LEDi]+time[2+LEDi])/2-time[0]<15 && (time[1+LEDi]+time[2+LEDi])/2-time[0]>-1) {
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

  rymin_l=mean1-1*rms1;
  rymax_l=mean1+0.8*rms1;
  rymin_r=mean2-1*rms2;
  rymax_r=mean2+0.8*rms2;
    
  

  tymin=mean3-1*rms3;
  tymax=mean3+0.8*rms3;
 


  
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


  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2F* h2_t= new TH2F("h2_t", "histo h2_t",nbinx,txmin,txmax,nbiny,tymin,tymax);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);


    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	h2_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]);
	if (amp_max[4]/max < 0.35)h2_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]);
	h2_t->Fill((time[1+LEDi]-time[2+LEDi]),(time[1+LEDi]+time[2+LEDi])/2-time[0]);

      }//chiudo if
	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4+LEDi]-time[0] <<endl;
      

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

  TF1* hyp_r = new TF1("hyp_r","[0]+[2]*log(x+[1])",0.8*fit_r->GetParameter(1),1.5*fit_r->GetParameter(1));
  TF1* hyp_l = new TF1("hyp_l","[0]+[2]*log(x+[1])",0.5*fit_l->GetParameter(1),1.2*fit_l->GetParameter(1));

  TF1* hyp_t = new TF1("hyp_t","[1]*x**2+[2]*x+[0]",-0.1,0.65);
  
  gStyle->SetOptStat("");


    /* SetParameters*/
  hyp_l->SetParameter(0, 8.51);
  hyp_l->SetParameter(1, 5);
  hyp_l->SetParameter(2, 1.2);
  /* hyp_l->SetParameter(3, -2.43e-2);
  */
  hyp_r->SetParameter(0, 8.51);
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
  graph_l->Draw("SAMEP");
  hyp_l->Draw("same");


  wf_c->cd(2);
  h2_r->GetYaxis()->SetTitle("t_right-t_MCP [ns]");
  h2_r->GetXaxis()->SetTitle("max.amplitude [mV]");
  h2_r->Draw("COLZ");
  graph_r->Fit("hyp_r","0L");
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

  rymin_lc=rymin_l-hyp_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))+hyp_l->GetParameter(0);
  rymax_lc=rymax_l-hyp_l->Eval(fit_l->GetParameter(1)+0.5*fit_l->GetParameter(2))+hyp_l->GetParameter(0);
rymin_rc=rymin_r-hyp_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))+hyp_r->GetParameter(0);
rymax_rc=rymax_r-hyp_r->Eval(fit_r->GetParameter(1)+0.5*fit_r->GetParameter(2))+hyp_r->GetParameter(0);
  tymin_c=tymin-(hyp_l->Eval(0.25)-hyp_l->GetParameter(0)+hyp_r->Eval(0.25)-hyp_r->GetParameter(0))/2;
  tymax_c=tymax-(hyp_l->Eval(0.25)-hyp_l->GetParameter(0)+hyp_r->Eval(0.25)-hyp_r->GetParameter(0))/2;

  
  TH2F* hc_l= new TH2F("hc_l", "histo hc_l",nbinx,rxmin,rxmax,nbiny,rymin_lc,rymax_lc);
  TH2F* hc_r= new TH2F("hc_r", "histo hc_r",nbinx,rxmin,rxmax,nbiny,rymin_rc,rymax_rc);
  TH2F* hc_t= new TH2F("hc_t", "histo hc_t",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  TH2F* hc_tdiff= new TH2F("hc_tdiff", "histo hc_tdiff",nbinx,txmin,txmax,nbiny,tymin_c,tymax_c);
  

  
   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	hc_l->Fill(amp_max[3]/max,time[1+LEDi]-time[0]-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0));
        hc_r->Fill(amp_max[4]/max,time[2+LEDi]-time[0]-hyp_r->Eval(amp_max[4]/max)+hyp_r->GetParameter(0));
	hc_t->Fill(time[1+LEDi]-time[2+LEDi]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[3]/max)-hyp_l->GetParameter(0))/2);

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
   TF1* fit_tdiff = new TF1("fit_tdiff","[0]+[1]*x",-0.1,0.65);
  
   graph_tc->Fit("fit_tdiff");
   
   for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);
    if (0.8*(fit_l->GetParameter(1)) < (amp_max[4]/max) && (amp_max[4]/max) < (3*fit_l->GetParameter(1)) && amp_max[0]/max > 0.4 && amp_max[0]/max < 0.75)
      {
	hc_tdiff->Fill(time[1+LEDi]-time[2+LEDi]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)),(time[1+LEDi]+time[2+LEDi])/2-time[0]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)+hyp_l->Eval(amp_max[3]/max)-hyp_l->GetParameter(0))/2-fit_tdiff->Eval(time[1+LEDi]-time[2+LEDi]-(hyp_r->Eval(amp_max[4]/max)-hyp_r->GetParameter(0)-hyp_l->Eval(amp_max[3]/max)+hyp_l->GetParameter(0)))+fit_tdiff->GetParameter(0)); //GIUSTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
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

   TCanvas * timeres = new TCanvas("timeres","plot_timeres",600,550);

   TLegend* l1=new TLegend(0.1,0.7,0.48,0.9);
   l1->SetHeader("time stamps");
   l1->AddEntry(histo_cl,"t_left-t_MCP");
   l1->AddEntry(histo_cr,"t_right-t_MCP");
   l1->AddEntry(histo_cr,"t_ave-t_MCP");
   
   gStyle->SetOptStat("");
   histo_ct->Draw();
   gaus_ct->SetParameter(0,500);
   histo_ct->GetYaxis()->SetTitle("counts");
   histo_ct->GetXaxis()->SetTitle("t_{ave}-t_{MCP}(ns)");
   histo_cl->SetLineColor(kBlue);

   histo_cl->Fit("gaus_cl");
   histo_cr->Fit("gaus_cr");
   if (blind==true) histo_ct->Fit("gaus_ct");

   histo_cr->Draw("same");
  
   histo_cl->Draw("same");

   l1->Draw();

   TCanvas* tdiff = new TCanvas("tdiff","plot_tdiff",600,550);
   TLegend* l2=new TLegend(0.1,0.7,0.48,0.9);
   TH1D* histotemp_t[(Int_t)nbinx/13];
   TF1* fit[(Int_t)nbinx/13];
   l2->SetHeader("time stamps");
   l2->AddEntry(histo_ct,"t_ave-t_MCP");
   l2->AddEntry(histo_ctdiff,"t_ave-t_MCP(tdiff corr)");
   gaus_ctdiff->SetLineColor(kGreen);
   histo_ctdiff->Fit("gaus_ctdiff");
   histo_ctdiff->SetLineColor(kGreen);
   histo_ct->Draw();
   histo_ctdiff->Draw("same");
   l2->Draw();

   *sig = gaus_ctdiff->GetParameter(2);
   *ssig = gaus_ctdiff->GetParError(2);
  
   cout << "########################### "<< gaus_ct->GetParameter(2)/gaus_ct->GetParameter(1) << "_________" << gaus_ctdiff->GetParameter(2)/gaus_ctdiff->GetParameter(1) << endl;
   TCanvas* rest_gaussine = new TCanvas("rest_gaussine","rest_plotgaus",1800,1100);
   rest_gaussine->Divide(nbinx/(13*4),4);
   for (i=0;i<nbinx/13;i++){
     cut[i] =txmin+(Float_t)(txmax-txmin)*(i*13)/nbinx;
   }
   for (i=0;i<nbinx/13;i++){
     rest_gaussine->cd(i+1);
     fit[i] = new TF1(((string)("fit"+to_string(i))).c_str(),"gaus",0,10);
    
     histotemp_t[i]=hc_tdiff->ProjectionY(((string)("histoY"+to_string(i))).c_str(),hc_tdiff->GetXaxis()->FindBin(cut[i]-(cut[i+1]-cut[i])/2),hc_tdiff->GetXaxis()->FindBin(cut[i]+(cut[i+1]-cut[i])/2));
     cout << "____________" << hc_tdiff->GetXaxis()->FindBin(cut[i]-(cut[i+1]-cut[i])/2) << endl;
     gStyle->SetOptFit(00010);
     histotemp_t[i]->Fit(("fit"+to_string(i)).c_str(),"R");
     histotemp_t[i]->Draw();
     
     sigma[i]=sqrt(fit[i]->GetParameter(2)*fit[i]->GetParameter(2)-0.015*0.015);
     erry[i]=fit[i]->GetParError(2);
     errx[i]= (txmax-txmin)*13/(2*nbinx);

   
   
     // delete fit;


   }

   histotemp_t[4]->Draw();
   TCanvas* rest_plot = new TCanvas("rest","rest_plot",600,550);
   TGraphErrors* rest = new TGraphErrors(nbinx/13,cut,sigma,errx,erry);
  
   rest->GetXaxis()->SetTitle("t_{left}-t_{right}(ns)");
   rest->GetYaxis()->SetTitle("#sigma_{t_{ave}}(ns)");
   rest->SetMarkerStyle(8);
   rest->SetMarkerSize(.8);
   rest->Draw("AP");
   


   bool ConfrontoTdiff=true;
   
   if(ConfrontoTdiff){
     TCanvas* ConfCanv = new TCanvas("confronto","",1200,800);
     
     ConfCanv->Divide(3,1);
     
     ConfCanv->cd(1);
     h2_t->Draw("COLZ");
     graph_t->Draw("P");

     ConfCanv->cd(2);
     hc_t->Draw("COLZ");
     graph_tc->Draw("P");

     ConfCanv->cd(3);
     hc_tdiff->Draw("COLZ");
     graph_tcdiff->Draw("P");

   }

   
}


void ResConf(){
  string filename,temp;
  Double_t conf[5];
  Int_t i,n;
  Double_t sigma[5], ssigma[5];

  conf[0]=1.1;
  conf[1]= 1.2;
  conf[2]= 1.3;
  conf[3]=4.1;
  conf[4]=4.2;
  cout << "hello!" << endl;
 
 
  
  for(i=0;i<5;i++){
    temp = to_string(conf[i]);
    temp.resize(3);
    cout << temp << endl;
    filename =" Pd/ConfT100-B72-"+ temp+".root";
    plotWF_tdiff(filename.c_str(),&sigma[i],&ssigma[i]);
	

}
  TGraphErrors* gconf = new TGraphErrors(5,conf,sigma,0,ssigma);
  TCanvas* plot_conf = new TCanvas("plot_conf","plot_conf",600,550);
  gconf->SetMarkerStyle(8);
  gconf->SetMarkerSize(.5);
  gconf->Draw("AP");
}
