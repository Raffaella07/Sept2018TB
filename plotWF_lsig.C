//to run with ranges [0.125;2]right, [0.13;2]left on 4.1
//to run with ranges [0.125;2]right, [0.13;2]left on 1.3



void plotWF_lsig(const char * filename){


  TFile*  file= TFile::Open(filename);
  //TTree * WFTree = (TTree*)file->Get("wf");
  TTree* digiTree = (TTree*)file->Get("digi");



  Float_t amp_max[54], time[54];
  int k,maxbin_l,maxbin_r,maxbin_t,LED300;
  Float_t rxmin,rxmax,rymin_l,rymax_l,rymin_r,rymax_r,tymin,tymax;
  bool debug=false;
  Double_t max=0,tmax=0;
  rxmin=0;
  rxmax=0.5;


  const Int_t  nbinx=100,nbiny=550;

  rymin_l=5;
  rymax_l=21;
  rymin_r=5;
  rymax_r=21;
  tymin=6;
  tymax=17;



  Float_t x_r[nbinx],y_r[nbiny], x_l[nbinx],y_l[nbiny],rmsy_l[nbiny],rmsy_r[nbiny];
  Float_t xt[nbinx],yt[nbinx],rmsyt[nbinx];


  TH1F *hr_amp =new TH1F("hr_amp","histos_ampr",nbinx,0.0,1);
  TH1F *hl_amp =new TH1F("hl_amp","histos_ampl",nbinx,0.0,1);

  TF1 *fit_r = new TF1("f_r","landau",0.10,1);
  TF1 *fit_l = new TF1("f_l","landau",0.10,1);


  digiTree->SetBranchAddress("amp_max",&amp_max);
  digiTree->SetBranchAddress("time",&time);
    digiTree->SetBranchAddress("LED300",&LED300);
  // digiTree->SetBranchAddress("LED30",&LED30);
  //digiTree->SetBranchAddress("LED50",&LED50);

  /*  for(k=0; k<digiTree->GetEntries(); k++){
    digiTree->GetEntry(k);
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    if(amp_max[3]>max) {max=amp_max[3];}
    if(amp_max[4]>max) {max=amp_max[4];}
    if(time[3]-time[4]>tmax && time[3]-time[4]<10) {tmax = time[3]-time[4];}

  }
*/

  max=4096;
  for(k=0;k<digiTree->GetEntries();k++){
    if (k%10000==0) cout<<"On entry " <<k<<endl;
    digiTree->GetEntry(k);

    hr_amp->Fill(amp_max[3]/max);
    hl_amp->Fill(amp_max[4]/max);
  }/*chiudo for */

  hr_amp->Scale(1/(hr_amp->Integral()));
  hl_amp->Scale(1/(hl_amp->Integral()));

  cout << tmax <<endl;
  cout<< max << endl;

  hr_amp->Fit("f_r","RQ0");
  hl_amp->Fit("f_l","RQ0");


  TH2F* h2_l= new TH2F("h2_l", "histo h2_l",nbinx,rxmin,rxmax,nbiny,rymin_l,rymax_l);
  TH2F* h2_r= new TH2F("h2_r", "histo h2_r",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);
  TH2F* h2_m= new TH2F("h2_m", "histo h2_m",nbinx,rxmin,rxmax,nbiny,rymin_r,rymax_r);

  for(k=0;k<digiTree->GetEntries();k++){

    digiTree->GetEntry(k);

    if (0.8*(fit_l->GetParameter(1)) < (amp_max[3]/max) && (amp_max[3]/max) < (3*fit_l->GetParameter(1)))
      {
	h2_l->Fill(amp_max[3]/max,time[1+LED300]-time[0]);
	h2_r->Fill(amp_max[4]/max,time[2+LED300]-time[0]);
	h2_m->Fill((amp_max[3]+amp_max[4])/(2*max),(time[1+LED300]+time[2+LED300])/2-time[0]);


	if(debug) cout << 0.8*fit_l->GetParameter(1) << " < " << amp_max[3]/max << " < " << 3*fit_l->GetParameter(1) << " ////  " << time[4]-time[0] <<endl;
      }

  }//chiudo for k


    // for(k=0;k<nbinx;k++){
    TH1D* histotemp_l;
    TH1D* histotemp_r;
    TH1D* histotemp_m;

    histotemp_l=h2_l->ProjectionY("h2_lprojY",0,nbinx);
    histotemp_r=h2_r->ProjectionY("h2_rprojY",0,nbinx);
    histotemp_m=h2_m->ProjectionY("h2_tprojY",0,nbinx);





  TCanvas* wf_c =new TCanvas("wf","Plot wf",600,550);
  // TGraphErrors* graph_r=new TGraphErrors(nbinx-1,x_r,y_r,0,rmsy_r);
  // TGraphErrors* graph_l=new TGraphErrors(nbinx-1,x_l,y_l,0,rmsy_l);
  // TGraphErrors* graph_t=new TGraphErrors(nbinx-1,xt,yt,0,rmsyt);
  TF1* g_r = new TF1("g_r","gaus",0,21);
  TF1* g_l = new TF1("g_l","gaus",0,21);
  TF1* g_m = new TF1("g_m","gaus",0,21);
  //  hyp_r->SetParameter(0,10);
  // hyp_l->SetParameter(0,10);
  histotemp_m->Fit("g_m","0");
  histotemp_l->Fit("g_l","0");
  histotemp_r->Fit("g_r","0");

  // hyp_r->SetParLimits(0,1,8);
  gStyle->SetOptStat("");
  gStyle->SetOptFit();
  histotemp_r->SetLineColor(kRed);
  histotemp_m->SetLineColor(kBlack);
  g_m->SetLineColor(kBlack);

 

  histotemp_m->GetXaxis()->SetTitle("t_ave-t_MCP");
  histotemp_m->GetYaxis()->SetTitle("counts");
  
  histotemp_m->Draw();
  // histotemp_r->Draw("same");
  // histotemp_l->Draw("same");
  // g_l->Draw("same");
  // g_r->Draw("same");
  g_m->Draw("same");

 

}

