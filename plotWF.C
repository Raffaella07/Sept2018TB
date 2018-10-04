//pioooo
//pippo
void plotWF(const char * filename,int i){
  
  TFile *  file= new TFile(filename);
  TTree * WFTree = (TTree*)file->Get("wf");
  TTree * digiTree = (TTree*)file->Get("digi");
 
  Float_t tot_WF[14336],wf_time[14336],digi_time[28];
  Float_t ch[14][1024], time[1024],startime[2];
  int wf_sample,step[14336];
  int j,k,LED300;

  
  WFTree->SetBranchAddress("WF_time",&wf_time);
  WFTree->SetBranchAddress("WF_val",&tot_WF);
  // WFTree->SetBranchAddress("Instance",&step);
  WFTree->SetBranchAddress("WF_samples",&wf_sample);
  digiTree->SetBranchAddress("time",&digi_time);
  digiTree->SetBranchAddress("LED",&LED300);
  int p;
  p=i;
  
  for(k=0;k<WFTree->GetEntries();k++){

    if (i==0) p=k;
    cout<<p<<endl;
    WFTree->GetEntry(p);
    digiTree->GetEntry(p);
    
    for(j=0;j<wf_sample;j++){
      step[j]= j;
      ch[(int)j/1024][j%1024]=tot_WF[j];
      if(j<1024) time[j]= wf_time[j];
    }//chiudo for j

    startime[0]=digi_time[1+LED300];
    startime[1]=digi_time[2+LED300];
    TCanvas* wf_c =new TCanvas("wf","Plot wf",1000,650);
    TGraph* graph[14];
    TLine* start[2];
    wf_c->Divide(7,2);
    
    for (j=0;j<14;j++){  
      graph[j] = new TGraph(wf_sample/14,wf_time,ch[j]);
      wf_c->cd(j+1);
      //graph[j]->GetXaxis()->SetTitleSize(10);
      graph[j]->GetXaxis()->SetTitle("time [ns]");
      //graph[j]->GetYaxis()->SetTitleSize(10);
      graph[j]->GetYaxis()->SetTitle("amplitude [mV]");
      graph[j]->Draw();
      if(j==1){
	start[0]=new TLine(startime[0],0,startime[0],500);
	start[0]->SetLineColor(kRed);
	//	start[0]->Draw("same");
      }
      if(j==2){
	start[1]=new TLine(startime[1],0,startime[1],500);
	start[1]->SetLineColor(kRed);
	//start[1]->Draw("same");
      }
    }//chiudo for j
    
    
    wf_c->Update();  
    if(i!=0) break;
    gPad->WaitPrimitive();
   
  }//chiudo for k
}

