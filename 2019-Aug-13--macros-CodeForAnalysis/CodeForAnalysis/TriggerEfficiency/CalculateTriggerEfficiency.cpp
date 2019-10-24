#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TKey.h>

#include <TH1D.h>

TList* GetQAListGRID(const TString str,const TString dir);
TList* GetResultListGRID(const TString str,const TString dir);

void CalculateTriggerEfficiency(){

  const TString dirnameData = "TriggerEff";

  TString mcsource = "AnalysisResults.root";

  TList *ResultListMC = GetResultListGRID(mcsource,dirnameData);

  TH1D* hGeneratedEvents = (TH1D*) ResultListMC->FindObject("hGeneratedEvents");
  TH1D* hTriggeredEvents = (TH1D*) ResultListMC->FindObject("hTriggeredEvents");
  TH1D* hSelectedEvents = (TH1D*) ResultListMC->FindObject("hSelectedEvents");

  Double_t Binning[] = {0,10,20,40,100};
  hGeneratedEvents = (TH1D*) hGeneratedEvents->Rebin(4, "hGeneratedEvents",Binning);
  hTriggeredEvents = (TH1D*) hTriggeredEvents->Rebin(4, "hTriggeredEvents",Binning);
  hSelectedEvents = (TH1D*) hSelectedEvents->Rebin(4, "hSelectedEvents",Binning);

  TH1D* hTriggerEfficiency = (TH1D*) hTriggeredEvents->Clone("hTriggerEfficiency");
  hTriggerEfficiency->Divide(hTriggeredEvents, hGeneratedEvents, 1, 1, "B");
  hTriggerEfficiency->SetMarkerStyle(20);
  
  TH1D* hEventSelectionEfficiency = (TH1D*) hSelectedEvents->Clone("hEventSelectionEfficiency");
  hTriggerEfficiency->Divide(hSelectedEvents, hGeneratedEvents, 1, 1, "B");
  hEventSelectionEfficiency->SetMarkerStyle(21);
  
  TCanvas* cEventEfficiencies = new TCanvas("cEventEfficiencies","Event efficienies",800,800);
  hTriggerEfficiency->DrawCopy("E");
  hEventSelectionEfficiency->DrawCopy("Esame");
  
  TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(hTriggerEfficiency,"Trigger efficiency","f");
  leg->AddEntry(hEventSelectionEfficiency,"Event selection efficiency","l");
  leg->Draw();
  
  
}

//_________________________________________________________________________
TList* GetResultListGRID(const TString str,const TString dir){

  TFile *f = TFile::Open(str.Data());
  if(!f || f->IsZombie()){
    printf("Could not read file %s\n",str.Data());
    return NULL ;
  }
  if(f->TestBit(TFile::kRecovered)){
    printf("File \"%s\" is corrupt!\n",str.Data());
  }
  TDirectory *d = f->GetDirectory(dir.Data());
  if(!d || d->IsZombie()){
    printf("Could not read file %s\n",dir.Data());
    return NULL ;
  }
  if(d->TestBit(TFile::kRecovered)){
    printf("File \"%s\" is corrupt!\n",dir.Data());
  }

  TKey *k;
  TIter next(d->GetListOfKeys());
  while ((k = dynamic_cast<TKey *>(next()))){TString s(k->GetName()); if(s.Contains("Results")) break;}
  if(!k){
    printf("Output container not found\n");
    f->Close(); delete f;
    return NULL;
  }
  TList *returnlist = dynamic_cast<TList *>(k->ReadObj());
  f->Close(); delete f;
  return returnlist;
}

//_________________________________________________________________________
TList* GetQAListGRID(const TString str,const TString dir){

  TFile *f = TFile::Open(str.Data());
  if(!f || f->IsZombie()){
    printf("Could not read file %s\n",str.Data());
    return NULL ;
  }
  if(f->TestBit(TFile::kRecovered)){
    printf("File \"%s\" is corrupt!\n",str.Data());
  }
  TDirectory *d = f->GetDirectory(dir.Data());
  if(!d || d->IsZombie()){
    printf("Could not read file %s\n",dir.Data());
    return NULL ;
  }
  if(d->TestBit(TFile::kRecovered)){
    printf("File \"%s\" is corrupt!\n",dir.Data());
  }

  TKey *k;
  TIter next(d->GetListOfKeys());
  while ((k = dynamic_cast<TKey *>(next()))){
    TString s(k->GetName());
    if(s.Contains("QA")) break;
  }
  if(!k){
    printf("Output container not found\n");
    f->Close(); delete f;
    return NULL;
  }
  TList *returnlist = dynamic_cast<TList *>(k->ReadObj());
  f->Close(); delete f;
  return returnlist;
}
