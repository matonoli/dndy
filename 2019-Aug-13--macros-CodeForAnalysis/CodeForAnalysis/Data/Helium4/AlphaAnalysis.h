#ifndef AlphaAnalysis_h
#define AlphaAnalysis_h

class TH1;
class TH2;
class TTree;

//_____________________________________________________________________________________
class AlphaAnalysis : public TObject {

  
  public :
  AlphaAnalysis();
  virtual ~AlphaAnalysis();
  
  void SetFileNames (const char *InputFileName, const char *InputListName, const char *OutputFileName) {
    fInputFileName  = InputFileName;
    fInputListName  = InputListName;
    fOutputFileName = OutputFileName;
  }
  
  void GetInputTree();
  void SetTreeBranches();
  void CreateHistograms();
  void ExecuteAnalysis();
  void WriteOutputFile();
  
private:
  Double_t CalculateRapidity(Double_t y);
  Double_t CorrectPt(Double_t px, Double_t py);
  
  const char *fInputFileName;
  const char *fInputListName;
  const char *fOutputFileName;
  TTree *tree;
  Bool_t SaveToFile;
  
  Int_t magFieldSign;
  Double_t multPercentile_V0M;
  Double_t multPercentile_V0A;
  Double_t multPercentile_V0C;
  Double_t multPercentile_OnlineV0M;
  Double_t multPercentile_OnlineV0A;
  Double_t multPercentile_OnlineV0C;
  Double_t multPercentile_ADM;
  Double_t multPercentile_ADA;
  Double_t multPercentile_ADC;
  Double_t multPercentile_SPDClusters;
  Double_t multPercentile_SPDTracklets;
  Double_t multPercentile_RefMult05;
  Double_t multPercentile_RefMult08;
  Double_t multPercentile_CL1;
  Double_t multPercentile_ZNA;
  Int_t Ntrk_V0M;
  Int_t Ntrk_V0A;
  Int_t Ntrk_V0C;
  Int_t Ntrk_OnlineV0M;
  Int_t Ntrk_OnlineV0A;
  Int_t Ntrk_OnlineV0C;
  Int_t Ntrk_ADM;
  Int_t Ntrk_ADA;
  Int_t Ntrk_ADC;
  Int_t Ntrk_SPDClusters;
  Int_t Ntrk_SPDTracklets;
  Int_t Ntrk_RefMult05;
  Int_t Ntrk_RefMult08;
  Int_t nVertexContributors;
  Double_t xVertex;
  Double_t yVertex;
  Double_t zVertex;
  Double_t px;
  Double_t py;
  Double_t pz;
  Double_t TPCmomentum;
  Double_t TRDmomentum;
  Double_t integratedLength;
  Double_t timeOfFlight;
  Double_t beta;
  Double_t gamma;
  Double_t mass;
  Int_t trackID;
  Double_t eta;
  Double_t phi;
  Double_t theta;
  Double_t y;
  Int_t q;
  Double_t dcaxy;
  Double_t dcaz;
  Int_t nTPC_Clusters;
  Int_t nTRD_Clusters;
  Int_t nITS_Clusters;
  Int_t nTPC_FindableClusters;
  Int_t nTPC_CrossedRows;
  Int_t nTPC_Clusters_dEdx;
  Int_t HasPointOnITSLayer0;
  Int_t HasPointOnITSLayer1;
  Int_t HasPointOnITSLayer2;
  Int_t HasPointOnITSLayer3;
  Int_t HasPointOnITSLayer4;
  Int_t HasPointOnITSLayer5;
  Int_t HasSharedPointOnITSLayer0;
  Int_t HasSharedPointOnITSLayer1;
  Int_t HasSharedPointOnITSLayer2;
  Int_t HasSharedPointOnITSLayer3;
  Int_t HasSharedPointOnITSLayer4;
  Int_t HasSharedPointOnITSLayer5;
  Double_t chi2_TPC;
  Double_t chi2_NDF;
  Double_t chi2_ITS;
  Double_t ITSsignal;
  Double_t TPCsignal;
  Double_t TOFsignal;
  Double_t TRDsignal;
  Double_t HMPIDsignal;
  Double_t nSigmaITS_He4;
  Double_t nSigmaTPC_He4;
  Double_t nSigmaTOF_He4;
  Double_t nSigmaTRD_He4;
  Double_t nSigmaHMPID_He4;
  Double_t nSigmaITS_He3;
  Double_t nSigmaTPC_He3;
  Double_t nSigmaTOF_He3;
  Double_t nSigmaTRD_He3;
  Double_t nSigmaHMPID_He3;
  Double_t nSigmaITS_Trit;
  Double_t nSigmaTPC_Trit;
  Double_t nSigmaTOF_Trit;
  Double_t nSigmaTRD_Trit;
  Double_t nSigmaHMPID_Trit;
  Double_t nSigmaITS_Deut;
  Double_t nSigmaTPC_Deut;
  Double_t nSigmaTOF_Deut;
  Double_t nSigmaTRD_Deut;
  Double_t nSigmaHMPID_Deut;
  Double_t nSigmaITS_Prot;
  Double_t nSigmaTPC_Prot;
  Double_t nSigmaTOF_Prot;
  Double_t nSigmaTRD_Prot;
  Double_t nSigmaHMPID_Prot;
  Double_t nSigmaITS_Pion;
  Double_t nSigmaTPC_Pion;
  Double_t nSigmaTOF_Pion;
  Double_t nSigmaTRD_Pion;
  Double_t nSigmaHMPID_Pion;
  Double_t nSigmaITS_Kaon;
  Double_t nSigmaTPC_Kaon;
  Double_t nSigmaTOF_Kaon;
  Double_t nSigmaTRD_Kaon;
  Double_t nSigmaHMPID_Kaon;
  Double_t nSigmaITS_Elec;
  Double_t nSigmaTPC_Elec;
  Double_t nSigmaTOF_Elec;
  Double_t nSigmaTRD_Elec;
  Double_t nSigmaHMPID_Elec;
  
  TF1* PtCorr;
  Double_t NumberOfEvents;
  
  //User Histograms
  TH1D *fHistoRapidity;
  TH1D *fHistoMean;
  TH1D *fHistoWidth;
  
  TH1D *fHistoPtSpectrum;
  TH1D *fHistoAntiPtSpectrum;
  TH1D *fHistoRatioSpectra;

  ClassDef(AlphaAnalysis,1) //calculate the CL upper limit using the Feldman-Cousins method
};
#endif
