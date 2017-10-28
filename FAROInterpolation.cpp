#define NLine 44
#define NRow 93

void FAROInterpolation(){
  TString datapath = "$(MEG2SYS)/analyzer/macros/xec/survey/";
  TString filename = "XECFAROMPPC_well_fitted.root";
  TString filepath= datapath + filename;
  TFile* fin = new TFile(filepath.Data(),"READ");
  TTree* tin = (TTree*)fin->Get("faro");
  //
  Int_t MPPCchannel;
  Double_t MPPCXPos;
  Double_t MPPCYPos;
  Double_t MPPCZPos;
  Bool_t DataQual;
  tin->SetBranchAddress("channel",&MPPCchannel);
  tin->SetBranchAddress("XPos",&MPPCXPos);
  tin->SetBranchAddress("YPos",&MPPCYPos);
  tin->SetBranchAddress("ZPos",&MPPCZPos);
  tin->SetBranchAddress("DataQual",&DataQual);

  Int_t Nwf=tin->GetEntries();
  std::cout<<"total channel: "<<Nwf<<std::endl;
  int NMppcRow[NRow]={};
  for(int i=0;i<Nwf;i++ ){
	tin->GetEntry(i);
	Int_t Row = MPPCchannel/NLine;
	Int_t Column = MPPCchannel%NLine;
	NMppcRow[Row]+=1;
  }
  for(int Row=0; Row<NRow;Row++){
	std::cout<<"Row: "<<Row<<" MPPC: "<<NMppcRow[Row]<<std::endl;
  }
}
