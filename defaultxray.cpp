#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#ifndef __CINT__
#include "ROMETreeInfo.h"
#endif

void xray_ana(Int_t run, Int_t scan);
Double_t XYZ2Phi(Double_t x, Double_t y, Double_t /* z */);


/*
Analysis of X-ray meas.
xra_ana(Int_t run, Int_t scan)
run : run number
scan : drection of scan. 0: Phi scan, 1 : Z scan
*/

void defaultxray(int runnum, int scan) {
   gStyle->SetFuncColor(kRed);
xray_ana(runnum, scan);

  return;
}

//________________________________________________________________________________
void xray_ana(Int_t run, Int_t scan) {

   /*-----Initialization----*/
   const Int_t nMPPC=4092;
   const Int_t kNchforXray = 272;

   Int_t valid[kNchforXray] ={};
   Int_t MPPCindex[kNchforXray] ={};
   Float_t MPPCXYZ[kNchforXray][3] ={};
   Float_t MPPCPhi[kNchforXray]={};

   gStyle->SetTitleW(2);
   TGraph* grBeamPosition= new TGraph();
   TGraph* grScaler[kNchforXray];
   for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
      grScaler[iScalerCh] = new TGraph();
   }

   TFile *fout =new TFile("xray.root","RECREATE");

   /*-----Define rec tree to be read----*/
   TFile *frec = new TFile(Form("$(MEG2SYS)/analyzer/rec%06d.root", run),"READ");
   TTree *rec = (TTree*)frec->Get("rec");
   Int_t nEvent = rec->GetEntries();
   std::cout<<"Read rec "<<run<<". "<<nEvent<<" scaler events."<<std::endl;

   /*-----Read PM RunHeader----*/
   TClonesArray* pmrhArray = (TClonesArray*)frec->Get("XECPMRunHeader");

   MEGXECPMRunHeader *pmrh = 0;
   Int_t ScalerCh =0;
   for (Int_t iPM = 0; iPM < nMPPC; iPM++) {
      pmrh = (MEGXECPMRunHeader*)(pmrhArray->At(iPM));
      if (pmrh->GetDRSAddress() >= 0) {
         // std::cout<<pmrh->GetDRSChipID()<<" "<<pmrh->GetDRSChannelOnChip()<<std::endl;
         ScalerCh = pmrh->GetDRSChipID() * 8 + pmrh->GetDRSChannelOnChip();
         if (ScalerCh < 0 || ScalerCh > kNchforXray) continue;

         valid[ScalerCh] = 1;
         MPPCindex[ScalerCh] = iPM;
         MPPCXYZ[ScalerCh][0] = pmrh->GetXYZAt(0) * 10;
         MPPCXYZ[ScalerCh][1] = pmrh->GetXYZAt(1) * 10;
         MPPCXYZ[ScalerCh][2] = pmrh->GetXYZAt(2) * 10;
         MPPCPhi[ScalerCh] = XYZ2Phi(MPPCXYZ[ScalerCh][0], MPPCXYZ[ScalerCh][1], MPPCXYZ[ScalerCh][2]);
      }
   }

   /*-----Define rec tree to be read----*/
   TClonesArray* recTRGScaler = new TClonesArray("MEGTRGScaler");
   TBranch *branchrectrgscaler = 0;
   branchrectrgscaler = rec->GetBranch("trgscaler.");
   if (branchrectrgscaler) branchrectrgscaler->SetAddress(&recTRGScaler);

   MEGXRAYData* recXRAYData = new MEGXRAYData();
   TBranch *branchrecxraydata = 0;
   branchrecxraydata = rec->GetBranch("xraydata.");
   if (branchrecxraydata) branchrecxraydata->SetAddress(&recXRAYData);

   /*-----Read each event----*/
   Int_t Scaler_buf;
   Float_t BeamZ_buf, BeamPhi_buf;
   for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {
      branchrectrgscaler->GetEntry(iEvent);
      branchrecxraydata->GetEntry(iEvent);

      if (recXRAYData->GetIsGood()) continue;
      BeamZ_buf = recXRAYData->GetBeamZ();
      BeamPhi_buf = recXRAYData->GetBeamPhi();
      // std::cout<<BeamZ_buf<<" "<<BeamPhi_buf<<std::endl;
      grBeamPosition->SetPoint(grBeamPosition->GetN(), BeamZ_buf, BeamPhi_buf);

      Int_t nTRGScalerEntries = recTRGScaler->GetEntriesFast();
      if (nTRGScalerEntries < kNchforXray) return;
      // Int_t iPnt = 0;
      for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
         Scaler_buf = ((MEGTRGScaler*)(recTRGScaler->At(iScalerCh)))->GetScaler();
         // grScaler[iScalerCh]->SetPoint(iEvent, iEvent, Scaler_buf);
         if (scan == 0) {
            if (TMath::Abs(BeamZ_buf - MPPCXYZ[iScalerCh][2]) < 18. ){
               grScaler[iScalerCh]->SetPoint(grScaler[iScalerCh]->GetN(), BeamPhi_buf, Scaler_buf);
            }
         } else if (scan == 1) {
            if (TMath::Abs(BeamPhi_buf - MPPCPhi[iScalerCh]) < 1.5 && TMath::Abs(MPPCXYZ[iScalerCh][2])< 157 ){
               grScaler[iScalerCh]->SetPoint(grScaler[iScalerCh]->GetN(), BeamZ_buf, Scaler_buf);
            }
         }
      }
   }

   /*-----Draw----*/
   // TCanvas *cgr0=new TCanvas("cgr0", "cgr0");
   // grScaler[0]->Draw();

   /*-----Save to file----*/
   fout->cd();
   grBeamPosition->Write("BeamPosition");

   for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
      if (!valid[iScalerCh]) {
         grScaler[iScalerCh]->SetTitle("No MPPC connected.");
      } else {
         grScaler[iScalerCh]->SetTitle(Form("WD: %d, MPPC: %04d (@phi: %2.1lf deg, z: %2.1lf mm)", iScalerCh, MPPCindex[iScalerCh], MPPCPhi[iScalerCh], MPPCXYZ[iScalerCh][2]));
      }
   }

   TCanvas* c1 = new TCanvas("c1", "c1", 0, 0, 1200, 1000);
   c1->Divide(5,4);
   
   int iCh = 0;
   for (Int_t iScalerCh = 0; iScalerCh < kNchforXray; iScalerCh++) {
      if (valid[iScalerCh] && grScaler[iScalerCh]->GetN() >10){
         grScaler[iScalerCh]->Write(Form("WDch%03d",iScalerCh));
         iCh++;
         c1->cd(iCh);
         grScaler[iScalerCh]->Draw("APL");
      }
   }
   if (iCh > 20) {
      std::cout << "Too many channels to draw. Please use TBrowser." << std::endl;
   }

   return;
}


//________________________________________________________________________________
Double_t XYZ2Phi(Double_t x, Double_t y, Double_t /* z */) {
   return x == 0.0 && y == 0.0 ? 0.0 : 180 * degree + TMath::ATan2(-y, -x) * radian;
}
