#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include <iostream>

void treeCloner(TString fname){

  TFile *f_input  = TFile::Open(fname,"READ");
  
  TObjArray* tokens=fname.Tokenize("/");
  TString input_name=(((TObjString*)tokens->Last())->String());
 
  TString output_name = input_name.ReplaceAll(".root","");
  std::cout<<"the output file name will be:" << output_name << std::endl;

  TFile *f_output = new TFile("./slimmed_"+output_name+".root","recreate");

  f_input->cd();

  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(f_input->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {

      std::cout <<"entering directory:"<< key->GetName() << std::endl;
      
      if( ((TString)(key->GetName())).Contains("gainCalibrationTreeIsoMuon") ) continue;
      if( ((TString)(key->GetName())).Contains("gainCalibrationTreeStdBunch") ) continue;
      
      f_input->cd(key->GetName());
      TDirectory *savdir = gDirectory;

      std::cout<< "output directory is:" << savdir->GetName() << std::endl;

      TTree *tinput;
      if((TString)(savdir->GetName())=="anEff") {
	tinput = (TTree*)savdir->Get("traj");
      } else {
	tinput = (TTree*)savdir->Get("tree");
      }
      f_output->cd();
      TDirectory *adir;
      adir = gDirectory->mkdir(savdir->GetName());
      adir->cd();
      tinput->CloneTree(-1,"fast")->Write();
      adir->SaveSelf(kTRUE);
      std::cout << " the tree is written" << std::endl;
    }
  }

  f_output->Close();
  f_input->Close();
  delete f_input;
  delete f_output;

}
