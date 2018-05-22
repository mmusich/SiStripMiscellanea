#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
void CopyDir(TDirectory *source) {
   //copy all objects and subdirs of directory source as a subdir of the current directory
  source->ls();
  TDirectory *savdir = gDirectory;
  TDirectory *adir;

  if(!((TString)(source->GetName())).Contains("root")){
    adir = savdir->mkdir(source->GetName());
  } else {
    adir = savdir;
  } 
  adir->cd();
  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {

      if( ((TString)(key->GetName())).Contains("gainCalibrationTreeIsoMuon") ) continue;
      if( ((TString)(key->GetName())).Contains("gainCalibrationTreeStdBunch") ) continue;
      
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir);
      adir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      adir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write();
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

void CopyFile(const char *fname) {
  //Copy all objects and subdirs of file fname as a subdir of the current directory
  TDirectory *target = gDirectory;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("Cannot copy file: %s\n",fname);
    target->cd();
    return;
  }
  target->cd();
  CopyDir(f);
  delete f;
  target->cd();
}

void copyFiles() {
  TFile *f = new TFile("result.root","recreate");
  f->SetCompressionSettings(6);
  CopyFile("calibTree_307073_50.root");
  f->ls();
  delete f;
}
