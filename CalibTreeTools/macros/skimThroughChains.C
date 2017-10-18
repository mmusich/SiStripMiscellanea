void skimThroughChains(UInt_t run_number=283946) {
  
  gROOT->Reset();

  // See how long the regression takes.
  TStopwatch timer;
  timer.Start(kTRUE);

  //Get old file, old tree and set top branch address
  // an eff tree
  TChain oldtree("anEff/traj");
  TSystemDirectory dir("root://eoscms.cern.ch//store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR16_ReReco_old/","root://eoscms.cern.ch//store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR16_ReReco_old/");
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(Form("%i",run_number))) {
  	cout <<" adding "<< fname.Data() << endl;
  	oldtree.Add((dir.GetName()+fname).Data());
      }
    }
    delete file;
  }

  //oldtree.Add("root://eoscms.cern.ch//store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR16_ReReco_old/calibTree_283946_25.root");
  //oldtree.Add("root://eoscms.cern.ch//store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR16_ReReco_old/calibTree_283946_50.root");

  Int_t nentries = (Int_t)oldtree.GetEntries();

  std::cout<<"there are: "<<nentries << " entries in total in the HitEff tree"<< std::endl;
  UInt_t run;
  oldtree.SetBranchAddress("run",&run);
  
  //Create a new file + a clone of old tree in new file
  TFile *newfile = new TFile(("skimmed_"+std::to_string(run_number)+".root").c_str(),"recreate");
  TTree *newtree = oldtree.CloneTree(0);

  for (Int_t i=0;i<nentries; i++) {
    oldtree.GetEntry(i);
    for(Int_t j=0;j<100;j++){
      if(i%nentries==j*(nentries/100))std::cout<<i<<" |fraction of hit eff tree: "<< j  <<" % completed" << std::endl;
    }
    if (run==run_number) newtree->Fill();
  }
  newtree->Print();
  newtree->AutoSave();

  TChain oldtree2("testTree/tree");
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root") && fname.Contains(Form("%i",run_number))) {
  	cout <<" adding "<< fname.Data() << endl;
  	oldtree2.Add((dir.GetName()+fname).Data());
      }
    }
    delete file;
  }
  
  //oldtree2.Add("root://eoscms.cern.ch//store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR16_ReReco_old/calibTree_283946_25.root");
  //oldtree2.Add("root://eoscms.cern.ch//store/group/dpg_tracker_strip/comm_tracker/Strip/Calibration/calibrationtree/GR16_ReReco_old/calibTree_283946_50.root");

  Int_t nentries2 = (Int_t)oldtree2.GetEntries();

  // lumi info tree
  std::cout<<"there are: "<<nentries2 << " entries in total in the lumi info tree"<< std::endl;
  UInt_t run2;
  oldtree2.SetBranchAddress("run",&run2);
  
  //Create a new file + a clone of old tree in new file
  TTree *newtree2 = oldtree2.CloneTree(0);

  for (Int_t k=0;k<nentries2; k++) {
    oldtree2.GetEntry(k);
    for(Int_t l=0;l<100;l++){
      if(k%nentries==l*(nentries2/100))std::cout<<k<<" |fraction of lumi tree : "<< l  <<" % completed" << std::endl;
    }
    if (run==run_number) newtree2->Fill();
  }
  newtree2->Print();
  newtree2->AutoSave();

  delete newfile;

  std::cout << std::endl << "Total calculation time: " << timer.RealTime() << std::endl;

}
