1,27c1,4
< //////////////////////////////////////////////////////////
< // This class has been automatically generated on
< // Thu Jul 12 13:10:55 2018 by ROOT version 6.10/09
< // from TChain tree/
< //////////////////////////////////////////////////////////
< 
< #ifndef TreeMkrTemplate_V12_h
< #define TreeMkrTemplate_V12_h
< 
< #include <TROOT.h>
< #include <TChain.h>
< #include <TFile.h>
< 
< // Header file for the classes stored in the TTree if any.
< #include "vector"
< #include "vector"
< #include "vector"
< #include "vector"
< #include "vector"
< #include "vector"
< 
< class TreeMkrTemplate_V12 {
< public :
<    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
<    Int_t           fCurrent; //!current Tree number in a TChain
< 
< // Fixed size dimensions of array or collections stored in the TTree if any.
---
> // Include file for TreeMaker V12
> //
> // Declare leaf variables and provide a function to set the branch addresses
> // Generated via TChain::MakeClass("TreeMkrTemplate")
459,534c436
<    TreeMkrTemplate_V12(TTree *tree=0);
<    virtual ~TreeMkrTemplate_V12();
<    virtual Int_t    Cut(Long64_t entry);
<    virtual Int_t    GetEntry(Long64_t entry);
<    virtual Long64_t LoadTree(Long64_t entry);
<    virtual void     Init(TTree *tree);
<    virtual void     Loop();
<    virtual Bool_t   Notify();
<    virtual void     Show(Long64_t entry = -1);
< };
< 
< #endif
< 
< #ifdef TreeMkrTemplate_V12_cxx
< TreeMkrTemplate_V12::TreeMkrTemplate_V12(TTree *tree) : fChain(0) 
< {
< // if parameter tree is not specified (or zero), connect the file
< // used to generate this class and read the Tree.
<    if (tree == 0) {
< 
< #ifdef SINGLE_TREE
<       // The following code should be used if you want this class to access
<       // a single tree instead of a chain
<       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
<       if (!f || !f->IsOpen()) {
<          f = new TFile("Memory Directory");
<       }
<       f->GetObject("tree",tree);
< 
< #else // SINGLE_TREE
< 
<       // The following code should be used if you want this class to access a chain
<       // of trees.
<       TChain * chain = new TChain("tree","");
<       chain->Add("/nfs/data38/cms/mulholland/lpcTrees/Skims/Run2ProductionV12/tree_signal/tree_TTZToLLNuNu.root/tree");
<       tree = chain;
< #endif // SINGLE_TREE
< 
<    }
<    Init(tree);
< }
< 
< TreeMkrTemplate_V12::~TreeMkrTemplate_V12()
< {
<    if (!fChain) return;
<    delete fChain->GetCurrentFile();
< }
< 
< Int_t TreeMkrTemplate_V12::GetEntry(Long64_t entry)
< {
< // Read contents of entry.
<    if (!fChain) return 0;
<    return fChain->GetEntry(entry);
< }
< Long64_t TreeMkrTemplate_V12::LoadTree(Long64_t entry)
< {
< // Set the environment to read one entry
<    if (!fChain) return -5;
<    Long64_t centry = fChain->LoadTree(entry);
<    if (centry < 0) return centry;
<    if (fChain->GetTreeNumber() != fCurrent) {
<       fCurrent = fChain->GetTreeNumber();
<       Notify();
<    }
<    return centry;
< }
< 
< void TreeMkrTemplate_V12::Init(TTree *tree)
< {
<    // The Init() function is called when the selector needs to initialize
<    // a new tree or chain. Typically here the branch addresses and branch
<    // pointers of the tree will be set.
<    // It is normally not necessary to make changes to the generated
<    // code, but the routine can be extended by the user if needed.
<    // Init() will be called many times when running on PROOF
<    // (once per file to be processed).
---
>   void setBranchAddress(TChain* fChain) {
685a588
> 
687,689c590,592
<    if (!tree) return;
<    fChain = tree;
<    fCurrent = -1;
---
>    /* if (!tree) return; */
>    /* fChain = tree; */
>    /* fCurrent_ = -1; */
905,917c808
<    Notify();
< }
< 
< Bool_t TreeMkrTemplate_V12::Notify()
< {
<    // The Notify() function is called when a new file is opened. This
<    // can be either for a new TTree in a TChain or when when a new TTree
<    // is started when using PROOF. It is normally not necessary to make changes
<    // to the generated code, but the routine can be extended by the
<    // user if needed. The return value is currently not used.
< 
<    return kTRUE;
< }
---
>    /* Notify(); */
919,933c810
< void TreeMkrTemplate_V12::Show(Long64_t entry)
< {
< // Print contents of entry.
< // If entry is not specified, print current entry
<    if (!fChain) return;
<    fChain->Show(entry);
< }
< Int_t TreeMkrTemplate_V12::Cut(Long64_t entry)
< {
< // This function may be called from Loop.
< // returns  1 if entry is accepted.
< // returns -1 otherwise.
<    return 1;
< }
< #endif // #ifdef TreeMkrTemplate_V12_cxx
---
>   };
