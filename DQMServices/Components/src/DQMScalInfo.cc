/*
 * \file DQMDcsInfo.cc
 * \author A.Meyer - DESY
 * Last Update:
 *
 */

#include "DQMScalInfo.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

using namespace std;

const static int XBINS=100;
const static int YBINS=23;

// Framework

DQMScalInfo::DQMScalInfo(const edm::ParameterSet& ps)
{
  std::cout << "Parameter setting" << std::endl;
  parameters_ = ps;

  scalfolder_          = parameters_.getUntrackedParameter<std::string>("dqmScalFolder", "Scal") ;
  gtCollection_        = consumes<L1GlobalTriggerReadoutRecord>(parameters_.getUntrackedParameter<std::string>("gtCollection","gtDigis"));
  dcsStatusCollection_ = consumes<DcsStatusCollection>(parameters_.getUntrackedParameter<std::string>("dcsStatusCollection","scalersRawToDigi"));
  l1tscollectionToken_ = consumes<Level1TriggerScalersCollection>(parameters_.getUntrackedParameter<std::string>("l1TSCollection", "scalersRawToDigi"));

  for (int i=0;i<25;i++) dcs25[i]=true;
  lastlumi_=0;
}

DQMScalInfo::~DQMScalInfo(){
}

void DQMScalInfo::bookHistograms(DQMStore::IBooker & ibooker,
                                edm::Run const & /* iRun */,
                                edm::EventSetup const & /* iSetup */) {
  std::cout << "histo booking" << std::endl;
  // Fetch GlobalTag information and fill the string/ME.
  ibooker.cd();
  ibooker.setCurrentFolder(scalfolder_ +"/L1TriggerScalers/");
  const int fracLS = 16;
  const int maxLS  = 250;
  const int maxNbins = 2001;
  //  MySCALhisto_ = ibooker.book1D("MyScalHisto","MyScalHisto",10,0.,10.);
  hlresync_    = ibooker.book1D("lresync","Orbit of last resync",fracLS*maxLS,0,maxLS*262144);
  hlOC0_       = ibooker.book1D("lOC0","Orbit of last OC0",fracLS*maxLS,0,maxLS*262144);
  hlTE_        = ibooker.book1D("lTE","Orbit of last TestEnable",fracLS*maxLS,0,maxLS*262144);
  hlstart_     = ibooker.book1D("lstart","Orbit of last Start",fracLS*maxLS,0,maxLS*262144);
  hlEC0_       = ibooker.book1D("lEC0","Orbit of last EC0",fracLS*maxLS,0,maxLS*262144);
  hlHR_        = ibooker.book1D("lHR","Orbit of last HardReset",fracLS*maxLS,0,maxLS*262144);

  hphysTrig_   = ibooker.book1D("Physics_Triggers", "Physics Triggers", maxNbins, -0.5, double(maxNbins)-0.5);
  hphysTrig_->setAxisTitle("Lumi Section", 1);

  ibooker.cd();
  ibooker.setCurrentFolder(scalfolder_ +"/DcsStatus/");
  //  DCS_dummy_ = ibooker.book1D("DCS_dummy","DCS_dummy",10,0,10);
  DCSbyLS_     = ibooker.book2D("DCSbyLS", "DCS status vs Lumi", XBINS, 1., XBINS+1, YBINS+1, 0., YBINS+1);
  DCSbyLS_->setBinLabel(1,"CSC+",2);
  DCSbyLS_->setBinLabel(2,"CSC-",2);
  DCSbyLS_->setBinLabel(3,"DT0",2);
  DCSbyLS_->setBinLabel(4,"DT+",2);
  DCSbyLS_->setBinLabel(5,"DT-",2);
  DCSbyLS_->setBinLabel(6,"EB+",2);
  DCSbyLS_->setBinLabel(7,"EB-",2);
  DCSbyLS_->setBinLabel(8,"EE+",2);
  DCSbyLS_->setBinLabel(9,"EE-",2);
  DCSbyLS_->setBinLabel(10,"ES+",2);
  DCSbyLS_->setBinLabel(11,"ES-",2);
  DCSbyLS_->setBinLabel(12,"HBHEa",2);
  DCSbyLS_->setBinLabel(13,"HBHEb",2);
  DCSbyLS_->setBinLabel(14,"HBHEc",2);
  DCSbyLS_->setBinLabel(15,"HF",2);
  DCSbyLS_->setBinLabel(16,"HO",2);
  DCSbyLS_->setBinLabel(17,"BPIX",2);
  DCSbyLS_->setBinLabel(18,"FPIX",2);
  DCSbyLS_->setBinLabel(19,"RPC",2);
  DCSbyLS_->setBinLabel(20,"TIBTID",2);
  DCSbyLS_->setBinLabel(21,"TOB",2);
  DCSbyLS_->setBinLabel(22,"TECp",2);
  DCSbyLS_->setBinLabel(23,"TECm",2);
  DCSbyLS_->setBinLabel(24,"CASTOR",2);
  //  DCSbyLS_->setBinLabel(25,"ZDC",2);
  DCSbyLS_->setAxisTitle("Luminosity Section");
  DCSbyLS_->getTH2F()->SetCanExtend(TH1::kAllAxes);

  for (int i=0;i<25;i++) dcs25[i]=true;
  lastlumi_=0;

  ibooker.cd();
  ibooker.setCurrentFolder(scalfolder_ +"/BeamSpotOnline/");
  BeamSpot_dummy_ = ibooker.book1D("BeamSpot_dummy","BeamSpot_dummy",10,0.,10.);

  ibooker.cd();
  ibooker.setCurrentFolder(scalfolder_ +"/L1AcceptBunchCrossing/");
  L1AcceptBunchCrossing_dummy_ = ibooker.book1D("L1AcceptBunchCrossing_dummy","L1AcceptBunchCrossing_dummy",10,0.,10.);

  ibooker.cd();
  ibooker.setCurrentFolder(scalfolder_ +"/LumiScalers/");
  LumiScale_dummy_ = ibooker.book1D("LumiScale_dummy","LumiScale_dummy",10,0.,10.);

}

void DQMScalInfo::analyze(const edm::Event& e, const edm::EventSetup& c){
  makeL1Scalars(e);
  makeDcsInfo(e);
  //  makeGtInfo(e);
  return;
}

void
DQMScalInfo::endLuminosityBlock(const edm::LuminosityBlock& l, const edm::EventSetup& c)
{
  std::cout << "Ending LS blocks" << std::endl;

  int nlumi = l.id().luminosityBlock();
  cout << "nlumi = " << nlumi << endl;

  if (nlumi <= lastlumi_ ) return;

  TH2F * h_DCSbyLS_ = DCSbyLS_->getTH2F();
  // set to previous in case there was a jump or no previous fill
  for (int l=lastlumi_+1;l<nlumi;l++){
    // setting valid flag to zero for missed LSs
    //  DCSbyLS_->setBinContent(l,YBINS+1,0.);
    h_DCSbyLS_->SetBinContent(l,YBINS+1,0.);
    //    h_DCSbyLS_->SetFillColor(6);

    // setting all other bins to -1 for missed LSs
    for (int i=0;i<YBINS;i++){
      //      DCSbyLS_->setBinContent(l,i+1,-1.);
      h_DCSbyLS_->SetBinContent(l,i+1,-1);
      //      h_DCSbyLS_->SetFillColor(9);
    }
  }

  // fill dcs vs lumi
  //  DCSbyLS_->setBinContent(nlumi,YBINS+1,1.);
  h_DCSbyLS_->SetBinContent(nlumi,YBINS+1,1.);
  //  h_DCSbyLS_->SetFillColor(7);
  for(int i=0; i<25;i++){
    if(dcs25[i]){
      //       DCSbyLS_->setBinContent(nlumi,i+1,1.);
      h_DCSbyLS_->SetBinContent(nlumi,i+1,1.);
      //      h_DCSbyLS_->SetFillColor(4);
    } else {
      //       DCSbyLS_->setBinContent(nlumi,i+1,0.);
      h_DCSbyLS_->SetBinContent(nlumi,i+1,0.);
      //      h_DCSbyLS_->SetFillColor(2);
    }
    // set next lumi to -1 for better visibility
    if (nlumi < XBINS)
      //      DCSbyLS_->setBinContent(nlumi+1,i+1,-1.);
      h_DCSbyLS_->SetBinContent(nlumi+1,i+1,-1.);  
    
    dcs25[i]=true;
  }

  //Set info to be plotted vs LS
  //  MySCALhisto_->setBinContent(1,1.);
  //  DCS_dummy_->setBinContent(1,1.);


  BeamSpot_dummy_->setBinContent(1,1.);
  L1AcceptBunchCrossing_dummy_->setBinContent(1,1.);
  LumiScale_dummy_->setBinContent(1,1.);

  lastlumi_=nlumi;
  // cout << "lastLimi = " << lastlumi_ << endl;

  return;
}

void
DQMScalInfo::makeL1Scalars(const edm::Event& e)
{
  edm::Handle<Level1TriggerScalersCollection> l1ts;
  e.getByToken(l1tscollectionToken_,l1ts);
  if(l1ts->size()==0) return;
  hlresync_->Fill((*l1ts)[0].lastResync());
  hlOC0_->Fill((*l1ts)[0].lastOrbitCounter0());
  hlTE_->Fill((*l1ts)[0].lastTestEnable());
  hlstart_->Fill((*l1ts)[0].lastStart());
  hlEC0_->Fill((*l1ts)[0].lastEventCounter0());
  hlHR_->Fill((*l1ts)[0].lastHardReset());  

  unsigned int lumisection = (*l1ts)[0].lumiSegmentNr();
  hphysTrig_->setBinContent(lumisection+1, (*l1ts)[0].l1AsPhysics());
  //Access whatever collections needed and fill histos...
  return ;
}

void
DQMScalInfo::makeDcsInfo(const edm::Event& e){

  edm::Handle<DcsStatusCollection> dcsStatus;
  e.getByToken(dcsStatusCollection_, dcsStatus);

  for(DcsStatusCollection::const_iterator dcsStatusItr = dcsStatus->begin();
      dcsStatusItr != dcsStatus->end(); ++dcsStatusItr){

    //    cout << __LINE__ << endl;
    if(!dcsStatusItr->ready(DcsStatus::CSCp))   dcs25[0]=false;
    if(!dcsStatusItr->ready(DcsStatus::CSCm))   dcs25[1]=false;
    if(!dcsStatusItr->ready(DcsStatus::DT0))    dcs25[2]=false;
    if(!dcsStatusItr->ready(DcsStatus::DTp))    dcs25[3]=false;
    if(!dcsStatusItr->ready(DcsStatus::DTm))    dcs25[4]=false;
    if(!dcsStatusItr->ready(DcsStatus::EBp))    dcs25[5]=false;
    if(!dcsStatusItr->ready(DcsStatus::EBm))    dcs25[6]=false;
    if(!dcsStatusItr->ready(DcsStatus::EEp))    dcs25[7]=false;
    if(!dcsStatusItr->ready(DcsStatus::EEm))    dcs25[8]=false;
    if(!dcsStatusItr->ready(DcsStatus::ESp))    dcs25[9]=false;
    if(!dcsStatusItr->ready(DcsStatus::ESm))    dcs25[10]=false;
    if(!dcsStatusItr->ready(DcsStatus::HBHEa))  dcs25[11]=false;
    if(!dcsStatusItr->ready(DcsStatus::HBHEb))  dcs25[12]=false;
    if(!dcsStatusItr->ready(DcsStatus::HBHEc))  dcs25[13]=false;
    if(!dcsStatusItr->ready(DcsStatus::HF))     dcs25[14]=false;
    if(!dcsStatusItr->ready(DcsStatus::HO))     dcs25[15]=false;
    if(!dcsStatusItr->ready(DcsStatus::BPIX))   dcs25[16]=false;
    if(!dcsStatusItr->ready(DcsStatus::FPIX))   dcs25[17]=false;
    if(!dcsStatusItr->ready(DcsStatus::RPC))    dcs25[18]=false;
    if(!dcsStatusItr->ready(DcsStatus::TIBTID)) dcs25[19]=false;
    if(!dcsStatusItr->ready(DcsStatus::TOB))    dcs25[20]=false;
    if(!dcsStatusItr->ready(DcsStatus::TECp))   dcs25[21]=false;
    if(!dcsStatusItr->ready(DcsStatus::TECm))   dcs25[22]=false;
    if(!dcsStatusItr->ready(DcsStatus::CASTOR)) dcs25[23]=false;
    if(!dcsStatusItr->ready(DcsStatus::ZDC))    dcs25[24]=false;
  }
  //  cout << __LINE__ << endl;
  return;
}

// void
// DQMScalInfo::makeGtInfo(const edm::Event& e){

//   edm::Handle<L1GlobalTriggerReadoutRecord> gtrr_handle;
//   if ( ! e.getByToken(gtCollection_, gtrr_handle) ){
//     dcs25[24]=false;
//     return;
//   }

//   if ( ! gtrr_handle.isValid() ){
//     edm::LogWarning("DQMDcsInfo") << " gtDigis not found" ;
//     dcs25[24]=false;
//     return;
//   }

//   L1GlobalTriggerReadoutRecord const* gtrr = gtrr_handle.product();
//   L1GtFdlWord fdlWord ;
//   if (gtrr)
//     fdlWord = gtrr->gtFdlWord();
//   else{
//     edm::LogWarning ("DQMDcsInfo") << " phys decl. bit not accessible !!!";
//     dcs25[24]=false;
//     return;
//   }

//   if (fdlWord.physicsDeclared() !=1) dcs25[24]=false;

//   return;
// }
