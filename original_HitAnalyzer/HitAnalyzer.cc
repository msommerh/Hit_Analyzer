
// system include files
#include <memory>
#include <string> 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// For event setup
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

// ROOT includes
#include "TLorentzVector.h" 
#include "TMath.h" 
#include "TTree.h" 

// For output
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// For btag
#include "DataFormats/BTauReco/interface/JetTag.h"


// For jet
#include "DataFormats/JetReco/interface/PFJet.h" 
#include "DataFormats/JetReco/interface/PFJetCollection.h" 

// #include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
// #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
// #include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"

// For pixel clusters
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"

// Pixel topology
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#define NEW_ID
#ifdef NEW_ID
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#else 
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#endif 


#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"


//  For gen particle
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// For vertices
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// For triggers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//for MET
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

class HitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit HitAnalyzer(const edm::ParameterSet& conf);
  ~HitAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  //const std::vector<reco::Candidate> findDaughters(const reco::GenParticle *particle); 

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void ClusterMatcher(const TLorentzVector &jvector, 
  const std::vector<int> &n_clusters,
  const std::vector<std::vector<double>>  &cluster_x,
  const std::vector<std::vector<double>>  &cluster_y,
  const std::vector<std::vector<double>>  &cluster_z,
  const std::vector< int> &Unit_layer,
  const std::vector<double> &pv,
  std::vector<int> &_nClusters_L1004, 
  std::vector<int> &_nClusters_L1006,
  std::vector<int> &_nClusters_L1008,
  std::vector<int> &_nClusters_L1010,
  std::vector<int> &_nClusters_L1016,
  std::vector<int> &_nClusters_L2004,
  std::vector<int> &_nClusters_L2006,
  std::vector<int> &_nClusters_L2008,
  std::vector<int> &_nClusters_L2010,
  std::vector<int> &_nClusters_L2016,
  std::vector<int> &_nClusters_L3004,
  std::vector<int> &_nClusters_L3006,
  std::vector<int> &_nClusters_L3008,
  std::vector<int> &_nClusters_L3010,
  std::vector<int> &_nClusters_L3016,
  std::vector<int> &_nClusters_L4004,
  std::vector<int> &_nClusters_L4006,
  std::vector<int> &_nClusters_L4008,
  std::vector<int> &_nClusters_L4010,
  std::vector<int> &_nClusters_L4016

  );
  double GetPhi(const double X, const double Y);
  double GetTheta(const double X, const double Y, const double Z);
  double dR_theta_phi(const double &theta1, const double &theta2, const double &phi1, const double &phi2);
  void AddMatchedClusters(double DR, int &nC004, int &nC006, int &nC008, int &nC010, int &nC016);
  bool LooseJetsID(reco::PFJetCollection::const_iterator _jet);

  void reset( void );
  const reco::GenParticle* findMother(const reco::GenParticle *particle);
  //for Zprime fraction   
  const reco::GenParticle* findMother(std::vector<reco::GenParticle>::const_iterator& particle);   
  //end Zprime fraction
  edm::ParameterSet conf_;
  edm::InputTag src_;
  edm::InputTag HLTtriggers_;
  bool printLocal;
  bool phase1_;
  bool isMC_;
  double pT_cut_;
  int nJets_cut_;
  double leading_jet_eta_;
  bool loose_jets_cut_;
  bool tight_jets_cut_;
  double MET_over_sumEt_cut_;

       
  edm::EDGetTokenT< reco::GenParticleCollection>          genPtoken;
  edm::EDGetTokenT< reco::PFJetCollection >               ak4CHStoken;
  edm::EDGetTokenT< reco::PFJetCollection >               ak8CHStoken;
  edm::EDGetTokenT< edmNew::DetSetVector<SiPixelCluster>> clusterToken;
  edm::EDGetTokenT< reco::VertexCollection >              svToken;
  edm::EDGetTokenT< reco::VertexCollection >              pvToken;
  edm::EDGetTokenT< reco::TrackCollection >               trackToken;
  edm::EDGetTokenT< reco::JetTagCollection >              csv2Token;
  edm::EDGetTokenT< edm::TriggerResults >		  HLTtriggersToken;
  edm::EDGetTokenT< reco::PFMETCollection >		  METToken;

  //trying to implement deepCSV
  edm::EDGetTokenT< reco::JetTagCollection >		deepCSVToken_probb;
  edm::EDGetTokenT< reco::JetTagCollection >		deepCSVToken_probbb;
  //end deepCSV
 

  // ----------member data ---------------------------
  
  double	       MET_over_sumEt;
 
  int                  nJets;
  double	       dijet_mass;
  std::vector<double>  jet_pt;
  std::vector<double>  jet_eta;
  std::vector<double>  jet_phi;
  std::vector<double>  jet_mass;
  std::vector<int>     jet_pdgId;
  std::vector<double>  jet_bTag;

  //trying to implement deepCSV:
  std::vector<double>  jet_deepCSV_probb;
  std::vector<double>  jet_deepCSV_probbb;
  //end deepCSV

  std::vector<int>     jet_MC_bTag;

  std::vector<int>     nClusters_L1004;
  std::vector<int>     nClusters_L1006;
  std::vector<int>     nClusters_L1008;
  std::vector<int>     nClusters_L1010;
  std::vector<int>     nClusters_L1016;
  std::vector<int>     nClusters_L2004;
  std::vector<int>     nClusters_L2006;
  std::vector<int>     nClusters_L2008;
  std::vector<int>     nClusters_L2010;
  std::vector<int>     nClusters_L2016;
  std::vector<int>     nClusters_L3004;
  std::vector<int>     nClusters_L3006;
  std::vector<int>     nClusters_L3008;
  std::vector<int>     nClusters_L3010;
  std::vector<int>     nClusters_L3016;
  std::vector<int>     nClusters_L4004;
  std::vector<int>     nClusters_L4006;
  std::vector<int>     nClusters_L4008;
  std::vector<int>     nClusters_L4010;
  std::vector<int>     nClusters_L4016;

  int 		       nPV;
  std::vector<double>  PV_x;
  std::vector<double>  PV_y;
  std::vector<double>  PV_z; 
  
  int                  nGenParticles;
  std::vector<double>  genParticle_pt;
  std::vector<double>  genParticle_eta;
  std::vector<double>  genParticle_phi;
  std::vector<double>  genParticle_mass;
  std::vector<int>     genParticle_pdgId;
  std::vector<int>     genParticle_mother_pdgId;
  std::vector<int>     genParticle_status;
  std::vector<double>  genParticle_vx_x   ;
  std::vector<double>  genParticle_vx_y   ;
  std::vector<double>  genParticle_vx_z   ;
  std::vector<double>  genParticle_vx_eta ;
  std::vector<double>  genParticle_vx_phi ;
  std::vector<double>  genParticle_vx_r   ;
  std::vector<double>  genParticle_decayvx_x   ;
  std::vector<double>  genParticle_decayvx_y   ;
  std::vector<double>  genParticle_decayvx_z   ;
  std::vector<double>  genParticle_decayvx_eta ;
  std::vector<double>  genParticle_decayvx_phi ;
  std::vector<double>  genParticle_decayvx_r   ;
  //for Zprime fraction
  std::vector<int>	genParticle_motherID;
  std::vector<int>      jet_BMotherID;
  //end Zprime fraction
  int                          nDetUnits;
  std::vector<unsigned int>    detUnit_subdetId;
  std::vector< int>           detUnit_layer ; // 1-3
  std::vector<unsigned int>   detUnit_disk;   //1,2,3
  std::vector<unsigned int>   detUnit_side;   //size=1 for -z, 2 for +z
  std::vector< double >           detUnit_X     ;
  std::vector< double >           detUnit_Y     ;
  std::vector< double >           detUnit_Z     ;
  std::vector< double >           detUnit_R     ;
  std::vector< double >           detUnit_Phi   ;
       
  std::vector<int>                  nClusters;
  std::vector<std::vector<double>>  cluster_sizex;
  std::vector<std::vector<double>>  cluster_sizey;
  std::vector<std::vector<double>>  cluster_globalz;
  std::vector<std::vector<double>>  cluster_globalx;
  std::vector<std::vector<double>>  cluster_globaly;
  std::vector<std::vector<double>>  cluster_globalPhi;
  std::vector<std::vector<double>>  cluster_globalR;
  std::vector<std::vector<double>>  cluster_charge;     
       
  TTree *tree;
       
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HitAnalyzer::HitAnalyzer(const edm::ParameterSet& conf)
: conf_(conf), src_(conf.getParameter<edm::InputTag>( "src" )), HLTtriggers_(conf.getParameter<edm::InputTag>( "HLTtriggers" )) { 
//: conf_(conf), src_(conf.getParameter<edm::InputTag>( "src" )) { 
 
  printLocal = conf.getUntrackedParameter<bool>("Verbosity",false);
  phase1_ = conf.getUntrackedParameter<bool>("phase1",false);
  isMC_ = conf.getUntrackedParameter<bool>("isMC",false);
  pT_cut_ = conf.getUntrackedParameter<double>("pT_cut",0);
  nJets_cut_ = conf.getUntrackedParameter<int>("nJets_cut",0);
  leading_jet_eta_ = conf.getUntrackedParameter<double>("leading_jet_eta", 2.5);
  loose_jets_cut_ = conf.getUntrackedParameter<bool>("loose_jets_cut", false);
  tight_jets_cut_ = conf.getUntrackedParameter<bool>("tight_jets_cut", false);
  MET_over_sumEt_cut_ = conf.getUntrackedParameter<double>("MET_over_sumEt_cut", 1.);
     
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>( "tree", "tree" );

  tree->Branch( "MET_over_sumEt"    , &MET_over_sumEt);

  tree->Branch( "nJets"             , &nJets );
  tree->Branch( "dijet_mass"        , &dijet_mass );
  tree->Branch( "jet_pt"            , &jet_pt );
  tree->Branch( "jet_eta"           , &jet_eta );
  tree->Branch( "jet_phi"           , &jet_phi );
  tree->Branch( "jet_mass"          , &jet_mass );
  tree->Branch( "jet_bTag"          , &jet_bTag );

  //trying to implement deepCSV:
  tree->Branch( "jet_deepCSV_probb"  , &jet_deepCSV_probb );
  tree->Branch( "jet_deepCSV_probbb"  , &jet_deepCSV_probbb );
  //end deepCSV

  tree->Branch( "nClusters_L1004", &nClusters_L1004 ); 
  tree->Branch( "nClusters_L1006", &nClusters_L1006 );
  tree->Branch( "nClusters_L1008", &nClusters_L1008 );
  tree->Branch( "nClusters_L1010", &nClusters_L1010 );
  tree->Branch( "nClusters_L1016", &nClusters_L1016 );
  tree->Branch( "nClusters_L2004", &nClusters_L2004 );
  tree->Branch( "nClusters_L2006", &nClusters_L2006 );
  tree->Branch( "nClusters_L2008", &nClusters_L2008 );
  tree->Branch( "nClusters_L2010", &nClusters_L2010 );
  tree->Branch( "nClusters_L2016", &nClusters_L2016 );
  tree->Branch( "nClusters_L3004", &nClusters_L3004 );
  tree->Branch( "nClusters_L3006", &nClusters_L3006 );
  tree->Branch( "nClusters_L3008", &nClusters_L3008 );
  tree->Branch( "nClusters_L3010", &nClusters_L3010 );
  tree->Branch( "nClusters_L3016", &nClusters_L3016 );
  tree->Branch( "nClusters_L4004", &nClusters_L4004 );
  tree->Branch( "nClusters_L4006", &nClusters_L4006 );
  tree->Branch( "nClusters_L4008", &nClusters_L4008 );
  tree->Branch( "nClusters_L4010", &nClusters_L4010 );
  tree->Branch( "nClusters_L4016", &nClusters_L4016 );
  
  tree->Branch("nPV",  &nPV);
  tree->Branch("PV_x", &PV_x);
  tree->Branch("PV_y", &PV_y);
  tree->Branch("PV_z", &PV_z); 


  if (isMC_) {
    tree->Branch( "jet_MC_bTag"       , &jet_MC_bTag );
    tree->Branch( "nGenParticles"     , &nGenParticles );
    tree->Branch( "genParticle_pt"    , &genParticle_pt );
    tree->Branch( "genParticle_eta"   , &genParticle_eta );
    tree->Branch( "genParticle_phi"   , &genParticle_phi );
    tree->Branch( "genParticle_mass"  , &genParticle_mass );
    tree->Branch( "genParticle_pdgId" , &genParticle_pdgId );
    tree->Branch( "genParticle_mother_pdgId" , &genParticle_mother_pdgId );
    tree->Branch( "genParticle_status", &genParticle_status );
    tree->Branch( "genParticle_vx_x"   , &genParticle_vx_x   );
    tree->Branch( "genParticle_vx_y"   , &genParticle_vx_y   );
    tree->Branch( "genParticle_vx_z"   , &genParticle_vx_z   );
    tree->Branch( "genParticle_vx_eta" , &genParticle_vx_eta );
    tree->Branch( "genParticle_vx_phi" , &genParticle_vx_phi );
    tree->Branch( "genParticle_vx_r"   , &genParticle_vx_r   );
//    tree->Branch( "genParticle_decayvx_x"   , &genParticle_decayvx_x   );
//    tree->Branch( "genParticle_decayvx_y"   , &genParticle_decayvx_y   );
//    tree->Branch( "genParticle_decayvx_z"   , &genParticle_decayvx_z   );
//    tree->Branch( "genParticle_decayvx_eta" , &genParticle_decayvx_eta );
//    tree->Branch( "genParticle_decayvx_phi" , &genParticle_decayvx_phi );
//    tree->Branch( "genParticle_decayvx_r"   , &genParticle_decayvx_r   );
    //for Zprime fraction
    tree->Branch( "genParticle_motherID" , &genParticle_motherID );
    tree->Branch( "jet_BMotherID"	 , &jet_BMotherID	);
    //end Zprime fraction
    
  }
  
  std::string labelgenP("genParticles");
  std::string labelAK8s("ak8PFJetsCHS");
  std::string labelAK4s("ak4PFJetsCHS");
  std::string labelClusters("siPixelClusters");
  std::string labelSVs("inclusiveSecondaryVertices");
  std::string labelPVs("offlinePrimaryVertices");
  std::string labelTracks("generalTracks");
  std::string labelCSV("pfCombinedSecondaryVertexV2BJetTags");
  std::string labelHLTtriggers("TriggerResults");	
  std::string labelMET("pfMet");

  genPtoken      = consumes<reco::GenParticleCollection         > (edm::InputTag(labelgenP));
  ak8CHStoken    = consumes<reco::PFJetCollection               > (edm::InputTag(labelAK8s));
  ak4CHStoken    = consumes<reco::PFJetCollection               > (edm::InputTag(labelAK4s));
  clusterToken   = consumes<edmNew::DetSetVector<SiPixelCluster>> (src_);
  svToken        = consumes<reco::VertexCollection              > (edm::InputTag(labelSVs));
  pvToken        = consumes<reco::VertexCollection              > (edm::InputTag(labelPVs));
  csv2Token      = consumes<reco::JetTagCollection              > (edm::InputTag(labelCSV));
  trackToken     = consumes<reco::TrackCollection               > (edm::InputTag(labelTracks));
  HLTtriggersToken=consumes<edm::TriggerResults	    		> (HLTtriggers_);
  METToken	 = consumes<reco::PFMETCollection		> (edm::InputTag(labelMET));

  //trying to implement filters
  //std::string HCALlaserNoiseFilter_Selector_ = "Flag_hcalLaserEventFilter";
  //end filters

  //trying to implement deepCSV
  //std::string labelDeepCSV("pfDeepCSVJetTags");//, "probb", "RECO");
  deepCSVToken_probb  = consumes<reco::JetTagCollection		> (edm::InputTag("pfDeepCSVJetTags", "probb", "RECO"));
  deepCSVToken_probbb  = consumes<reco::JetTagCollection	> (edm::InputTag("pfDeepCSVJetTags", "probbb", "RECO"));
  //end deepCSV
}


HitAnalyzer::~HitAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
  HitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  using namespace edm;
   
  // Make sure vectors are clear
  reset();
   
  // // Get event setup
  edm::ESHandle<TrackerGeometry> geom;
  iSetup.get<TrackerDigiGeometryRecord>().get( geom );
  const TrackerGeometry& theTracker(*geom);
  
#ifdef NEW_ID
  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopoH;
  iSetup.get<TrackerTopologyRcd>().get(tTopoH);
  const TrackerTopology *tTopo=tTopoH.product();
#endif
   
  // Get handles
  Handle<edmNew::DetSetVector<SiPixelCluster> > clusters ; iEvent.getByToken( clusterToken , clusters );
  Handle<reco::PFJetCollection                > ak8CHS   ; iEvent.getByToken( ak8CHStoken  , ak8CHS   );
  Handle<reco::PFJetCollection                > ak4CHS   ; iEvent.getByToken( ak4CHStoken  , ak4CHS   );
  Handle<reco::GenParticleCollection          > genPs    ; iEvent.getByToken( genPtoken    , genPs    );
  Handle<reco::JetTagCollection               > CSVs     ; iEvent.getByToken( csv2Token    , CSVs     );
  Handle<reco::VertexCollection		      > PVs	 ; iEvent.getByToken( pvToken	   , PVs      );
  Handle<reco::VertexCollection               > SVs      ; iEvent.getByToken( svToken      , SVs      );
  Handle<reco::TrackCollection                > tracks   ; iEvent.getByToken( trackToken   , tracks   );
  Handle<edm::TriggerResults		      > HLTtriggers; iEvent.getByToken(HLTtriggersToken, HLTtriggers);
  Handle<reco::PFMETCollection		      > METs	 ; iEvent.getByToken( METToken	   , METs     );
  const reco::JetTagCollection & bTags = *(CSVs.product()); 

  //trying to implement filters:
  //event.getByToken(noiseFilterToken_, noiseFilterBits_);
  //const edm::TriggerNames &names = event.triggerNames(*noiseFilterBits_);
  //end filters

  //trying to implement deepCSV:
  Handle<reco::JetTagCollection               > deepCSVs_probb     ; iEvent.getByToken( deepCSVToken_probb, deepCSVs_probb );
  const reco::JetTagCollection & deepbTags_probb = *(deepCSVs_probb.product());
  Handle<reco::JetTagCollection               > deepCSVs_probbb     ; iEvent.getByToken( deepCSVToken_probbb, deepCSVs_probbb );
  const reco::JetTagCollection & deepbTags_probbb = *(deepCSVs_probbb.product());
  //end deepCSV

  // to print all available HLTs:
  //const edm::TriggerNames& trigNames = iEvent.triggerNames(*HLTtriggers);
  //for (unsigned int i=0; i<HLTtriggers->size(); ++i) {
  // std::cout<<"trigger number "<<i<<": "<<trigNames.triggerName(i)<<std::endl;
  //}

  MET_over_sumEt = (METs->front()).et()/(METs->front()).sumEt();

  if (ak4CHS->begin() == ak4CHS->end()) return; //filter out events without any jets

  //Loop over PVs
  for (reco::VertexCollection::const_iterator pv=PVs->begin(); pv != PVs->end(); ++pv){
    //if (pv->ndof() <= 4 || fabs(pv->z()) >= 24 || fabs(pv->position().rho()) > 2) continue;
    PV_x.push_back(pv->x());
    PV_y.push_back(pv->y());
    PV_z.push_back(pv->z());
    ++nPV;
  }

  //Jet selection loop:

  std::vector<reco::PFJetCollection::const_iterator> selected_jets;
  
  reco::PFJetCollection::const_iterator leading_jet1;
  reco::PFJetCollection::const_iterator leading_jet2;
  if (ak4CHS->begin()->pt() > (ak4CHS->begin()+1)->pt()){ //initialize search for the two leading jets by comparing the first two
    leading_jet1 = ak4CHS->begin();
    leading_jet2 = ak4CHS->begin()+1;
  }
  else {
    leading_jet2 = ak4CHS->begin();
    leading_jet1 = ak4CHS->begin()+1;
  }

  for ( reco::PFJetCollection::const_iterator jet = ak4CHS->begin(); jet != ak4CHS->end(); ++jet ) { 
    if (jet->pt()<pT_cut_ || (loose_jets_cut_ && !LooseJetsID(jet))) continue;
    if (jet != leading_jet1 && jet != leading_jet2){  //compare two the leading jets to see if the new jet has higher pt
      if (jet->pt() >= leading_jet1->pt()){
        leading_jet2 = leading_jet1;
        leading_jet1 = jet;
      }
      else if (jet->pt() > leading_jet2->pt()){
        leading_jet2 = jet;
      }
    }

    selected_jets.push_back(jet);
    ++nJets;
  }

  //event selection
  if (MET_over_sumEt>=MET_over_sumEt_cut_) return;
  if (nJets < nJets_cut_ || nPV < 1 || fabs(leading_jet1->eta()) >= leading_jet_eta_ || fabs(leading_jet2->eta()) >= leading_jet_eta_) return;
  if (tight_jets_cut_){
    if (fabs(leading_jet1->eta()) < 2.4 && (leading_jet1->neutralHadronEnergyFraction() >= 0.9 || leading_jet1->neutralEmEnergyFraction() >= 0.9)) return;
    if (fabs(leading_jet2->eta()) < 2.4 && (leading_jet2->neutralHadronEnergyFraction() >= 0.9 || leading_jet2->neutralEmEnergyFraction() >= 0.9)) return;
  }

   if (isMC_){ 
    // Loop opver gen particles
    TLorentzVector selectedGenP ;
    for ( reco::GenParticleCollection::const_iterator gp = genPs->begin(); gp != genPs->end(); ++gp ) {
        if (gp->status() != 23 && gp->status() !=2) continue;
        //if (fabs(gp->pdgId()) == 5 || (!(fabs(gp->pdgId()) > 500 && fabs(gp->pdgId())<600) && !(fabs(gp->pdgId()) > 5000 && fabs(gp->pdgId())<6000))) continue;
        //if (gp->pt() < 350) continue;	
	if (gp->status() == 23 && fabs(gp->pdgId()) == 5) continue;
	if (gp->status() == 2 && !(fabs(gp->pdgId()) > 500 && fabs(gp->pdgId())<600) && !(fabs(gp->pdgId()) > 5000 && fabs(gp->pdgId())<6000)) continue;
        TLorentzVector genP;
        genP.SetPtEtaPhiM(gp->pt(),gp->eta(),gp->phi(),gp->mass());
        int  pdgId = gp->pdgId();
        int  status = gp->status();
        double vx_x   = gp->vertex().Coordinates().X(); 
        double vx_y   = gp->vertex().Coordinates().Y(); 
        double vx_z   = gp->vertex().Coordinates().Z(); 
        double vx_eta = gp->vertex().Coordinates().Eta(); 
        double vx_phi = gp->vertex().Coordinates().Phi(); 
        double vx_r   = gp->vertex().Coordinates().R(); 
        genParticle_pt  .push_back(genP.Pt());
        genParticle_eta .push_back(genP.Eta());
        genParticle_phi .push_back(genP.Phi());
        genParticle_mass.push_back(genP.M());
        genParticle_pdgId.push_back(pdgId);
        genParticle_status.push_back(status);
        genParticle_vx_x   .push_back(vx_x   );
        genParticle_vx_y   .push_back(vx_y   );
        genParticle_vx_z   .push_back(vx_z   );
        genParticle_vx_eta .push_back(vx_eta );
        genParticle_vx_phi .push_back(vx_phi );
        genParticle_vx_r   .push_back(vx_r   );
	//for Zprime fraction
	if (gp->status() == 2) genParticle_motherID.push_back(findMother(findMother(gp))->pdgId());
	//end Zprime fraction
    }
    nGenParticles = genParticle_pt.size();
  }
  
  // Get vector of detunit ids and loop
  const edmNew::DetSetVector<SiPixelCluster>& input = *clusters;
  
  int numberOfDetUnits = 0;
  
  for ( edmNew::DetSetVector<SiPixelCluster>::const_iterator detUnit = input.begin(); detUnit != input.end(); ++detUnit ) {
    if (detUnit->size()<1) continue;
    unsigned int detid = detUnit->detId();
    DetId detId = DetId(detid);       // Get the Detid object
    unsigned int detType=detId.det(); // det type, pixel=1
    if(detType!=1) continue; // look only at pixels
    unsigned int subid=detId.subdetId(); //subdetector type, pix barrel=1, forward=2
    // Subdet id, pix barrel=1, forward=2
    
    // Get the geom-detector
    const PixelGeomDetUnit * theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detId) );
    const PixelTopology * topol = &(theGeomDet->specificTopology());
    double detX = theGeomDet->surface().position().x();
    double detY = theGeomDet->surface().position().y();
    double detZ = theGeomDet->surface().position().z();
    double detR = theGeomDet->surface().position().perp();
    double detPhi = theGeomDet->surface().position().phi();
    // barrel ids
    // unsigned int layerC=0;
    // unsigned int ladderC=0;
    // unsigned int zindex=0;
    // int shell  = 0; // shell id // Shell { mO = 1, mI = 2 , pO =3 , pI =4 };
    // int sector = 0; // 1-8
    // int ladder = 0; // 1-22
    int layer  = 0; // 1-3
    // int module = 0; // 1-4
    // bool half  = false; //

    // Endcap ids
    unsigned int disk=0; //1,2,3
    // unsigned int blade=0; //1-24
    // unsigned int moduleF=0; // plaquette 1,2,3,4
    unsigned int side=0; //size=1 for -z, 2 for +z
    // unsigned int panel=0; //panel=1
   
   
    if(subid==2) {  // forward
#ifdef NEW_ID
      disk=tTopo->pxfDisk(detid); //1,2,3
      // blade=tTopo->pxfBlade(detid); //1-24
      // zindex=tTopo->pxfModule(detid); //
      side=tTopo->pxfSide(detid); //size=1 for -z, 2 for +z
      // panel=tTopo->pxfPanel(detid); //panel=1
      PixelEndcapName pen(detid,tTopo,phase1_);
#else 
      PXFDetId pdetId = PXFDetId(detid); 
      disk=pdetId.disk();      //1,2,3
      // blade=pdetId.blade();    //1-24
      // moduleF=pdetId.module(); // plaquette
      side=pdetId.side();      //size=1 for -z, 2 for +z
      // panel=pdetId.panel();    //panel=1
#endif
    }
    else if (subid==1) {  // barrel
#ifdef NEW_ID
      // layerC=tTopo->pxbLayer(detid);
      // ladderC=tTopo->pxbLadder(detid);
      // zindex=tTopo->pxbModule(detid);
      PixelBarrelName pbn(detid,tTopo,phase1_);
#else      
      PXBDetId pdetId = PXBDetId(detid);
     
      // layerC=pdetId.layer();       // Barell layer = 1,2,3
      // ladderC=pdetId.ladder();     // Barrel ladder id 1-20,32,44.
      // zindex=pdetId.module();      // Barrel Z-index=1,8
      PixelBarrelName pbn(pdetId); // Convert to online     
#endif

      // PixelBarrelName::Shell sh = pbn.shell(); //enum
      // sector = pbn.sectorName();
      // ladder = pbn.ladderName();
      layer  = pbn.layerName();
      // module = pbn.moduleName();
      // half  = pbn.isHalfModule();
      // shell = int(sh);
    
      // if(shell==1 || shell==2) module = -module; // change the module sign for z<0
      // if(shell==1 || shell==3) ladder = -ladder; // change ladeer sign for Outer )x<0)
    }
    
    
    std::vector<double>  _cluster_sizex;
    std::vector<double>  _cluster_sizey;
    std::vector<double>  _cluster_globalz;
    std::vector<double>  _cluster_globalx;
    std::vector<double>  _cluster_globaly;
    std::vector<double>  _cluster_globalPhi;
    std::vector<double>  _cluster_globalR;
    std::vector<double>  _cluster_charge;
    int numberOfClusters = 0;
    for ( edmNew::DetSet<SiPixelCluster>::const_iterator clustIt = detUnit->begin(); clustIt != detUnit->end(); ++clustIt ) {
      numberOfClusters++;

      // get global position of the cluster
      int sizeX = clustIt->sizeX(); //x=row=rfi, 
      int sizeY = clustIt->sizeY(); //y=col=z_global
      float x = clustIt->x(); // row, cluster position in pitch units, as float (int+0.5);
      float y = clustIt->y(); // column, analog average
      LocalPoint lp = topol->localPosition(MeasurementPoint(x,y));
      // float lx = lp.x(); // local cluster position in cm
      // float ly = lp.y();

      double charge = clustIt->charge();
	
      GlobalPoint clustgp = theGeomDet->surface().toGlobal( lp );
      double gZ = clustgp.z();  // global z
      double gX = clustgp.x();
      double gY = clustgp.y();
      TVector3 v(gX,gY,gZ);
      float gPhi = v.Phi(); // phi of the hit
      float gR = v.Perp(); // r of the hit
      
      _cluster_sizex.push_back(sizeX);
      _cluster_sizey.push_back(sizeY);
      _cluster_globalz.push_back(gZ);
      _cluster_globalx.push_back(gX);
      _cluster_globaly.push_back(gY);
      _cluster_globalPhi.push_back(gPhi);
      _cluster_globalR.push_back(gR);
      _cluster_charge.push_back(charge);
    }
    if( numberOfClusters < 1.) continue;
    
    detUnit_subdetId.push_back(subid);  
    detUnit_layer .push_back( layer ); 
    detUnit_disk.push_back(disk);   
    detUnit_side.push_back(side);   
    detUnit_X         .push_back(detX     );
    detUnit_Y         .push_back(detY     );
    detUnit_Z         .push_back(detZ     );
    detUnit_R         .push_back(detR     );
    detUnit_Phi       .push_back(detPhi   );
    
    nClusters.push_back(numberOfClusters);
    cluster_sizex    .push_back(_cluster_sizex);
    cluster_sizey    .push_back(_cluster_sizey);
    cluster_globalz  .push_back(_cluster_globalz);
    cluster_globalx  .push_back(_cluster_globalx);
    cluster_globaly  .push_back(_cluster_globaly);
    cluster_globalPhi.push_back(_cluster_globalPhi);
    cluster_globalR  .push_back(_cluster_globalR);
    cluster_charge   .push_back(_cluster_charge);

    numberOfDetUnits++;
  }
  nDetUnits = numberOfDetUnits;

 
  // Loop over jets
  for ( reco::PFJetCollection::const_iterator jet : selected_jets ) {
  //for ( reco::PFJetCollection::const_iterator jet = ak4CHS->begin(); jet != ak4CHS->end(); ++jet ) {
    //if (jet->pt()<pT_cut_) continue;
    TLorentzVector TVjet;
    TVjet.SetPtEtaPhiM(jet->pt(),jet->eta(),jet->phi(),jet->mass());
    jet_pt.push_back(jet->pt());
    jet_eta.push_back(jet->eta());
    jet_phi.push_back(jet->phi());
    jet_mass.push_back(jet->mass());

    const std::vector<double> PV{PV_x[0], PV_y[0], PV_z[0]};
    ClusterMatcher(TVjet, nClusters, cluster_globalx, cluster_globaly, cluster_globalz, detUnit_layer, PV, nClusters_L1004, nClusters_L1006, nClusters_L1008, nClusters_L1010, nClusters_L1016, nClusters_L2004, nClusters_L2006, nClusters_L2008, nClusters_L2010, nClusters_L2016, nClusters_L3004, nClusters_L3006, nClusters_L3008, nClusters_L3010, nClusters_L3016, nClusters_L4004, nClusters_L4006, nClusters_L4008, nClusters_L4010, nClusters_L4016);
    
    //true MC btag
    if (isMC_){
      int true_bTag = -10;
      //Zprime fraction
      int bJet_mother = 0;
      //end Zprime fraction
      for (int p=0; p != nGenParticles ; ++p){
        //if (!(fabs(genParticle_pdgId[p])>500 && fabs(genParticle_pdgId[p])<600) && !(fabs(genParticle_pdgId[p])>5000 && fabs(genParticle_pdgId[p])<6000)) continue;
        //if (genParticle_status[p] != 2) continue;
	if (genParticle_pt[p] < 350) continue;
	if (genParticle_status[p] == 2){
          TLorentzVector pvector;
          pvector.SetPtEtaPhiM(genParticle_pt[p],genParticle_eta[p],genParticle_phi[p],genParticle_mass[p]);
          if (TVjet.DeltaR(pvector) < 0.3){
            true_bTag = 1;
	    //checking Zprime fraction
	    bJet_mother = genParticle_motherID[p];
	    //end Zprime fraction
            break;
          }
	}
	else if (genParticle_status[p] == 23){
          TLorentzVector pvector;
          pvector.SetPtEtaPhiM(genParticle_pt[p],genParticle_eta[p],genParticle_phi[p],genParticle_mass[p]);
          if (TVjet.DeltaR(pvector) < 0.3){
            true_bTag = 0;
            break;
            }	
	}
      }
    jet_MC_bTag.push_back(true_bTag);
    //Zprime fraction
    jet_BMotherID.push_back(bJet_mother);
    //end Zprime fraction
    }

    // CSV b-tag infos
    double match = 0.4;
    double csv2 = -99.;
    for (unsigned int i = 0; i != bTags.size(); ++i) {
      //if (bTags[i].first->pt()<170.) continue; //this sets all values in current sample to -99
      TLorentzVector bTagJet;
      bTagJet.SetPtEtaPhiM(bTags[i].first->pt(),bTags[i].first->eta(),bTags[i].first->phi(),bTags[i].first->mass());
      float dR = TVjet.DeltaR(bTagJet);
      if (dR > match ) continue;
      match = dR;
      csv2 = bTags[i].second; 
    }
    jet_bTag.push_back(csv2);

    // trying to implement deepCSV
    double deepmatch = 0.4;
    double deepcsv_probb = -99.;
    double deepcsv_probbb = -99.;
    for (unsigned int i = 0; i != deepbTags_probb.size(); ++i) {
      TLorentzVector deepbTagJet;
      deepbTagJet.SetPtEtaPhiM(deepbTags_probb[i].first->pt(),deepbTags_probb[i].first->eta(),deepbTags_probb[i].first->phi(),deepbTags_probb[i].first->mass());
      float dR = TVjet.DeltaR(deepbTagJet);
      if (dR > deepmatch ) continue;
      deepmatch = dR;
      deepcsv_probb = deepbTags_probb[i].second;
      deepcsv_probbb = deepbTags_probbb[i].second;
    }
    jet_deepCSV_probb.push_back(deepcsv_probb);
    jet_deepCSV_probbb.push_back(deepcsv_probbb);
    //end deepCSV

  }

  //dijet mass
  TLorentzVector TVLeadingJet1, TVLeadingJet2;
  TVLeadingJet1.SetPtEtaPhiM(leading_jet1->pt(),leading_jet1->eta(),leading_jet1->phi(),leading_jet1->mass());
  TVLeadingJet2.SetPtEtaPhiM(leading_jet2->pt(),leading_jet2->eta(),leading_jet2->phi(),leading_jet2->mass());
  dijet_mass = (TVLeadingJet1 + TVLeadingJet2).M();

  tree->Fill();
}

// Private methods
void HitAnalyzer::reset( void ){

  MET_over_sumEt = 0;

  nJets = 0;
  dijet_mass = 0.;
  jet_pt.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_mass.clear();
  jet_pdgId.clear();
  jet_bTag.clear();

  //trying to implement deepCSV:
  jet_deepCSV_probb.clear();
  jet_deepCSV_probbb.clear();
  //end deepCSV
  
  jet_MC_bTag.clear();

  nClusters_L1004.clear(); 
  nClusters_L1006.clear();
  nClusters_L1008.clear();
  nClusters_L1010.clear();
  nClusters_L1016.clear();
  nClusters_L2004.clear();
  nClusters_L2006.clear();
  nClusters_L2008.clear();
  nClusters_L2010.clear();
  nClusters_L2016.clear();
  nClusters_L3004.clear();
  nClusters_L3006.clear();
  nClusters_L3008.clear();
  nClusters_L3010.clear();
  nClusters_L3016.clear();
  nClusters_L4004.clear();
  nClusters_L4006.clear();
  nClusters_L4008.clear();
  nClusters_L4010.clear();
  nClusters_L4016.clear();

  nPV = 0;
  PV_x.clear();
  PV_y.clear();
  PV_z.clear();
  
  nGenParticles = 0;
  genParticle_pt.clear();
  genParticle_eta.clear();
  genParticle_phi.clear();
  genParticle_mass.clear();
  genParticle_pdgId.clear();
  genParticle_mother_pdgId.clear();
  genParticle_status.clear();
  genParticle_vx_x  .clear(); 
  genParticle_vx_y  .clear(); 
  genParticle_vx_z  .clear(); 
  genParticle_vx_eta.clear(); 
  genParticle_vx_phi.clear(); 
  genParticle_vx_r  .clear(); 
  genParticle_decayvx_x  .clear(); 
  genParticle_decayvx_y  .clear(); 
  genParticle_decayvx_z  .clear(); 
  genParticle_decayvx_eta.clear(); 
  genParticle_decayvx_phi.clear(); 
  genParticle_decayvx_r  .clear();
  //for Zprime fraction
  genParticle_motherID.clear();
  jet_BMotherID.clear();
  //end Zprime fraction 
  
  nClusters.clear();
  cluster_sizex.clear();
  cluster_sizey.clear();
  cluster_globalz.clear();
  cluster_globalx.clear();
  cluster_globaly.clear();
  cluster_globalPhi.clear();
  cluster_globalR.clear();
  cluster_charge.clear();
  
  nDetUnits = 0.;
  detUnit_subdetId.clear();
  detUnit_layer .clear();
  detUnit_disk.clear();
  detUnit_side.clear();
  detUnit_layer .clear();
  detUnit_X     .clear();
  detUnit_Y     .clear();
  detUnit_Z     .clear();
  detUnit_R     .clear();
  detUnit_Phi   .clear();
  
}

// Find mother
//for Zprime fraction
const reco::GenParticle *HitAnalyzer::findMother(std::vector<reco::GenParticle>::const_iterator& particle) {
    const reco::GenParticle *tmp = &(*particle);
    int pdgId = particle->pdgId();
    while(const reco::GenParticle *mother = dynamic_cast<const reco::GenParticle *>(tmp->mother())) {
      if(mother->pdgId() == pdgId) {
        tmp = mother;
        continue;
      }
        return mother;
      }
      return 0;
    }
//end Zprime fraction
const reco::GenParticle *HitAnalyzer::findMother(const reco::GenParticle *particle) {
      const reco::GenParticle *tmp = particle;
      int pdgId = particle->pdgId();
      while(const reco::GenParticle *mother = dynamic_cast<const reco::GenParticle *>(tmp->mother())) {
        if(mother->pdgId() == pdgId) {
          tmp = mother;
          continue;
        }
        return mother;
      }
      return 0;
    }
              
// ------------ method called once each job just before starting event loop  ------------
void 
  HitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
  HitAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ method matches clusters to a jet  ------------
void 
HitAnalyzer::ClusterMatcher(const TLorentzVector &jvector, 
  const std::vector<int> &n_clusters,
  const std::vector<std::vector<double>>  &cluster_x,
  const std::vector<std::vector<double>>  &cluster_y,
  const std::vector<std::vector<double>>  &cluster_z,
  const std::vector<int> &Unit_layer, 
  const std::vector<double> &pv,
  std::vector<int> &_nClusters_L1004, 
  std::vector<int> &_nClusters_L1006,
  std::vector<int> &_nClusters_L1008,
  std::vector<int> &_nClusters_L1010,
  std::vector<int> &_nClusters_L1016,
  std::vector<int> &_nClusters_L2004,
  std::vector<int> &_nClusters_L2006,
  std::vector<int> &_nClusters_L2008,
  std::vector<int> &_nClusters_L2010,
  std::vector<int> &_nClusters_L2016,
  std::vector<int> &_nClusters_L3004,
  std::vector<int> &_nClusters_L3006,
  std::vector<int> &_nClusters_L3008,
  std::vector<int> &_nClusters_L3010,
  std::vector<int> &_nClusters_L3016,
  std::vector<int> &_nClusters_L4004,
  std::vector<int> &_nClusters_L4006,
  std::vector<int> &_nClusters_L4008,
  std::vector<int> &_nClusters_L4010,
  std::vector<int> &_nClusters_L4016
 
 ){
const double jtheta = jvector.Theta();
const double jphi = jvector.Phi();
int nCl_L1004 = 0; 
int nCl_L1006 = 0;
int nCl_L1008 = 0;
int nCl_L1010 = 0;
int nCl_L1016 = 0;
int nCl_L2004 = 0;
int nCl_L2006 = 0;
int nCl_L2008 = 0;
int nCl_L2010 = 0;
int nCl_L2016 = 0;
int nCl_L3004 = 0;
int nCl_L3006 = 0;
int nCl_L3008 = 0;
int nCl_L3010 = 0;
int nCl_L3016 = 0;
int nCl_L4004 = 0;
int nCl_L4006 = 0;
int nCl_L4008 = 0;
int nCl_L4010 = 0;
int nCl_L4016 = 0;

for ( unsigned int i = 0; i != n_clusters.size(); ++i ){ //loop over detUnits
  for (int j = 0; j != n_clusters[i]; ++j){  //loop over clusters in each derUnit
    double cluster_phi = GetPhi(cluster_x[i][j] - pv[0], cluster_y[i][j] - pv[1]);
    double cluster_theta = GetTheta(cluster_x[i][j] - pv[0], cluster_y[i][j] - pv[1], cluster_z[i][j] - pv[2]);
    double dR = dR_theta_phi(jtheta, cluster_theta, jphi, cluster_phi);
    switch (Unit_layer[i]) {
      case 1:
        AddMatchedClusters(dR, nCl_L1004, nCl_L1006, nCl_L1008, nCl_L1010, nCl_L1016);
        break;
      case 2:
        AddMatchedClusters(dR, nCl_L2004, nCl_L2006, nCl_L2008, nCl_L2010, nCl_L2016);
        break;
      case 3:
        AddMatchedClusters(dR, nCl_L3004, nCl_L3006, nCl_L3008, nCl_L3010, nCl_L3016);
        break;
      case 4:
        AddMatchedClusters(dR, nCl_L4004, nCl_L4006, nCl_L4008, nCl_L4010, nCl_L4016);
        break;
    }
  }
}

_nClusters_L1004.push_back(nCl_L1004); 
_nClusters_L1006.push_back(nCl_L1006);
_nClusters_L1008.push_back(nCl_L1008);
_nClusters_L1010.push_back(nCl_L1010);
_nClusters_L1016.push_back(nCl_L1016);
_nClusters_L2004.push_back(nCl_L2004);
_nClusters_L2006.push_back(nCl_L2006);
_nClusters_L2008.push_back(nCl_L2008);
_nClusters_L2010.push_back(nCl_L2010);
_nClusters_L2016.push_back(nCl_L2016);
_nClusters_L3004.push_back(nCl_L3004);
_nClusters_L3006.push_back(nCl_L3006);
_nClusters_L3008.push_back(nCl_L3008);
_nClusters_L3010.push_back(nCl_L3010);
_nClusters_L3016.push_back(nCl_L3016);
_nClusters_L4004.push_back(nCl_L4004);
_nClusters_L4006.push_back(nCl_L4006);
_nClusters_L4008.push_back(nCl_L4008);
_nClusters_L4010.push_back(nCl_L4010);
_nClusters_L4016.push_back(nCl_L4016);
}

double HitAnalyzer::GetPhi(const double X, const double Y){
    if (X>0) return TMath::ATan(Y/X);
    else if (X<0){
      if (Y >= 0) return TMath::ATan(Y/X) + TMath::Pi();
      else return TMath::ATan(Y/X) - TMath::Pi();
      }
    else{
      if (Y>=0) return 0.5*TMath::Pi();
      else return -0.5*TMath::Pi();
      }
}

double HitAnalyzer::GetTheta(const double X, const double Y, const double Z){
  return 0.5*TMath::Pi()-TMath::ATan(Z/TMath::Sqrt(X*X + Y*Y));
}

double HitAnalyzer::dR_theta_phi(const double &theta1, const double &theta2, const double &phi1, const double &phi2){
  return TMath::Sqrt(TMath::Sq(TMath::Log(TMath::Tan(0.5*theta1))-TMath::Log(TMath::Tan(0.5*theta2))) + TMath::Sq(phi1-phi2));
}

void HitAnalyzer::AddMatchedClusters(double DR, int &nC004, int &nC006, int &nC008, int &nC010, int &nC016){
  if (DR<0.16){
    ++nC016;
    if (DR<0.10){
      ++nC010;
      if (DR<0.08){
        ++nC008;
        if (DR<0.06){
          ++nC006;
          if (DR<0.04){
            ++nC004;
          }
        }
      }
    }
  }

}

bool HitAnalyzer::LooseJetsID(reco::PFJetCollection::const_iterator _jet){
    if (_jet->neutralHadronEnergyFraction() < 0.99 && _jet->neutralEmEnergyFraction() < 0.99 && _jet->chargedMultiplicity()+_jet->neutralMultiplicity() > 1){
      if (fabs(_jet->eta())<2.4 && (_jet->chargedHadronEnergyFraction() == 0 || _jet->chargedEmEnergyFraction() >= 0.99 || _jet->chargedMultiplicity() == 0)){
        return false;  
      }
      else {
        return true;
      }
    }
    else {
      return false;
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HitAnalyzer);
