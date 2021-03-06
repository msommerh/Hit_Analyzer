
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

// For jet
#include "DataFormats/JetReco/interface/PFJet.h" 
#include "DataFormats/JetReco/interface/PFJetCollection.h" 

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

class HitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit HitAnalyzer(const edm::ParameterSet& conf);
  ~HitAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

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

  void reset( void );
  edm::ParameterSet conf_;
  edm::InputTag src_;
  bool printLocal;
  bool phase1_;
  bool isMC_;
  double pT_cut_;
       
  edm::EDGetTokenT< reco::GenParticleCollection>          genPtoken;
  edm::EDGetTokenT< reco::PFJetCollection >               ak4CHStoken;
  edm::EDGetTokenT< edmNew::DetSetVector<SiPixelCluster>> clusterToken;
  edm::EDGetTokenT< reco::VertexCollection >              pvToken;
 

  // ----------member data ---------------------------
  
  int                  nJets;
  std::vector<double>  jet_pt;
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
  std::vector<int>     genParticle_status;
  
  std::vector< int>           detUnit_layer ; // 1-4
  std::vector< double >           detUnit_X     ;
  std::vector< double >           detUnit_Y     ;
  std::vector< double >           detUnit_Z     ;
       
  std::vector<int>                  nClusters;
  std::vector<std::vector<double>>  cluster_globalz;
  std::vector<std::vector<double>>  cluster_globalx;
  std::vector<std::vector<double>>  cluster_globaly;
  std::vector<std::vector<double>>  cluster_charge;     
       
  TTree *tree;
       
};


HitAnalyzer::HitAnalyzer(const edm::ParameterSet& conf)
: conf_(conf), src_(conf.getParameter<edm::InputTag>( "src" )) { 
 
  printLocal = conf.getUntrackedParameter<bool>("Verbosity",false);
  phase1_ = conf.getUntrackedParameter<bool>("phase1",false);
  isMC_ = conf.getUntrackedParameter<bool>("isMC",false);
  pT_cut_ = conf.getUntrackedParameter<double>("pT_cut",0);
     
  // initialization
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>( "tree", "tree" );

  tree->Branch( "nJets"             , &nJets );
  tree->Branch( "jet_pt"            , &jet_pt );

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

  if (isMC_) {
    tree->Branch( "jet_MC_bTag"       , &jet_MC_bTag );
  }
  
  std::string labelgenP("genParticles");
  std::string labelAK4s("ak4PFJetsCHS");
  std::string labelClusters("siPixelClusters");
  std::string labelPVs("offlinePrimaryVertices");

  genPtoken      = consumes<reco::GenParticleCollection         > (edm::InputTag(labelgenP));
  ak4CHStoken    = consumes<reco::PFJetCollection               > (edm::InputTag(labelAK4s));
  clusterToken   = consumes<edmNew::DetSetVector<SiPixelCluster>> (src_);
  pvToken        = consumes<reco::VertexCollection              > (edm::InputTag(labelPVs));
}


HitAnalyzer::~HitAnalyzer()
{
 
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
  Handle<reco::PFJetCollection                > ak4CHS   ; iEvent.getByToken( ak4CHStoken  , ak4CHS   );
  Handle<reco::GenParticleCollection          > genPs    ; iEvent.getByToken( genPtoken    , genPs    );
  Handle<reco::VertexCollection		      > PVs	 ; iEvent.getByToken( pvToken	   , PVs      );

  if (ak4CHS->begin() == ak4CHS->end()) return; //filter out events without any jets

  //Loop over PVs
  for (reco::VertexCollection::const_iterator pv=PVs->begin(); pv != PVs->end(); ++pv){
    PV_x.push_back(pv->x());
    PV_y.push_back(pv->y());
    PV_z.push_back(pv->z());
    ++nPV;
  }


  if (isMC_){ 
    // Loop opver gen particles for later truth matching
    TLorentzVector selectedGenP ;
    for ( reco::GenParticleCollection::const_iterator gp = genPs->begin(); gp != genPs->end(); ++gp ) {
        if (gp->status() != 23 && gp->status() !=2) continue;
	if (gp->status() == 23 && fabs(gp->pdgId()) == 5) continue;
	if (gp->status() == 2 && !(fabs(gp->pdgId()) > 500 && fabs(gp->pdgId())<600) && !(fabs(gp->pdgId()) > 5000 && fabs(gp->pdgId())<6000)) continue;
        TLorentzVector genP;
        genP.SetPtEtaPhiM(gp->pt(),gp->eta(),gp->phi(),gp->mass());
        int  pdgId = gp->pdgId();
        int  status = gp->status();
        genParticle_pt  .push_back(genP.Pt());
        genParticle_eta .push_back(genP.Eta());
        genParticle_phi .push_back(genP.Phi());
        genParticle_mass.push_back(genP.M());
        genParticle_pdgId.push_back(pdgId);
        genParticle_status.push_back(status);
    }
    nGenParticles = genParticle_pt.size();
  }
  
  // Get vector of detunit ids and loop
  const edmNew::DetSetVector<SiPixelCluster>& input = *clusters;
  
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
    int layer  = 0; // 1-4

    if(subid==2) {  // forward
#ifdef NEW_ID
      PixelEndcapName pen(detid,tTopo,phase1_);
#else 
      PXFDetId pdetId = PXFDetId(detid); 
#endif
    }
    else if (subid==1) {  // barrel
#ifdef NEW_ID
      PixelBarrelName pbn(detid,tTopo,phase1_);
#else      
      PXBDetId pdetId = PXBDetId(detid);
      PixelBarrelName pbn(pdetId); // Convert to online     
#endif
      layer  = pbn.layerName();
    }
    
    std::vector<double>  _cluster_globalz;
    std::vector<double>  _cluster_globalx;
    std::vector<double>  _cluster_globaly;
    std::vector<double>  _cluster_charge;
    int numberOfClusters = 0;
    for ( edmNew::DetSet<SiPixelCluster>::const_iterator clustIt = detUnit->begin(); clustIt != detUnit->end(); ++clustIt ) {
      numberOfClusters++;

      // get global position of the cluster
      float x = clustIt->x(); // row, cluster position in pitch units, as float (int+0.5);
      float y = clustIt->y(); // column, analog average
      LocalPoint lp = topol->localPosition(MeasurementPoint(x,y));

      double charge = clustIt->charge();
	
      GlobalPoint clustgp = theGeomDet->surface().toGlobal( lp );
      double gZ = clustgp.z();  // global z
      double gX = clustgp.x();
      double gY = clustgp.y();
      
      _cluster_globalz.push_back(gZ);
      _cluster_globalx.push_back(gX);
      _cluster_globaly.push_back(gY);
      _cluster_charge.push_back(charge);
    }
    if( numberOfClusters < 1.) continue;
    
    detUnit_layer .push_back( layer ); 
    detUnit_X         .push_back(detX     );
    detUnit_Y         .push_back(detY     );
    detUnit_Z         .push_back(detZ     );
    
    nClusters.push_back(numberOfClusters);
    cluster_globalz  .push_back(_cluster_globalz);
    cluster_globalx  .push_back(_cluster_globalx);
    cluster_globaly  .push_back(_cluster_globaly);
    cluster_charge   .push_back(_cluster_charge);
  }

  // Loop over jets
  for ( reco::PFJetCollection::const_iterator jet = ak4CHS->begin(); jet != ak4CHS->end(); ++jet ) {
    if (jet->pt()<pT_cut_) continue;
    ++nJets;
    jet_pt.push_back(jet->pt());
    TLorentzVector TVjet;
    TVjet.SetPtEtaPhiM(jet->pt(),jet->eta(),jet->phi(),jet->mass());

    const std::vector<double> PV{PV_x[0], PV_y[0], PV_z[0]};
    ClusterMatcher(TVjet, nClusters, cluster_globalx, cluster_globaly, cluster_globalz, detUnit_layer, PV, nClusters_L1004, nClusters_L1006, nClusters_L1008, nClusters_L1010, nClusters_L1016, nClusters_L2004, nClusters_L2006, nClusters_L2008, nClusters_L2010, nClusters_L2016, nClusters_L3004, nClusters_L3006, nClusters_L3008, nClusters_L3010, nClusters_L3016, nClusters_L4004, nClusters_L4006, nClusters_L4008, nClusters_L4010, nClusters_L4016);
    
    //true MC btag (truth matching)
    if (isMC_){
      int true_bTag = -10;
      for (int p=0; p != nGenParticles ; ++p){
	if (genParticle_pt[p] < 350) continue;
	if (genParticle_status[p] == 2){
          TLorentzVector pvector;
          pvector.SetPtEtaPhiM(genParticle_pt[p],genParticle_eta[p],genParticle_phi[p],genParticle_mass[p]);
          if (TVjet.DeltaR(pvector) < 0.3){
            true_bTag = 1;
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
    }

  }

  tree->Fill();

}

// Private methods
void HitAnalyzer::reset( void ){

  nJets = 0;
  jet_pt.clear();
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
  genParticle_status.clear();
  
  nClusters.clear();
  cluster_globalz.clear();
  cluster_globalx.clear();
  cluster_globaly.clear();
  cluster_charge.clear();	
  
  detUnit_layer .clear();
  
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

// ------------ method matches pixel clusters to a jet  ------------
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
    double cluster_phi = GetPhi(cluster_x[i][j] - pv[0], cluster_y[i][j] - pv[1]); // subtract PV coordinates to use them as the origin of the spherical coordinate system
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

// ------------ helper method that converts cartesian X,Y to sperical phi  ------------
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

// ------------ helper method that converts cartesian X,Y to spherical theta  ------------
double HitAnalyzer::GetTheta(const double X, const double Y, const double Z){
  return 0.5*TMath::Pi()-TMath::ATan(Z/TMath::Sqrt(X*X + Y*Y));
}

// ------------ helper method that returns deltaR from spherical coordinates  ------------
double HitAnalyzer::dR_theta_phi(const double &theta1, const double &theta2, const double &phi1, const double &phi2){
  return TMath::Sqrt(TMath::Sq(TMath::Log(TMath::Tan(0.5*theta1))-TMath::Log(TMath::Tan(0.5*theta2))) + TMath::Sq(phi1-phi2));
}

// ------------ helper method that fills pixel clusters into the right dR slots  ------------
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

//define this as a plug-in
DEFINE_FWK_MODULE(HitAnalyzer);
