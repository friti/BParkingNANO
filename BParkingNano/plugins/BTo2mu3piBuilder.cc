#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TLorentzVector.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h"


#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

constexpr bool debugGen = false;

class BTo2mu3piBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  typedef std::vector<reco::GenParticle> GenParticleCollection;

  explicit BTo2mu3piBuilder(const edm::ParameterSet &cfg):
    k_selection_{cfg.getParameter<std::string>("kaonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    dz_value{cfg.getParameter<double>("dz_value")},
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    kaons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("kaons") )},
    kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    
    //GEN
    srcToken_(consumes<GenParticleCollection>(cfg.getParameter<edm::InputTag>("srcGen"))),
    //TRIGGER
    triggerBits_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(cfg.getParameter<edm::InputTag>("objects"))), 

    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>();
    }

  ~BTo2mu3piBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> k_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const double dz_value;

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> kaons_;
  const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;

  //GEN
  edm::EDGetTokenT<reco::GenParticleCollection> srcToken_;


  //TRIGGER                                                              
  const edm::EDGetTokenT<edm::TriggerResults> triggerBits_;                         
  const edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;

  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BTo2mu3piBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  evt.getByToken(dileptons_, dileptons);
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> kaons;
  evt.getByToken(kaons_, kaons);
  
  edm::Handle<TransientTrackCollection> kaons_ttracks;
  evt.getByToken(kaons_ttracks_, kaons_ttracks);  

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);

  //GEN
  int flag_jpsi_mu = -99. ;
  int flag_psi2s_mu = -99.;
  int flag_chic0_mu = -99.;
  int flag_chic1_mu = -99.;
  int flag_chic2_mu = -99.;
  int flag_hc_mu = -99.;
  int flag_jpsi_tau = -99.;
  int flag_psi2s_tau = -99.;
  int flag_jpsi_pi = -99.;
  int flag_jpsi_3pi = -99.;
  int flag_jpsi_hc = -99.;
  int flag_error = -99.;
  
  if(evt.eventAuxiliary().run() == 1){ //only if it is a MC sample
    //GEN
    if(debugGen) std::cout<<"In the flag code..."<<std::endl;
    edm::Handle<GenParticleCollection> src;
    evt.getByToken(srcToken_, src);
    const size_t n = src->size();
    std::vector<int> final_daus;

    if(debugGen) std::cout<<"Number of gen particles: "<<n<<std::endl;
    
    for(unsigned int  i = 0; i < n; ++i) {  //loop on gen particles
      const reco::GenParticle & gen = (*src)[i];
      const reco::Candidate* daughter; 
      const reco::Candidate* the_b;
      int is_doublemu = 0;
      int is_b = 0;
      final_daus.clear();
    
      if(abs(gen.pdgId()) == 443){  // looking for jpsi      

	if(debugGen) std::cout<<"There is a jpsi and she has "<<gen.numberOfDaughters()<<" daughters"<<std::endl;
	for(unsigned int dau = 0; dau < gen.numberOfDaughters(); dau++){  //loop of jpsi daughters
	  if(debugGen) std::cout<<"Jpsi daughter: "<<gen.daughter(dau)->pdgId()<<std::endl;
	  is_b = 0;
	  if (abs(gen.daughter(dau)->pdgId())==13){
	    is_doublemu += 1;
	  }
	} //end loop on daughters
	if(is_doublemu>=2){  // jpsi -> mu mu
	  if(debugGen) std::cout<<"The daughters are muons"<<std::endl;
	  the_b = gen.mother(0); // jpsi mother
	  if(abs(the_b->pdgId()) == 541){ //Bc->jpsi
	    if(debugGen) std::cout<<"The direct mother is a Bc"<<std::endl;
	    is_b = 1;
	  }  
	  else if(the_b->numberOfMothers() > 0){
	    the_b = gen.mother(0)->mother(0); // Bc->X->jpsi
	    if(abs(the_b->pdgId()) == 541 ){
	      if(debugGen) std::cout<<"The non direct mother is a Bc"<<std::endl;
	      is_b = 1;
	    }
	  }
	  if(is_b == 1){
	    
	    if(debugGen) std::cout<<"The Bc has "<<the_b->numberOfDaughters()<<"daughters"<<std::endl;

	    for(unsigned int bdau=0; bdau < the_b->numberOfDaughters(); bdau ++){
	      daughter = the_b->daughter(bdau);
	      if(abs(daughter->pdgId())!= 541 and abs(daughter->pdgId())!= 22){    //not gamma
		final_daus.push_back(abs(daughter->pdgId()));
		//	      cout<<daughter->pdgId()<<endl;
		if(debugGen) std::cout<<"The Bc daughters are "<< daughter->pdgId()<<std::endl;

	      }
	    }
      
	    std::sort(final_daus.begin(), final_daus.end());  //sort the pdgIds of the daughters
	    
	    flag_jpsi_mu = 0;
	    flag_psi2s_mu = 0;
	    flag_chic0_mu = 0;
	    flag_chic1_mu = 0;
	    flag_chic2_mu = 0;
	    flag_hc_mu = 0;
	    flag_jpsi_tau = 0;
	    flag_psi2s_tau = 0;
	    flag_jpsi_pi = 0;
	    flag_jpsi_3pi = 0;
	    flag_jpsi_hc = 0;
	    flag_error = 0;
	    
	    if(final_daus[0] == 13){  //muon
	      if(final_daus[1] == 14){
		if(final_daus[2] == 443)  flag_jpsi_mu=1;
		else if (final_daus[2] == 100443) flag_psi2s_mu = 1;
		else if (final_daus[2] == 10441) flag_chic0_mu = 1;
		else if (final_daus[2] == 20443) flag_chic1_mu = 1;
		else if (final_daus[2] == 445) flag_chic2_mu = 1;
		else if (final_daus[2] == 10443) flag_hc_mu = 1;
	      }
	      
	    }
	    else if(final_daus[0] == 15){ //tau
	      if (final_daus[1] == 16){
		if(final_daus[2] == 443) flag_jpsi_tau = 1;
		else if(final_daus[2] == 100443) flag_psi2s_tau = 1;
	      }
	    }
	    
	    else if(final_daus[0] == 211){
	      if (final_daus[1] == 443) flag_jpsi_pi = 1;
	      if (final_daus[1] == 211 && final_daus[2] ==211 && final_daus[1] == 443) flag_jpsi_3pi = 1;
	    }
	    
	    else if ((final_daus[0] == 431 && final_daus[1] == 443) || (final_daus[0] == 433 && final_daus[1] == 443)) flag_jpsi_hc = 1;  
	    else{
	      flag_error = 1;
	    }
	
	  } // if(is_b == 1)
	  
	} //if(is_doublemu>=2)
      }//if(abs(gen.pdgId()) == 443)
    }//for(unsigned int  i = 0; i < n; ++i)
  }//if(evt.eventAuxiliary().run() == 1)


  //TRIGGER                                                
  edm::Handle<edm::TriggerResults> triggerBits;                                     
  evt.getByToken(triggerBits_, triggerBits);               
  const edm::TriggerNames &names = evt.triggerNames(*triggerBits);  
  std::vector<pat::TriggerObjectStandAlone> triggeringMuons;        
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;   
  evt.getByToken(triggerObjects_, triggerObjects);                              

  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  std::vector<int> used_lep1_id, used_lep2_id, used_pi1_id,used_pi2_id,used_pi3_id;


  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  
  unsigned int index = names.triggerIndex("HLT_DoubleMu4_JpsiTrk_Displaced_v14");
  unsigned int index2 = names.triggerIndex("HLT_DoubleMu4_JpsiTrk_Displaced_v15");
  bool pass_hlt;
  if(index==triggerBits->size() && index2==triggerBits->size()){
    evt.put(std::move(ret_val));
  }                                                                               
  else{
    if(index!=triggerBits->size()){ 
      pass_hlt=triggerBits->accept(index);
    }
    else if (index2!=triggerBits->size()){
      pass_hlt=triggerBits->accept(index2);
    }
    std::vector<pat::TriggerObjectStandAlone> pass_jpsi;
    std::vector<pat::TriggerObjectStandAlone> pass_trk;
    if(pass_hlt){
      for (pat::TriggerObjectStandAlone obj : *triggerObjects){
	obj.unpackFilterLabels(evt, *triggerBits);
	obj.unpackPathNames(names);
	
	if(obj.hasFilterLabel("hltDisplacedmumuFilterDoubleMu4Jpsi")) {
	  pass_jpsi.push_back(obj);
	}
	else if(obj.hasFilterLabel("hltJpsiTkVertexFilter")) {
	  pass_trk.push_back(obj);
	}
	
      }
   

      //loop on jpsi muons (I need it here because I need dz condition on tracks)
      for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
	edm::Ptr<pat::CompositeCandidate> ll_prt(dileptons, ll_idx);
	edm::Ptr<reco::Candidate> l1_ptr = ll_prt->userCand("l1");
	edm::Ptr<reco::Candidate> l2_ptr = ll_prt->userCand("l2");
	int l1_idx = ll_prt->userInt("l1_idx");
	int l2_idx = ll_prt->userInt("l2_idx");
	
	int flag = 0;
	for(size_t pi1_idx = 0; pi1_idx < kaons->size(); ++pi1_idx) {
	  edm::Ptr<pat::CompositeCandidate> pi1_ptr(kaons, pi1_idx);
	  if( !k_selection_(*pi1_ptr) ) continue;
	  //dz requirement
	  if ( fabs(kaons_ttracks->at(pi1_idx).track().dz() - l1_ptr->bestTrack()->dz()) > 0.4 ||  fabs(kaons_ttracks->at(pi1_idx).track().dz() - l2_ptr->bestTrack()->dz()) > 0.4) continue;
	  //dR requirement
	  if(deltaR(leptons_ttracks->at(l2_idx).track().eta(),leptons_ttracks->at(l2_idx).track().phi(),kaons_ttracks->at(pi1_idx).track().eta(),kaons_ttracks->at(pi1_idx).track().phi()) < 0.1 ) continue;
	  if(deltaR(leptons_ttracks->at(l1_idx).track().eta(),leptons_ttracks->at(l1_idx).track().phi(),kaons_ttracks->at(pi1_idx).track().eta(),kaons_ttracks->at(pi1_idx).track().phi()) < 0.1 ) continue;

	  //matching online-offline della trk
	  for(pat::TriggerObjectStandAlone obj: pass_trk){
	    if(deltaR(obj,*pi1_ptr)<0.02) flag=1;
	  }    
	  
	  if(flag == 1){ 

	
	    math::PtEtaPhiMLorentzVector pi1_p4(
						pi1_ptr->pt(), 
						pi1_ptr->eta(),
						pi1_ptr->phi(),
						PI_MASS
						);
	    
	    for(size_t pi2_idx = 0; pi2_idx < kaons->size(); ++pi2_idx) {
	      edm::Ptr<pat::CompositeCandidate> pi2_ptr(kaons, pi2_idx);
	      if( !k_selection_(*pi2_ptr) ) continue;
	      if(pi2_idx == pi1_idx) continue;
	      // dz between track and leptons
	      if ( fabs(kaons_ttracks->at(pi2_idx).track().dz() - l1_ptr->bestTrack()->dz()) > 0.4 ||  fabs(kaons_ttracks->at(pi2_idx).track().dz() - l2_ptr->bestTrack()->dz()) > 0.4) continue;
	      //DR between tracks and leptons
	      if(deltaR(leptons_ttracks->at(l2_idx).track().eta(),leptons_ttracks->at(l2_idx).track().phi(),kaons_ttracks->at(pi2_idx).track().eta(),kaons_ttracks->at(pi2_idx).track().phi()) < 0.1 ) continue;
	      if(deltaR(leptons_ttracks->at(l1_idx).track().eta(),leptons_ttracks->at(l1_idx).track().phi(),kaons_ttracks->at(pi2_idx).track().eta(),kaons_ttracks->at(pi2_idx).track().phi()) < 0.1 ) continue;

	      math::PtEtaPhiMLorentzVector pi2_p4(
						  pi2_ptr->pt(), 
						  pi2_ptr->eta(),
						  pi2_ptr->phi(),
						  PI_MASS
						  );
	      for(size_t pi3_idx = 0; pi3_idx < kaons->size(); ++pi3_idx) {
		edm::Ptr<pat::CompositeCandidate> pi3_ptr(kaons, pi3_idx);
		if( !k_selection_(*pi3_ptr) ) continue;
		if(pi3_idx == pi1_idx or pi3_idx == pi2_idx) continue;
		//dz requirement
		if ( fabs(kaons_ttracks->at(pi3_idx).track().dz() - l1_ptr->bestTrack()->dz()) > 0.4 ||  fabs(kaons_ttracks->at(pi3_idx).track().dz() - l2_ptr->bestTrack()->dz()) > 0.4) continue;
		//DR requirements
		if(deltaR(leptons_ttracks->at(l2_idx).track().eta(),leptons_ttracks->at(l2_idx).track().phi(),kaons_ttracks->at(pi3_idx).track().eta(),kaons_ttracks->at(pi3_idx).track().phi()) < 0.1 ) continue;
		if(deltaR(leptons_ttracks->at(l1_idx).track().eta(),leptons_ttracks->at(l1_idx).track().phi(),kaons_ttracks->at(pi3_idx).track().eta(),kaons_ttracks->at(pi3_idx).track().phi()) < 0.1 ) continue;

		
		math::PtEtaPhiMLorentzVector pi3_p4(
						    pi3_ptr->pt(), 
						    pi3_ptr->eta(),
						    pi3_ptr->phi(),
						    PI_MASS
						    );
		
		//matching online- offline di jpsi
		flag = 0;
		for(pat::TriggerObjectStandAlone obj: pass_jpsi){
		  for(pat::TriggerObjectStandAlone obj2: pass_jpsi){  
		    if((deltaR(obj,*l1_ptr)<0.02) && (deltaR(obj2,*l2_ptr)<0.02))  flag=1;
		  }
		}
		if(flag==0) continue;
		else{
		  pat::CompositeCandidate cand;
		  cand.setP4(ll_prt->p4() + pi1_p4 + pi2_p4 + pi3_p4);
		  cand.setCharge(ll_prt->charge() + pi1_ptr->charge() + pi2_ptr->charge() + pi3_ptr->charge()  );
	      // Use UserCands as they should not use memory but keep the Ptr itself
	      // Put the lepton passing the corresponding selection
		  cand.addUserCand("l1", l1_ptr);
		  cand.addUserCand("l2", l2_ptr);
		  cand.addUserCand("pi1", pi1_ptr);
		  cand.addUserCand("pi2", pi2_ptr);
		  cand.addUserCand("pi3", pi3_ptr);
		  cand.addUserCand("dilepton", ll_prt);
	  
		  cand.addUserInt("l1_idx", l1_idx);
		  cand.addUserInt("l2_idx", l2_idx);
		  cand.addUserInt("pi1_idx", pi1_idx);
		  cand.addUserInt("pi2_idx", pi2_idx);
		  cand.addUserInt("pi3_idx", pi3_idx);
	      
		  auto dr_info = min_max_dr({l1_ptr, l2_ptr, pi1_ptr, pi2_ptr, pi3_ptr});
		  cand.addUserFloat("min_dr", dr_info.first);
		  cand.addUserFloat("max_dr", dr_info.second);

		  
		  //GEN variables
		  cand.addUserInt("is_jpsi_mu", flag_jpsi_mu);
		  cand.addUserInt("is_psi2s_mu", flag_psi2s_mu);
		  cand.addUserInt("is_chic0_mu", flag_chic0_mu);
		  cand.addUserInt("is_chic1_mu", flag_chic1_mu);
		  cand.addUserInt("is_chic2_mu", flag_chic2_mu);
		  cand.addUserInt("is_hc_mu", flag_hc_mu);
		  cand.addUserInt("is_jpsi_tau", flag_jpsi_tau);
		  cand.addUserInt("is_psi2s_tau", flag_psi2s_tau);
		  cand.addUserInt("is_jpsi_pi", flag_jpsi_pi);
		  cand.addUserInt("is_jpsi_3pi", flag_jpsi_3pi);
		  cand.addUserInt("is_jpsi_hc", flag_jpsi_hc);
		  cand.addUserInt("is_error", flag_error);
		  
		  int weight;
		  
		  if(evt.eventAuxiliary().run() == 1){
		    
		    if(flag_jpsi_mu == 1) weight = 1.;
		    else if(flag_psi2s_mu == 1) weight = 0.5474;
		    else if(flag_chic0_mu == 1) weight = 0.0116;
		    else if(flag_chic1_mu == 1) weight = 0.3440;
		    else if(flag_chic2_mu == 1) weight = 0.1950;
		    else if(flag_hc_mu == 1) weight = 0.01;
		    else if(flag_jpsi_tau == 1) weight = 1.;
		    else if(flag_psi2s_tau == 1) weight = 0.5474;
		    else if(flag_jpsi_pi == 1) weight = 1.;
		    else if(flag_jpsi_3pi == 1) weight = 1.;
		    else if(flag_jpsi_hc == 1) weight = 1.;
		    else weight = -1.;
		    
		  }
		  else weight = -99.;
		  cand.addUserFloat("weightGen", weight);
		  
		  
		  if( !pre_vtx_selection_(cand) ) continue;
		  
		  
		  KinVtxFitter fitter(
				      {leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx), kaons_ttracks->at(pi1_idx), kaons_ttracks->at(pi2_idx), kaons_ttracks->at(pi3_idx) },
				      {l1_ptr->mass(), l2_ptr->mass(), PI_MASS, PI_MASS, PI_MASS},
				      {LEP_SIGMA, LEP_SIGMA, PI_SIGMA, PI_SIGMA, PI_SIGMA} //some small sigma for the lepton mass
				      );
		  if(!fitter.success()) continue; // hardcoded, but do we need otherwise?

		  cand.setVertex( 
				 reco::Candidate::Point( 
							fitter.fitted_vtx().x(),
							fitter.fitted_vtx().y(),
							fitter.fitted_vtx().z()
							 )  
				  );
		  used_lep1_id.emplace_back(l1_idx);
		  used_lep2_id.emplace_back(l2_idx);
		  used_pi1_id.emplace_back(pi1_idx);
		  used_pi2_id.emplace_back(pi2_idx);
		  used_pi3_id.emplace_back(pi3_idx);
		  cand.addUserInt("sv_OK" , fitter.success());
		  cand.addUserFloat("sv_chi2", fitter.chi2());
		  cand.addUserFloat("sv_ndof", fitter.dof()); // float??
		  cand.addUserFloat("sv_prob", fitter.prob());
		  cand.addUserFloat("fitted_mll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass());
		  auto fit_p4 = fitter.fitted_p4();
		  cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
		  cand.addUserFloat("fitted_eta" , fit_p4.eta());
		  cand.addUserFloat("fitted_phi" , fit_p4.phi());
		  cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());      
		  cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
		  cand.addUserFloat(
				    "cos_theta_2D", 
				    cos_theta_2D(fitter, *beamspot, cand.p4())
				    );
		  cand.addUserFloat(
				    "fitted_cos_theta_2D", 
				    cos_theta_2D(fitter, *beamspot, fit_p4)
				    );
		  auto lxy = l_xy(fitter, *beamspot);
		  cand.addUserFloat("l_xy", lxy.value());
		  cand.addUserFloat("l_xy_unc", lxy.error());
		  cand.addUserFloat("vtx_x", cand.vx());
		  cand.addUserFloat("vtx_y", cand.vy());
		  cand.addUserFloat("vtx_z", cand.vz());
		  cand.addUserFloat("vtx_ex", sqrt(fitter.fitted_vtx_uncertainty().cxx()));
		  cand.addUserFloat("vtx_ey", sqrt(fitter.fitted_vtx_uncertainty().cyy()));
		  cand.addUserFloat("vtx_ez", sqrt(fitter.fitted_vtx_uncertainty().czz()));
		  
		  cand.addUserFloat("fitted_l1_pt" , fitter.daughter_p4(0).pt()); 
		  cand.addUserFloat("fitted_l1_eta", fitter.daughter_p4(0).eta());
		  cand.addUserFloat("fitted_l1_phi", fitter.daughter_p4(0).phi());
		  cand.addUserFloat("fitted_l2_pt" , fitter.daughter_p4(1).pt()); 
		  cand.addUserFloat("fitted_l2_eta", fitter.daughter_p4(1).eta());
		  cand.addUserFloat("fitted_l2_phi", fitter.daughter_p4(1).phi());
		  cand.addUserFloat("fitted_pi1_pt"  , fitter.daughter_p4(2).pt()); 
		  cand.addUserFloat("fitted_pi1_eta" , fitter.daughter_p4(2).eta());
		  cand.addUserFloat("fitted_pi1_phi" , fitter.daughter_p4(2).phi());
		  cand.addUserFloat("fitted_pi2_pt"  , fitter.daughter_p4(3).pt()); 
		  cand.addUserFloat("fitted_pi2_eta" , fitter.daughter_p4(3).eta());
		  cand.addUserFloat("fitted_pi2_phi" , fitter.daughter_p4(3).phi());
		  cand.addUserFloat("fitted_pi3_pt"  , fitter.daughter_p4(4).pt()); 
		  cand.addUserFloat("fitted_pi3_eta" , fitter.daughter_p4(4).eta());
		  cand.addUserFloat("fitted_pi3_phi" , fitter.daughter_p4(4).phi());

		  const reco::BeamSpot &bm = *beamspot;
		  
		  cand.addUserFloat("beamspot_x", bm.x0());
		  cand.addUserFloat("beamspot_y", bm.y0());
		  cand.addUserFloat("beamspot_z", bm.z0());
		  
		  
		  //new variables
		  
		  TLorentzVector P_b;
		  P_b.SetPtEtaPhiM(cand.pt(),cand.eta(),cand.phi(),cand.mass());
		  
		  TLorentzVector P_pi1;
		  P_pi1.SetPtEtaPhiM(pi1_ptr->pt(),pi1_ptr->eta(),pi1_ptr->phi(),pi1_ptr->mass());

		  TLorentzVector P_pi2;
		  P_pi2.SetPtEtaPhiM(pi2_ptr->pt(),pi2_ptr->eta(),pi2_ptr->phi(),pi2_ptr->mass());
		  
		  TLorentzVector P_pi3;
		  P_pi3.SetPtEtaPhiM(pi3_ptr->pt(),pi3_ptr->eta(),pi3_ptr->phi(),pi3_ptr->mass());

		  TLorentzVector P_l1;
		  P_l1.SetPtEtaPhiM(l1_ptr->pt(),l1_ptr->eta(),l1_ptr->phi(),l1_ptr->mass());
	      
		  TLorentzVector P_l2;
		  P_l2.SetPtEtaPhiM(l2_ptr->pt(),l2_ptr->eta(),l2_ptr->phi(),l2_ptr->mass());
		  
		  
		  float m_miss_2=(P_b-P_pi1 - P_pi2 - P_pi3 -P_l1-P_l2)*(P_b-P_pi1 - P_pi2 - P_pi3-P_l1-P_l2);
		  
		  float Q_2=(P_b-P_l1-P_l2)*(P_b-P_l1-P_l2);
	      
		  float pt_miss=(P_b.Pt()-P_pi1.Pt() - P_pi2.Pt() - P_pi3.Pt()-P_l1.Pt()-P_l2.Pt());
	      
		  float pt_miss_vec=((P_b-P_pi1 - P_pi2 - P_pi3-P_l1-P_l2).Pt());
		  
		  float pt_var=((P_l1+P_l2).Pt()-(P_pi1 - P_pi2 - P_pi3).Pt());
		  
		  float DR=deltaR(P_l1.Eta(),P_l1.Phi(),P_l2.Eta(),P_l2.Phi());
		  
		  //float deta = P_l1.Eta() - P_l2.Eta();
		  //float dphi = P_l1.Phi() - P_l2.Phi();
		  
		  //float DR_2=std::sqrt(deta * deta + dphi * dphi);	      
		  
		  
		  float m_jpsi=sqrt((P_l1+P_l2)*(P_l1+P_l2));
	      
		  cand.addUserFloat("m_miss_2", m_miss_2);
		  cand.addUserFloat("Q_2",Q_2);
		  
		  cand.addUserFloat("pt_miss",pt_miss);
		  cand.addUserFloat("pt_miss_vec",pt_miss_vec);
		  
		  
		  cand.addUserFloat("pt_var",pt_var);
		  
		  cand.addUserFloat("DR",DR);
		  
		  cand.addUserFloat("m_jpsi",m_jpsi);

		  
		  //energia del mu unpaired in diversi sistemi di riferimento                  
		  /*TLorentzVector P_mu=P_k;	      
		  
		  TVector3 mu_beta_lab=P_b.BoostVector();
		  
		  P_mu.Boost(-mu_beta_lab);
	      
		  cand.addUserFloat("E_mu_star",P_mu.E());
		  
		  P_mu=P_k;	      
		  
		  TLorentzVector jpsi=P_l1+P_l2;
		  
		  TVector3 jpsi_beta_lab=jpsi.BoostVector();
		  
		  P_mu.Boost(-jpsi_beta_lab);
		  
		  cand.addUserFloat("E_mu_#",P_mu.E());
		  */

		  if( !post_vtx_selection_(cand) ) continue;        

		  //compute isolation
		  float l1_iso03 = 0;
		  float l1_iso04 = 0;
		  float l2_iso03 = 0;
		  float l2_iso04 = 0;
		  float pi1_iso03  = 0;
		  float pi1_iso04  = 0;
		  float pi2_iso03  = 0;
		  float pi2_iso04  = 0;
		  float pi3_iso03  = 0;
		  float pi3_iso04  = 0;
		  float b_iso03  = 0;
		  float b_iso04  = 0;
		  
		  for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
		    
		    const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
		    // define selections for iso tracks (pT, eta, ...)
		    if( !isotrk_selection_(trk) ) continue;

		    // only consider tracks originating close to the three bodies
		    if ( !l1_ptr->bestTrack() || fabs(trk.dz() - l1_ptr->bestTrack()->dz()) > 0.4 ) continue;
		    if ( !l2_ptr->bestTrack() || fabs(trk.dz() - l2_ptr->bestTrack()->dz()) > 0.4 ) continue;
		    if ( fabs(trk.dz() - kaons_ttracks->at(pi1_idx).track().dz()) > 0.4 ) continue;
		    if ( fabs(trk.dz() - kaons_ttracks->at(pi2_idx).track().dz()) > 0.4 ) continue;
		    if (  fabs(trk.dz() - kaons_ttracks->at(pi3_idx).track().dz()) > 0.4 ) continue;
		    
		    if(track_to_lepton_match(pi1_ptr, iso_tracks.id(), iTrk)  )  continue;
		    if(track_to_lepton_match(pi2_ptr, iso_tracks.id(), iTrk)  )  continue;
		    if(track_to_lepton_match(pi3_ptr, iso_tracks.id(), iTrk)  )  continue;

		    // check if the track is one of the two leptons
		    if (track_to_lepton_match(l1_ptr, iso_tracks.id(), iTrk) || 
			track_to_lepton_match(l2_ptr, iso_tracks.id(), iTrk) ) continue;
		
		    // add to final particle iso if dR < cone
		    float dr_to_l1 = deltaR(cand.userFloat("fitted_l1_eta"), cand.userFloat("fitted_l1_phi"), trk.eta(), trk.phi());
		    float dr_to_l2 = deltaR(cand.userFloat("fitted_l2_eta"), cand.userFloat("fitted_l2_phi"), trk.eta(), trk.phi());
		    float dr_to_pi1  = deltaR(cand.userFloat("fitted_pi1_eta") , cand.userFloat("fitted_pi1_phi") , trk.eta(), trk.phi());
		    float dr_to_pi2  = deltaR(cand.userFloat("fitted_pi2_eta") , cand.userFloat("fitted_pi2_phi") , trk.eta(), trk.phi());
		    float dr_to_pi3  = deltaR(cand.userFloat("fitted_pi3_eta") , cand.userFloat("fitted_pi3_phi") , trk.eta(), trk.phi());
		    float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());
		
		    if (dr_to_l1 < 0.4 && dr_to_l1>0.01){
		      l1_iso04 += trk.pt();
		      if ( dr_to_l1 < 0.3) l1_iso03 += trk.pt();
		    }
		    if (dr_to_l2 < 0.4 && dr_to_l2>0.01){
		      l2_iso04 += trk.pt();
		      if (dr_to_l2 < 0.3)  l2_iso03 += trk.pt();
		    }
		    if (dr_to_pi1 < 0.4 && dr_to_pi1>0.01){
		      pi1_iso04 += trk.pt();
		      if (dr_to_pi1 < 0.3) pi1_iso03 += trk.pt();
		    }
		    if (dr_to_pi2 < 0.4 && dr_to_pi2>0.01){
		      pi2_iso04 += trk.pt();
		      if (dr_to_pi2 < 0.3) pi2_iso03 += trk.pt();
		    }
		    if (dr_to_pi3 < 0.4 && dr_to_pi3>0.01){
		      pi3_iso04 += trk.pt();
		      if (dr_to_pi3 < 0.3) pi3_iso03 += trk.pt();
		    }
		    if (dr_to_b < 0.4){
		      b_iso04 += trk.pt();
		      if (dr_to_b < 0.3) b_iso03 += trk.pt();
		    }
		    }
		  cand.addUserFloat("l1_iso03", l1_iso03);
		  cand.addUserFloat("l1_iso04", l1_iso04);
		  cand.addUserFloat("l2_iso03", l2_iso03);
		  cand.addUserFloat("l2_iso04", l2_iso04);
		  cand.addUserFloat("pi1_iso03" , pi1_iso03 );
		  cand.addUserFloat("pi1_iso04" , pi1_iso04 );
		  cand.addUserFloat("pi2_iso03" , pi2_iso03 );
		  cand.addUserFloat("pi2_iso04" , pi2_iso04 );
		  cand.addUserFloat("pi3_iso03" , pi3_iso03 );
		  cand.addUserFloat("pi3_iso04" , pi3_iso04 );
		  cand.addUserFloat("b_iso03" , b_iso03 );
		  cand.addUserFloat("b_iso04" , b_iso04 );

		  ret_val->push_back(cand);
		}//trigger dileptons
	      } // loop pion3

	    } //loop p2
	  } // track trigger
	} // loop p1

      }// loop muons from jpsi

      for (auto & cand: *ret_val){
	cand.addUserInt("n_pi1_used", std::count(used_pi1_id.begin(),used_pi1_id.end(),cand.userInt("pi1_idx")));
	cand.addUserInt("n_pi2_used", std::count(used_pi2_id.begin(),used_pi2_id.end(),cand.userInt("pi2_idx")));
	cand.addUserInt("n_pi3_used", std::count(used_pi3_id.begin(),used_pi3_id.end(),cand.userInt("pi3_idx")));
	cand.addUserInt("n_l1_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l1_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l1_idx")));
	cand.addUserInt("n_l2_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l2_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l2_idx")));
	} 
	}//pass_hlt
    evt.put(std::move(ret_val));
  }//index ok
}//produce

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(BTo2mu3piBuilder);

