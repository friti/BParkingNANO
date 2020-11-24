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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "PhysicsTools/HepMCCandAlgos/interface/PdgEntryReplacer.h"
#include "DataFormats/Common/interface/Association.h"


#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

constexpr bool debugGen = false;

class BTommmBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;
  typedef std::vector<reco::GenParticle> GenParticleCollection;

  explicit BTommmBuilder(const edm::ParameterSet &cfg):
    k_selection_{cfg.getParameter<std::string>("kaonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    kaons_{consumes<pat::MuonCollection>( cfg.getParameter<edm::InputTag>("kaons") )},
    kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),


    //GEN
    srcToken_(consumes<GenParticleCollection>(cfg.getParameter<edm::InputTag>("srcGen"))), 

    //TRIGGER                       
    triggerBits_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>(\
    "bits"))),                                                                      
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(cfg.getParameter<edm::InputTag>("objects"))),                                                  



    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>();
    }

  ~BTommmBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::Muon> k_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<pat::MuonCollection> kaons_;
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

void BTommmBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  evt.getByToken(dileptons_, dileptons);
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<pat::MuonCollection> kaons;
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
	    /*
	    for(unsigned int item=0; item< final_daus.size(); item ++){
	      if(debugGen) std::cout<<final_daus[item]<<std::endl;
	      if(item == final_daus.size() -1) if(debugGen) std::cout<<" "<<std::endl;
	    }
	    */

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

  std::vector<int> used_lep1_id, used_lep2_id, used_trk_id;


  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  //choose the HLT path
  unsigned int index = names.triggerIndex("HLT_Dimuon0_Jpsi3p5_Muon2_v5"); 
  
  if(index==triggerBits->size()){
    //    std::cout<<"Non ha HLT path giusto"<<std::endl;
    evt.put(std::move(ret_val));
    
  }                                                                               
  else if(index!=triggerBits->size()){
    
    //salvo gli oggetti di trigger che mi interessano
    std::vector<pat::TriggerObjectStandAlone> pass_jpsi;
    std::vector<pat::TriggerObjectStandAlone> pass_3mu;
    //std::cout<<"ha HLT path giusto"<<std::endl;
    bool pass_hlt=triggerBits->accept(index);
    /*if(!pass_hlt) {
      std::cout<<"non in pass_hlt"<<std::endl;
      std::cout<<"ret_val size"<<ret_val->size()<<std::endl;
      evt.put(std::move(ret_val));
      }*/
    if(pass_hlt){
      //std::cout<<"entra in pass_hlt"<<std::endl;
      for (pat::TriggerObjectStandAlone obj : *triggerObjects){
	obj.unpackFilterLabels(evt, *triggerBits);
	obj.unpackPathNames(names);
	
	if(obj.hasFilterLabel("hltVertexmumuFilterJpsiMuon3p5")){
	  pass_jpsi.push_back(obj);
	  //std::cout<<"filtro jpsi="<<obj.pt()<<std::endl;
	}
	
	if (obj.hasFilterLabel("hltTripleMuL3PreFiltered222")){
	  pass_3mu.push_back(obj);
	  //std::cout<<"filtro 3mu="<<obj.pt()<<std::endl;
	}	
      }
      
    


      //riconoscere quale dei 3 muoni in pass_3mu e' quello displaced
      std::vector<pat::TriggerObjectStandAlone> displaced;
      int flag=0;
      //      std::cout<<"pass_3mu size="<<pass_3mu.size()<<std::endl;    
      for(unsigned int i=0;i<pass_3mu.size() ;i++){
	flag=0;
	for(unsigned int j=0;j<pass_jpsi.size();j++){
	  if(pass_3mu[i].pt()==pass_jpsi[j].pt()){
	    flag=1;
	  }
	}
	if(flag==0){
	  displaced.push_back(pass_3mu[i]);
	  //std::cout<<"il mu displaced e'="<<pass_3mu[i].pt()<<std::endl;
	}
      }
      
      
      //      std::cout<<"all muons="<<kaons->size()<<std::endl;
      //Loop  on displaced muons    
      for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx) {
	edm::Ptr<pat::Muon> k_ptr(kaons, k_idx);
	if( !k_selection_(*k_ptr) ) continue;

	//matching online-offline del mu displaced                   
	int flag=0;                 
	for(pat::TriggerObjectStandAlone obj: displaced){
	  if(deltaR(obj,*k_ptr)<0.02){
	    //std::cout<<"deltaR"<<deltaR(obj,*k_ptr)<<std::endl;
	    flag=1;
	  }
	}
	
	if(flag==0){
	  continue;
	}
	else{ //ha trovato il mu displaced
    
	  //std::cout<<"mu displaced triggerato="<<k_ptr->pt()<<std::endl;
	  math::PtEtaPhiMLorentzVector k_p4(
					    k_ptr->pt(), 
					    k_ptr->eta(),
					    k_ptr->phi(),
					    k_ptr->mass()
					    );

	  for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
	    edm::Ptr<pat::CompositeCandidate> ll_prt(dileptons, ll_idx);
	    edm::Ptr<reco::Candidate> l1_ptr = ll_prt->userCand("l1");
	    edm::Ptr<reco::Candidate> l2_ptr = ll_prt->userCand("l2");
	    int l1_idx = ll_prt->userInt("l1_idx");
	    int l2_idx = ll_prt->userInt("l2_idx");

	    //matching online-offline di jpsi                 
	    std::vector<pat::TriggerObjectStandAlone> trig_obj_jpsi;
	    std::vector<edm::Ptr<reco::Candidate>> muon_jpsi;
	    flag=0;
	    for(pat::TriggerObjectStandAlone obj: pass_jpsi){
	      for(pat::TriggerObjectStandAlone obj2: pass_jpsi){  
		if((deltaR(obj,*l1_ptr)<0.02) && (deltaR(obj2,*l2_ptr)<0.02)){
		  flag=1;
		}
	      }
	    }
	    
	    if(flag==0) continue;
	    else{

	      //std::cout<<"dileptons: 1="<<l1_ptr->pt()<<"; 2= "<<l2_ptr->pt()<<std::endl;
	      pat::CompositeCandidate cand;
	      cand.setP4(ll_prt->p4() + k_p4);
	      cand.setCharge(ll_prt->charge() + k_ptr->charge());
	      // Use UserCands as they should not use memory but keep the Ptr itself
	      // Put the lepton passing the corresponding selection
	      //std::cout<<cand.pt()<<" "<<std::endl;
	      cand.addUserCand("l1", l1_ptr);
	      cand.addUserCand("l2", l2_ptr);
	      cand.addUserCand("K", k_ptr);
	      cand.addUserCand("dilepton", ll_prt);
	      
	      cand.addUserInt("l1_idx", l1_idx);
	      cand.addUserInt("l2_idx", l2_idx);
	      cand.addUserInt("k_idx", k_idx);
	      
	      auto dr_info = min_max_dr({l1_ptr, l2_ptr, k_ptr});
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
	      
	      if( !pre_vtx_selection_(cand) ) {
		//std::cout<<"pre vtx selection dies"<<std::endl;
		continue;
	      }
	      
	      //	      std::cout<<"PRIMA"<<std::endl;
	      KinVtxFitter fitter(
				  {leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx), kaons_ttracks->at(k_idx)},
				  {l1_ptr->mass(), l2_ptr->mass(), k_ptr->mass()},
				  {LEP_SIGMA, LEP_SIGMA, LEP_SIGMA} //some small sigma for the lepton mass
				  );
	      //std::cout<<"DOPO"<<std::endl;
	      if(!fitter.success()) {
		//std::cout<<"fitter success dies"<<std::endl;
		continue; // hardcoded, but do we need otherwise?
	      }
		cand.setVertex( 
			     reco::Candidate::Point( 
						    fitter.fitted_vtx().x(),
						    fitter.fitted_vtx().y(),
						    fitter.fitted_vtx().z()
						     )  
			      );
	      used_lep1_id.emplace_back(l1_idx);
	      used_lep2_id.emplace_back(l2_idx);
	      used_trk_id.emplace_back(k_idx);
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
	      cand.addUserFloat("fitted_k_pt"  , fitter.daughter_p4(2).pt()); 
	      cand.addUserFloat("fitted_k_eta" , fitter.daughter_p4(2).eta());
	      cand.addUserFloat("fitted_k_phi" , fitter.daughter_p4(2).phi());
	      
	      //mie variabili                                                                                
	      //conto quanti mu totali aveva l'evento
	      cand.addUserInt("pass_3mu",pass_3mu.size());                                 
	      
	      
	      float B_pt=(Bc_MASS/cand.mass())*cand.pt();
	      
	      TLorentzVector P_b;
	      P_b.SetPtEtaPhiM(B_pt,cand.eta(),cand.phi(),Bc_MASS);
	      
	      TLorentzVector P_k;
	      P_k.SetPtEtaPhiM(k_ptr->pt(),k_ptr->eta(),k_ptr->phi(),k_ptr->mass());
	      
	      TLorentzVector P_l1;
	      P_l1.SetPtEtaPhiM(l1_ptr->pt(),l1_ptr->eta(),l1_ptr->phi(),l1_ptr->mass());
	      
	      TLorentzVector P_l2;
	      P_l2.SetPtEtaPhiM(l2_ptr->pt(),l2_ptr->eta(),l2_ptr->phi(),l2_ptr->mass());
	      
	      
	      float m_miss_2=(P_b-P_k-P_l1-P_l2)*(P_b-P_k-P_l1-P_l2);
	      
	      float Q_2=(P_b-P_l1-P_l2)*(P_b-P_l1-P_l2);
	      
	      float pt_miss=(P_b.Pt()-P_k.Pt()-P_l1.Pt()-P_l2.Pt());
	      
	      float pt_miss_vec=((P_b-P_k-P_l1-P_l2).Pt());

	      float pt_var=((P_l1+P_l2).Pt()-P_k.Pt());
	     
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
	      TLorentzVector P_mu=P_k;	      

	      TVector3 mu_beta_lab=P_b.BoostVector();
	      
	      P_mu.Boost(-mu_beta_lab);
	      
	      cand.addUserFloat("E_mu_star",P_mu.E());
	      
	      P_mu=P_k;	      
	      	      
	      TLorentzVector jpsi=P_l1+P_l2;
	      
	      TVector3 jpsi_beta_lab=jpsi.BoostVector();
	      
	      P_mu.Boost(-jpsi_beta_lab);
	    
	      cand.addUserFloat("E_mu_#",P_mu.E());
	      
	      //std::cout<<cos_theta_2D(fitter, *beamspot, fit_p4)<<std::endl;
	      if( !post_vtx_selection_(cand) ) {
		//std::cout<<"post vtx dies"<<std::endl;
		continue;        
	      }
	      //std::cout<<"EHIIII"<<std::endl;

	      //compute isolation
	      float l1_iso03 = 0;
	      float l1_iso04 = 0;
	      float l2_iso03 = 0;
	      float l2_iso04 = 0;
	      float k_iso03  = 0;
	      float k_iso04  = 0;
	      float b_iso03  = 0;
	      float b_iso04  = 0;
	      
	      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
		
		const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
		// define selections for iso tracks (pT, eta, ...)
		if( !isotrk_selection_(trk) ) 
		  continue;
		// check if the track is the kaon
		/*if (k_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) {
		  
		  std::cout<<"old"<<std::endl;
		  continue;
		}*/
		
		if(track_to_lepton_match(k_ptr, iso_tracks.id(), iTrk)  ) {
		  // std::cout<<"new"<<std::endl;
		  continue;
		}
		// check if the track is one of the two leptons
		if (track_to_lepton_match(l1_ptr, iso_tracks.id(), iTrk) || 
		    track_to_lepton_match(l2_ptr, iso_tracks.id(), iTrk) ) continue;
		
		// add to final particle iso if dR < cone
		float dr_to_l1 = deltaR(cand.userFloat("fitted_l1_eta"), cand.userFloat("fitted_l1_phi"), trk.eta(), trk.phi());
		float dr_to_l2 = deltaR(cand.userFloat("fitted_l2_eta"), cand.userFloat("fitted_l2_phi"), trk.eta(), trk.phi());
		float dr_to_k  = deltaR(cand.userFloat("fitted_k_eta") , cand.userFloat("fitted_k_phi") , trk.eta(), trk.phi());
		float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());
		
		if (dr_to_l1 < 0.4){
		  l1_iso04 += trk.pt();
		  if ( dr_to_l1 < 0.3) l1_iso03 += trk.pt();
		}
		if (dr_to_l2 < 0.4){
		  l2_iso04 += trk.pt();
		  if (dr_to_l2 < 0.3)  l2_iso03 += trk.pt();
		}
		if (dr_to_k < 0.4){
		  k_iso04 += trk.pt();
		  if (dr_to_k < 0.3) k_iso03 += trk.pt();
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
	      cand.addUserFloat("k_iso03" , k_iso03 );
	      cand.addUserFloat("k_iso04" , k_iso04 );
	      cand.addUserFloat("b_iso03" , b_iso03 );
	      cand.addUserFloat("b_iso04" , b_iso04 );

	      //std::cout<<cand.pt()<<" cand"<<std::endl;	      
	      //salvo il candidato
	      ret_val->push_back(cand);
	      //	      std::cout<<ret_val->size()<<" ret"<<std::endl;
	    } //trigger dileptons 
	  }//for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {

	  
	} //trigger muon displaced 
      }//for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx)
      
   


	  
      for (auto & cand: *ret_val){
	cand.addUserInt("n_k_used", std::count(used_trk_id.begin(),used_trk_id.end(),cand.userInt("k_idx")));
	cand.addUserInt("n_l1_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l1_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l1_idx")));
	cand.addUserInt("n_l2_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l2_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l2_idx")));
      }
    }//pass_hlt     
      //    std::cout<<ret_val->size()<<std::endl;
    evt.put(std::move(ret_val));
      
  
  } //index ok
  
}//produce	
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BTommmBuilder);
