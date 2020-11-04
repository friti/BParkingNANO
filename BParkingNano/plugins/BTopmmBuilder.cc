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

class BTopmmBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BTopmmBuilder(const edm::ParameterSet &cfg):
    k_selection_{cfg.getParameter<std::string>("kaonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    kaons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("kaons") )},
    kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),

    //TRIGGER
    triggerBits_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(cfg.getParameter<edm::InputTag>("objects"))), 

    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>();
    }

  ~BTopmmBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> k_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> kaons_;
  const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;

  //TRIGGER                                                              
  const edm::EDGetTokenT<edm::TriggerResults> triggerBits_;                         
  const edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;

  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

void BTopmmBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

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
  
  unsigned int index = names.triggerIndex("HLT_DoubleMu4_JpsiTrk_Displaced_v14");
  unsigned int index2 = names.triggerIndex("HLT_DoubleMu4_JpsiTrk_Displaced_v15");
  bool pass_hlt;
  if(index==triggerBits->size() && index2==triggerBits->size()){
    //    std::cout<<"Non ha HLT path giusto"<<std::endl;
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
      //std::cout<<"entra in pass_hlt"<<std::endl;
      for (pat::TriggerObjectStandAlone obj : *triggerObjects){
	obj.unpackFilterLabels(evt, *triggerBits);
	obj.unpackPathNames(names);
	
	//std::cout<<names<<std::endl;
	if(obj.hasFilterLabel("hltDisplacedmumuFilterDoubleMu4Jpsi")) {
	  pass_jpsi.push_back(obj);
	}
	else if(obj.hasFilterLabel("hltJpsiTkVertexFilter")) {
	  pass_trk.push_back(obj);
	}
	
      }
   

      //loop on trks
      for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx) {
	edm::Ptr<pat::CompositeCandidate> k_ptr(kaons, k_idx);
	if( !k_selection_(*k_ptr) ) continue;
	
	//matching online-offline della trk
	int flag = 0;
	for(pat::TriggerObjectStandAlone obj: pass_trk){
	  if(deltaR(obj,*k_ptr)<0.02) flag=1;
	}    
	if(flag == 0) continue;
	
	else{// found the trk matching the trigger
	  math::PtEtaPhiMLorentzVector k_p4(
					    k_ptr->pt(), 
					    k_ptr->eta(),
					    k_ptr->phi(),
					    PI_MASS
					    );
	  
	  for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
	    edm::Ptr<pat::CompositeCandidate> ll_prt(dileptons, ll_idx);
	    edm::Ptr<reco::Candidate> l1_ptr = ll_prt->userCand("l1");
	    edm::Ptr<reco::Candidate> l2_ptr = ll_prt->userCand("l2");
	    int l1_idx = ll_prt->userInt("l1_idx");
	    int l2_idx = ll_prt->userInt("l2_idx");
	    
	    
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
	      cand.setP4(ll_prt->p4() + k_p4);
	      cand.setCharge(ll_prt->charge() + k_ptr->charge());
	      // Use UserCands as they should not use memory but keep the Ptr itself
	      // Put the lepton passing the corresponding selection
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
	      // TODO add meaningful variables
	      
	      if( !pre_vtx_selection_(cand) ) continue;
	      
	      KinVtxFitter fitter(
				  {leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx), kaons_ttracks->at(k_idx)},
				  {l1_ptr->mass(), l2_ptr->mass(), PI_MASS},
				  {LEP_SIGMA, LEP_SIGMA, PI_SIGMA} //some small sigma for the lepton mass
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
	      

	      //new variables

	      TLorentzVector P_b;
	      P_b.SetPtEtaPhiM(cand.pt(),cand.eta(),cand.phi(),cand.mass());
	      
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
	      

	      if( !post_vtx_selection_(cand) ) continue;        
	      
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
		if( !isotrk_selection_(trk) ) continue;
		// check if the track is the kaon
		if (k_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
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
	      
	      ret_val->push_back(cand);
	    }//trigger dileptons
	  } // for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
	}//trigger trk
      } // for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx)
      
      for (auto & cand: *ret_val){
	cand.addUserInt("n_k_used", std::count(used_trk_id.begin(),used_trk_id.end(),cand.userInt("k_idx")));
	cand.addUserInt("n_l1_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l1_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l1_idx")));
	cand.addUserInt("n_l2_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l2_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l2_idx")));
      }
    }//pass_hlt
    evt.put(std::move(ret_val));
  }//index ok
}//produce

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BTopmmBuilder);
