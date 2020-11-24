import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.BParkingNano.common_cff import BParkCandVars, ufloat, uint, ubool

JpsiMuonPairs = cms.EDProducer(
    'DiMuonBuilder',
    src                = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection      = cms.string('pt > 2.5'),
    lep2Selection      = cms.string(''),
    preVtxSelection    = cms.string(' && '.join([
#         'abs(userCand("l1").dz - userCand("l2").dz) <= 0.4 ',
        'abs(userCand("l1").bestTrack.dz - userCand("l2").bestTrack.dz) <= 0.4 ',
        'mass() > 2',
        'mass() < 4',
        'userFloat("lep_deltaR") > 0.01',
        ])
    ),
    postVtxSelection   = cms.string(
        'pt > 3 '
#         '&& userFloat("sv_chi2") < 998 ' 
        '&& userFloat("sv_prob") > 1.e-5 '
    ),
)

BuilderDefaultCfg = cms.PSet(
    dileptons             = cms.InputTag('JpsiMuonPairs'),
    leptonTransientTracks = JpsiMuonPairs.transientTracksSrc,
    kaons                 = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks  = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot              = cms.InputTag("offlineBeamSpot"),
    tracks                = cms.InputTag("packedPFCandidates"),
    lostTracks            = cms.InputTag("lostTracks"),
    kaonSelection         = cms.string(''),
    isoTracksSelection    = cms.string('pt > 0.5 && abs(eta)<2.5'),
    preVtxSelection       = cms.string(''),
    postVtxSelection      = cms.string(' && '.join([
        'userInt("sv_OK") == 1',
        'userFloat("sv_prob") > 1e-8',
        'userFloat("fitted_cos_theta_2D") >= 0',
#         'userFloat("fitted_mass") > 4.5',
#         'userFloat("fitted_mass") < 8.',
        'mass > 4.5',
        'mass < 8.',
        ])
    ),
    bits                  = cms.InputTag("TriggerResults","","HLT"),               
    objects               = cms.InputTag("slimmedPatTrigger"), 
)

TableDefault = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("your_cand_name"),
    cut       = cms.string(""),
    name      = cms.string("your_cand_name"),
    doc       = cms.string("your_cand_name Variable"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(),
)

TableDefaultVariables = cms.PSet(
    # pre-fit quantities                                                      
    BParkCandVars,
    #nome branch= nome variabile del .cc
    l1Idx    = uint('l1_idx'),
    l2Idx    = uint('l2_idx'),
    kIdx     = uint('k_idx'),
    minDR    = ufloat('min_dr'),
    maxDR    = ufloat('max_dr'),
    # fit and vtx info                                                                                                    
    chi2     = ufloat('sv_chi2'),            
    svprob   = ufloat('sv_prob'),
    l_xy     = ufloat('l_xy'),
    l_xy_unc = ufloat('l_xy_unc'),
    vtx_x    = ufloat('vtx_x'),
    vtx_y    = ufloat('vtx_y'),
    vtx_z    = ufloat('vtx_z'),
    vtx_ex   = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix                                         
    vtx_ey   = ufloat('vtx_ey'),
    vtx_ez   = ufloat('vtx_ez'),
    # Mll                                                                                                                 
    mll_raw      = Var('userCand("dilepton").mass()', float),
    mll_llfit    = Var('userCand("dilepton").userFloat("fitted_mass")', float), # this might not work                        
    mllErr_llfit = Var('userCand("dilepton").userFloat("fitted_massErr")', float), # this might not work                  
    mll_fullfit  = ufloat('fitted_mll'),
    # Cos(theta)                                                                                                          
    cos2D     = ufloat('cos_theta_2D'),
    fit_cos2D = ufloat('fitted_cos_theta_2D'),
    # post-fit momentum                                                                                                   
    fit_mass    = ufloat('fitted_mass'),
    fit_massErr = ufloat('fitted_massErr'),
    fit_pt      = ufloat('fitted_pt'),
    fit_eta     = ufloat('fitted_eta'),
    fit_phi     = ufloat('fitted_phi'),
    fit_l1_pt   = ufloat('fitted_l1_pt'),
    fit_l1_eta  = ufloat('fitted_l1_eta'),
    fit_l1_phi  = ufloat('fitted_l1_phi'),
    fit_l2_pt   = ufloat('fitted_l2_pt'),
    fit_l2_eta  = ufloat('fitted_l2_eta'),
    fit_l2_phi  = ufloat('fitted_l2_phi'),
    fit_k_pt    = ufloat('fitted_k_pt'),
    fit_k_eta   = ufloat('fitted_k_eta'),
    fit_k_phi   = ufloat('fitted_k_phi'),
    l1_iso03    = ufloat('l1_iso03'),
    l1_iso04    = ufloat('l1_iso04'),
    l2_iso03    = ufloat('l2_iso03'),
    l2_iso04    = ufloat('l2_iso04'),
    k_iso03     = ufloat('k_iso03'),
    k_iso04     = ufloat('k_iso04'),
    b_iso03     = ufloat('b_iso03'),
    b_iso04     = ufloat('b_iso04'),

    # my variables
    m_miss_sq   = ufloat('m_miss_2'),
    Q_sq        = ufloat('Q_2'),
    pt_miss     = ufloat('pt_miss'),
    pt_miss_vec = ufloat('pt_miss_vec'),
    pt_var      = ufloat('pt_var'),
    DR          = ufloat('DR'),
    E_mu_star   = ufloat('E_mu_star'),
    E_mu_canc   = ufloat('E_mu_#'),
    m_jpsi      = ufloat('m_jpsi'),

    # Gen Variables
    is_jpsi_mu   = uint("is_jpsi_mu"),
    is_psi2s_mu  = uint("is_psi2s_mu"),
    is_chic0_mu  = uint("is_chic0_mu"),
    is_chic1_mu  = uint("is_chic1_mu"),
    is_chic2_mu  = uint("is_chic2_mu"),
    is_hc_mu     = uint("is_hc_mu"),
    is_jpsi_tau  = uint("is_jpsi_tau"),
    is_psi2s_tau = uint("is_psi2s_tau"),
    is_jpsi_pi   = uint("is_jpsi_pi"),
    is_jpsi_3pi  = uint("is_jpsi_3pi"),
    is_jpsi_hc   = uint("is_jpsi_hc"),
    is_error     = uint("is_error"),
    weightGen    = ufloat("weightGen"),
)
