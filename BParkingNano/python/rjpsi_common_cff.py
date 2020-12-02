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
    beamSpot              = cms.InputTag("offlineBeamSpot"),

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
#         'mass > 4.5',
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


    #3 particles vertex
    bodies3_chi2     = ufloat('sv_chi2'),            
    bodies3_svprob   = ufloat('sv_prob'),
    bodies3_l_xy     = ufloat('l_xy'),
    bodies3_l_xy_unc = ufloat('l_xy_unc'),
    bodies3_vtx_x    = ufloat('vtx_x'),
    bodies3_vtx_y    = ufloat('vtx_y'),
    bodies3_vtx_z    = ufloat('vtx_z'),
    bodies3_vtx_ex   = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix      
    bodies3_vtx_ey   = ufloat('vtx_ey'),
    bodies3_vtx_ez   = ufloat('vtx_ez'),
    bodies3_cos2D     = ufloat('cos_theta_2D'),

    # post-fit 3 particles vertex                                                                                          
    bodies3_fit_mass    = ufloat('fitted_mass'),
    bodies3_fit_massErr = ufloat('fitted_massErr'),
    bodies3_fit_pt      = ufloat('fitted_pt'),
    bodies3_fit_eta     = ufloat('fitted_eta'),
    bodies3_fit_phi     = ufloat('fitted_phi'),
    bodies3_fit_l1_pt   = ufloat('fitted_l1_pt'),
    bodies3_fit_l1_eta  = ufloat('fitted_l1_eta'),
    bodies3_fit_l1_phi  = ufloat('fitted_l1_phi'),
    bodies3_fit_l2_pt   = ufloat('fitted_l2_pt'),
    bodies3_fit_l2_eta  = ufloat('fitted_l2_eta'),
    bodies3_fit_l2_phi  = ufloat('fitted_l2_phi'),
    bodies3_fit_k_pt    = ufloat('fitted_k_pt'),
    bodies3_fit_k_eta   = ufloat('fitted_k_eta'),
    bodies3_fit_k_phi   = ufloat('fitted_k_phi'),
    bodies3_fit_cos2D = ufloat('fitted_cos_theta_2D'),

    #2 particles vertex
    jpsivtx_chi2 = Var('userCand("dilepton").userFloat("sv_chi2")', float),
    jpsivtx_svprob = Var('userCand("dilepton").userFloat("sv_prob")', float),
    jpsivtx_l_xy = Var('userCand("dilepton").userFloat("l_xy")', float),
    jpsivtx_l_xy_unc = Var('userCand("dilepton").userFloat("l_xy_unc")', float),
    jpsivtx_vtx_x = Var('userCand("dilepton").userFloat("vtx_x")', float),
    jpsivtx_vtx_y = Var('userCand("dilepton").userFloat("vtx_y")', float),
    jpsivtx_vtx_z = Var('userCand("dilepton").userFloat("vtx_z")', float),
    jpsivtx_vtx_ex = Var('userCand("dilepton").userFloat("vtx_ex")', float),
    jpsivtx_vtx_ey = Var('userCand("dilepton").userFloat("vtx_ey")', float),
    jpsivtx_vtx_ez = Var('userCand("dilepton").userFloat("vtx_ez")', float),
    jpsivtx_cos2D = Var('userCand("dilepton").userFloat("cos_theta_2D")', float),
    
    #post fit 2 particles vertex
    jpsivtx_fit_mass    = Var('userCand("dilepton").userFloat("fitted_mass")', float), 
    jpsivtx_fit_massErr = Var('userCand("dilepton").userFloat("fitted_massErr")', float), 
    jpsivtx_fit_pt    = Var('userCand("dilepton").userFloat("fitted_pt")', float), 
    jpsivtx_fit_eta    = Var('userCand("dilepton").userFloat("fitted_eta")', float), 
    jpsivtx_fit_phi   = Var('userCand("dilepton").userFloat("fitted_phi")', float), 
    jpsivtx_fit_l1_pt   = Var('userCand("dilepton").userFloat("fitted_l1_pt")', float), 
    jpsivtx_fit_l1_eta   = Var('userCand("dilepton").userFloat("fitted_l1_eta")', float), 
    jpsivtx_fit_l1_phi   = Var('userCand("dilepton").userFloat("fitted_l1_phi")', float), 
    jpsivtx_fit_l2_pt   = Var('userCand("dilepton").userFloat("fitted_l2_pt")', float), 
    jpsivtx_fit_l2_eta   = Var('userCand("dilepton").userFloat("fitted_l2_eta")', float), 
    jpsivtx_fit_l2_phi   = Var('userCand("dilepton").userFloat("fitted_l2_phi")', float), 
    jpsivtx_fit_cos2D = Var('userCand("dilepton").userFloat("fitted_cos_theta_2D")', float),
    

    
    # Mll
    mll_raw      = Var('userCand("dilepton").mass()', float),
    mll_fullfit  = ufloat('fitted_mll'),


    l1_iso03    = ufloat('l1_iso03'),
    l1_iso04    = ufloat('l1_iso04'),
    l2_iso03    = ufloat('l2_iso03'),
    l2_iso04    = ufloat('l2_iso04'),
    k_iso03     = ufloat('k_iso03'),
    k_iso04     = ufloat('k_iso04'),
    b_iso03     = ufloat('b_iso03'),
    b_iso04     = ufloat('b_iso04'),

    #beamspot
    beamspot_x     = ufloat('beamspot_x'),
    beamspot_y     = ufloat('beamspot_y'),
    beamspot_z     = ufloat('beamspot_z'),


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
