import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.BParkingNano.common_cff import BParkCandVars, ufloat, uint, ubool
from PhysicsTools.BParkingNano.rjpsi_common_cff import JpsiMuonPairs, BuilderDefaultCfg, TableDefaultVariables, TableDefault

BTommmCfg = BuilderDefaultCfg.clone()
BTommmCfg.dileptons             = cms.InputTag('JpsiMuonPairs')
BTommmCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc
BTommmCfg.kaons                 = cms.InputTag('muonTrgSelector', 'SelectedMuons')
BTommmCfg.kaonsTransientTracks  = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons')

BTommm = cms.EDProducer(
    'BTommmBuilder',
    BTommmCfg,
    srcGen = cms.InputTag("prunedGenParticles"),
)

BTommmTableVariables = TableDefaultVariables.clone()
# Gen Variables
BTommmTableVariables.is_jpsi_mu   = uint("is_jpsi_mu")
BTommmTableVariables.is_psi2s_mu  = uint("is_psi2s_mu")
BTommmTableVariables.is_chic0_mu  = uint("is_chic0_mu")
BTommmTableVariables.is_chic1_mu  = uint("is_chic1_mu")
BTommmTableVariables.is_chic2_mu  = uint("is_chic2_mu")
BTommmTableVariables.is_hc_mu     = uint("is_hc_mu")
BTommmTableVariables.is_jpsi_tau  = uint("is_jpsi_tau")
BTommmTableVariables.is_psi2s_tau = uint("is_psi2s_tau")
BTommmTableVariables.is_jpsi_pi   = uint("is_jpsi_pi")
BTommmTableVariables.is_jpsi_3pi  = uint("is_jpsi_3pi")
BTommmTableVariables.is_jpsi_hc   = uint("is_jpsi_hc")
BTommmTableVariables.is_error     = uint("is_error")
BTommmTableVariables.weightGen    = ufloat("weightGen")

BTommmTable = TableDefault.clone()
BTommmTable.src       = cms.InputTag("BTommm")
BTommmTable.name      = cms.string("BTommm")
BTommmTable.doc       = cms.string("BTommm Variable")
BTommmTable.variables = BTommmTableVariables

BTommmSequence = cms.Sequence(
    (JpsiMuonPairs * BTommm)
)

CountBTommm = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BTommm")
)    
