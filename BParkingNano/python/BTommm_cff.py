import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.BParkingNano.common_cff import BParkCandVars, ufloat, uint, ubool
from PhysicsTools.BParkingNano.rjpsi_common_cff import JpsiMuonPairs, BuilderDefaultCfg, TableDefaultVariables, TableDefault,Final3PartTableVariables

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

BTommmTableVariables = Final3PartTableVariables.clone()

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
