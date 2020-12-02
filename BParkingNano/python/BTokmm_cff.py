import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.BParkingNano.common_cff import BParkCandVars, ufloat, uint, ubool
from PhysicsTools.BParkingNano.rjpsi_common_cff import JpsiMuonPairs, BuilderDefaultCfg, TableDefaultVariables, TableDefault

BTokmmCfg = BuilderDefaultCfg.clone()
BTokmmCfg.dileptons = cms.InputTag('JpsiMuonPairs')
BTokmmCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc
BTokmmCfg.postVtxSelection = cms.string(' && '.join([
        BuilderDefaultCfg.postVtxSelection.value(),
        'mass > 4.5',
        ])
)

BTokmm = cms.EDProducer(
    'BToTrackmmBuilder',
    BTokmmCfg,
    srcGen = cms.InputTag("prunedGenParticles"),
    track_mass            = cms.double(0.493677)
)

BTokmmTableVariables = TableDefaultVariables.clone()

BTokmmTable = TableDefault.clone()
BTokmmTable.src       = cms.InputTag("BTokmm")
BTokmmTable.name      = cms.string("BTokmm")
BTokmmTable.doc       = cms.string("BTokmm Variable")
BTokmmTable.variables = BTokmmTableVariables

BTokmmSequence = cms.Sequence(
    (JpsiMuonPairs * BTokmm)
)


CountBTokmm = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BTokmm")
)    

