import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.BParkingNano.common_cff import BParkCandVars, ufloat, uint, ubool
from PhysicsTools.BParkingNano.rjpsi_common_cff import JpsiMuonPairs, BuilderDefaultCfg, TableDefaultVariables, TableDefault

BTopmmCfg = BuilderDefaultCfg.clone()
BTopmmCfg.dileptons = cms.InputTag('JpsiMuonPairs')
BTopmmCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc

BTopmm = cms.EDProducer(
    'BToTrackmmBuilder',
    BTopmmCfg,
    srcGen = cms.InputTag("prunedGenParticles"),
    track_mass            = cms.double(0.139571),
)

BTopmmTableVariables = TableDefaultVariables.clone()

BTopmmTable = TableDefault.clone()
BTopmmTable.src       = cms.InputTag("BTopmm")
BTopmmTable.name      = cms.string("BTopmm")
BTopmmTable.doc       = cms.string("BTopmm Variable")
BTopmmTable.variables = BTopmmTableVariables

BTopmmSequence = cms.Sequence(
    (JpsiMuonPairs * BTopmm)
)


CountBTopmm = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BTopmm")
)    

