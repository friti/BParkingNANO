import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.BParkingNano.common_cff import BParkCandVars, ufloat, uint, ubool
from PhysicsTools.BParkingNano.rjpsi_common_cff import JpsiMuonPairs, BuilderDefaultCfg, TableDefaultVariables, TableDefault

BTo2mu3piCfg = BuilderDefaultCfg.clone()
BTo2mu3piCfg.dileptons = cms.InputTag('JpsiMuonPairs')
BTo2mu3piCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc
BTo2mu3piCfg.postVtxSelection = cms.string(' && '.join([
        BuilderDefaultCfg.postVtxSelection.value(),
        'mass > 4.5',
        ])
)
BTo2mu3piCfg.kaonSelection = cms.string(' && '.join([
    #    BuilderDefaultCfg.kaonSelection.value(),
    'pt > 2.5',
    'eta < 2.5',
        ])
)

BTo2mu3pi = cms.EDProducer(
    'BTo2mu3piBuilder',
    BTo2mu3piCfg,
    srcGen = cms.InputTag("prunedGenParticles"),
    track_mass            = cms.double(0.493677)
)

BTo2mu3piTableVariables = TableDefaultVariables.clone(
    pi1Idx     = uint('pi1_idx'),
    pi2Idx     = uint('pi2_idx'),
    pi3Idx     = uint('pi3_idx'),
    bodies3_fit_pi1_pt    = ufloat('fitted_pi1_pt'),
    bodies3_fit_pi1_eta   = ufloat('fitted_pi1_eta'),
    bodies3_fit_pi1_phi   = ufloat('fitted_pi1_phi'),
    bodies3_fit_pi2_pt    = ufloat('fitted_pi2_pt'),
    bodies3_fit_pi2_eta   = ufloat('fitted_pi2_eta'),
    bodies3_fit_pi2_phi   = ufloat('fitted_pi2_phi'),
    bodies3_fit_pi3_pt    = ufloat('fitted_pi3_pt'),
    bodies3_fit_pi3_eta   = ufloat('fitted_pi3_eta'),
    bodies3_fit_pi3_phi   = ufloat('fitted_pi3_phi'),
    pi1_iso03     = ufloat('pi1_iso03'),
    pi1_iso04     = ufloat('pi1_iso04'),
    pi2_iso03     = ufloat('pi2_iso03'),
    pi2_iso04     = ufloat('pi2_iso04'),
    pi3_iso03     = ufloat('pi3_iso03'),
    pi3_iso04     = ufloat('pi3_iso04'),

)

BTo2mu3piTable = TableDefault.clone()
BTo2mu3piTable.src       = cms.InputTag("BTo2mu3pi")
BTo2mu3piTable.name      = cms.string("BTo2mu3pi")
BTo2mu3piTable.doc       = cms.string("BTo2mu3pi Variable")
BTo2mu3piTable.variables = BTo2mu3piTableVariables

BTo2mu3piSequence = cms.Sequence(
    (JpsiMuonPairs * BTo2mu3pi)
)


CountBTo2mu3pi = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BTo2mu3pi")
)    

