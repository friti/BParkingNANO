import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars

def ufloat(expr, precision = -1, doc = ''):
  return Var('userFloat("%s")' % expr, 
             float, precision = precision, doc = doc)

def uint(expr, doc = ''):
  return Var('userInt("%s")' % expr, int, doc = doc)

def ubool(expr, precision = -1, doc = ''):
  return Var('userInt("%s") == 1' % expr, bool, doc = doc)

# Explicitly set the desired precision
BParkCandVars = CandVars.clone()
BParkCandVars.charge.precision = cms.int32(-1)
BParkCandVars.eta   .precision = cms.int32(-1)
BParkCandVars.mass  .precision = cms.int32(-1)
BParkCandVars.pdgId .precision = cms.int32(-1)
BParkCandVars.phi   .precision = cms.int32(-1)
BParkCandVars.pt    .precision = cms.int32(-1)

# CandVars = cms.PSet(
#     charge = cms.PSet(
#         compression = cms.string('none'),
#         doc = cms.string('electric charge'),
#         expr = cms.string('charge'),
#         mcOnly = cms.bool(False),
#         precision = cms.int32(-1),
#         type = cms.string('int')
#     ),
#     eta = cms.PSet(
#         compression = cms.string('none'),
#         doc = cms.string('eta'),
#         expr = cms.string('eta'),
#         mcOnly = cms.bool(False),
#         precision = cms.int32(12),
#         type = cms.string('float')
#     ),
#     mass = cms.PSet(
#         compression = cms.string('none'),
#         doc = cms.string('mass'),
#         expr = cms.string('mass'),
#         mcOnly = cms.bool(False),
#         precision = cms.int32(10),
#         type = cms.string('float')
#     ),
#     pdgId = cms.PSet(
#         compression = cms.string('none'),
#         doc = cms.string('PDG code assigned by the event reconstruction (not by MC truth)'),
#         expr = cms.string('pdgId'),
#         mcOnly = cms.bool(False),
#         precision = cms.int32(-1),
#         type = cms.string('int')
#     ),
#     phi = cms.PSet(
#         compression = cms.string('none'),
#         doc = cms.string('phi'),
#         expr = cms.string('phi'),
#         mcOnly = cms.bool(False),
#         precision = cms.int32(12),
#         type = cms.string('float')
#     ),
#     pt = cms.PSet(
#         compression = cms.string('none'),
#         doc = cms.string('pt'),
#         expr = cms.string('pt'),
#         mcOnly = cms.bool(False),
#         precision = cms.int32(-1),
#         type = cms.string('float')
#     )
# )

