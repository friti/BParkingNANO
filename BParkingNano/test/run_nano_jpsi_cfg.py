from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('isMC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)

options.setDefault('maxEvents',10000)
options.setDefault('tag', '10215')
options.parseArguments()

globaltag = '106X_dataRun2_v28' if not options.isMC else '106X_upgrade2018_realistic_v11_L1v1'
if options._beenSet['globalTag']:
    globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['RJPsi', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['BParkFullEvt', extension[options.isMC], options.tag])+'.root')
if not options.inputFiles:
    options.inputFiles = ['root://cms-xrd-global.cern.ch//store/data/Run2018D/Charmonium/MINIAOD/12Nov2019_UL2018-v1/00000/00AC5CFC-5390-F947-911B-2074A7DFF23D.root'] if not options.isMC else \
                         ['file:/pnfs/psi.ch/cms/trivcat/store/user/manzoni/RJPsi_Bc_PMX_HLT_RECO_MINI_28oct20_v5/RJpsi-BcToXToJpsiMuMuSelected-RunIISummer19UL18MiniAOD_1000.root']                         

annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BParkNANO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.BParkingNano.nanoBPark_mmm_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileFEVT,
    outputCommands = (cms.untracked.vstring('keep *',
                                            'drop *_*_SelectedTransient*_*',
                     )),
    splitLevel = cms.untracked.int32(0)
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileNANO,
    outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')


from PhysicsTools.BParkingNano.nanoBPark_mmm_cff import *
process = nanoAOD_customizeMuonTriggerBPark(process)
process = nanoAOD_customizeTrackFilteredBPark(process)
process = nanoAOD_customizeBTommm(process)
process = nanoAOD_customizeBTopmm(process) #Bc-> JpsiPi -> mu mu pi
process = nanoAOD_customizeBTokmm(process) 
process = nanoAOD_customizeTriggerBitsBPark(process)

# Disable implicit track selection based on the proximity to a BPark trigger muon 
# (that may very well not be there, since we use different triggers)
tracksBPark.drTrg_Cleaning = cms.double(-1.)
tracksBPark.dzTrg_cleaning = cms.double(-1.)
tracksBPark.trkPtCut = cms.double(2.) # adjust track pt

# Path and EndPath definitions
process.nanoAOD_mmm_step = cms.Path(process.nanoSequence + process.nanoBmmmSequence + CountBTommm )
process.nanoAOD_pmm_step = cms.Path(process.nanoSequence + process.nanoBpmmSequence + CountBTopmm )
process.nanoAOD_kmm_step = cms.Path(process.nanoSequence + process.nanoBkmmSequence + CountBTokmm )

# customisation of the process.
if options.isMC:
    from PhysicsTools.BParkingNano.nanoBPark_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
    process.nanoAOD_mmm_step,
    process.nanoAOD_pmm_step,
    process.nanoAOD_kmm_step,
    process.endjob_step, 
    process.NANOAODoutput_step
)
if options.wantFullRECO:
    process.schedule = cms.Schedule(
        process.nanoAOD_mmm_step,
        process.nanoAOD_Kee_step, 
        process.nanoAOD_KstarMuMu_step,
        process.nanoAOD_KstarEE_step,
        process.endjob_step, 
        process.FEVTDEBUGHLToutput_step, 
        process.NANOAODoutput_step
    )
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring(
        'nanoAOD_mmm_step', 
        'nanoAOD_pmm_step',
        'nanoAOD_kmm_step',
    )
)


### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
