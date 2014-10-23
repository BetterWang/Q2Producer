import FWCore.ParameterSet.Config as cms

process = cms.Process("QWQ2")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.StandardSequences.SkimsHeavyIons_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_R_53_LV6::All', '')
#process.GlobalTag.globaltag='GR_R_53_LV6::All'

# Minimum bias trigger selection (later runs)
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltMinBiasHFOrBSC = process.hltHighLevel.clone()
process.hltMinBiasHFOrBSC.HLTPaths = ["HLT_HIMinBiasHfOrBSC_v1"]

# Common stuff if you haven't already loaded it
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
#process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")
#process.path = cms.Path(process.nameOfHltSelection*process.collisionEventSelection*process.nameOfAnalysisSequence)

# Coincidence of HF towers (threshold 3)
process.load("HeavyIonsAnalysis.Configuration.hfCoincFilter_cff")
#process.filter_step = cms.Path(process.hfCoincFilter3)

# selection of non-fake vertex (i.e. at least one pixel track)
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
                                           src = cms.InputTag("hiSelectedVertex"),
                                           cut = cms.string("!isFake && abs(z) <= 10 && position.Rho <= 2"),
                                           filter = cms.bool(True),   # otherwise it won't filter the events, instead making an empty vertex collection
                                           )

# Cluster-shape filter re-run offline
process.load("HLTrigger.special.hltPixelClusterShapeFilter_cfi")
process.hltPixelClusterShapeFilter.inputTag = "siPixelRecHits"
#process.filter_step = cms.Path(process.hfCoincFilter3*process.siPixelRecHits*process.hltPixelClusterShapeFilter)

# Reject BSC beam halo L1 technical bits
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff")
process.load("HLTrigger.HLTfilters.hltLevel1GTSeed_cfi")
process.noBSChalo = process.hltLevel1GTSeed.clone(
        L1TechTriggerSeeding = cms.bool(True),
        L1SeedsLogicalExpression = cms.string('NOT (36 OR 37 OR 38 OR 39)')
        )

process.filter_step = cms.Path(process.hltMinBiasHFOrBSC*process.hfCoincFilter3*process.siPixelRecHits*process.hltPixelClusterShapeFilter*process.noBSChalo*process.primaryVertexFilter)

#process.source = cms.Source("PoolSource",
 #                           fileNames = cms.untracked.vstring('/store/hidata/HIRun2011/HIMinBiasUPC/RECO/14Mar2014-v2/00000/0018A8E7-F9AF-E311-ADAB-FA163E565820.root')
  #                          )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/opt/data/2011/3EB2B76A-66B0-E311-B10A-FA163ED551F2.root'
        )
#    fileNames = cms.untracked.vstring(
#'/store/user/appeltel/HIMinBiasUPC/pixelTrackReco_devel_v0/a236e4501225ae15b3601563d612abb5/pixeltrackreco_1_1_sf8.root'
#  ),
#    secondaryFileNames = cms.untracked.vstring(
#'/store/hidata/HIRun2011/HIMinBiasUPC/RECO/14Mar2014-v2/00030/784B6517-F9B3-E311-B99A-FA163E54935D.root'
#  )
)

process.output = cms.OutputModule("PoolOutputModule",
                                  outputCommands = process.RECOEventContent.outputCommands,
                                  #outputCommands={ look up keep/drop
                                  fileName = cms.untracked.string('HIMinBias2011_umdskim.root'),
                                  SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('filter_step')),
                                  dataset = cms.untracked.PSet(
    dataTier = cms.untracked.string('RECO'))

                                  )
#Deciding what to keep

#Drop everything
process.output.outputCommands = ['drop *_*_*_*']

##Keeping tracks and pixel tracks
#process.output.outputCommands += ['keep *_hiGeneralTracks_*_*']
#process.output.outputCommands += ['keep *_hiConformalPixelTracks_*_*']
process.output.outputCommands +=['keep *_hiLowPtPixelTracks_*_*']
##Keep Vertex Info
process.output.outputCommands += ['keep *_hiSelectedVertex_*_*']
##Keep CaloTowers
process.output.outputCommands += ['keep *_towerMaker_*_*']
##Keep ZDC Digis
process.output.outputCommands += ['keep ZDCDataFramesSorted_hcalDigis_*_*']
##Keep Castor RecHits
process.output.outputCommands += ['keep *_castorreco_*_*']
##Keep Centrality
process.output.outputCommands += ['keep *_hiCentrality_*_*']

#Finally load the output step
#process.output_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(process.filter_step
#                            process.output_step
                            )

