import FWCore.ParameterSet.Config as cms

HiQ2 = cms.EDProducer('Q2Producer'
                      vtxCollection_ = cms.InputTag("hiSelectedVertex"),
                      caloCollection_ = cms.InputTag("towerMaker"),
                      trackCollection_ = cms.InputTag("hiGeneralAndPixelTracks")
)
