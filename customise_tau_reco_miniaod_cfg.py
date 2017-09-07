import FWCore.ParameterSet.Config as cms

# trick to make python know where to look for the imports
import sys
sys.path.append('..')

# for example: here
from tau_miniaod import process

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 200 )
)

# change the isolation cone
process.combinatoricRecoTaus.builders[0].isolationConeSize = cms.double(0.8) # originally 0.5

# change the output file name, don't overwrite the original file!
process.output.fileName = cms.untracked.string('outputFULL_isoCone0p8.root')
