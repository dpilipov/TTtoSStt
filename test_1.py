from TTtoSStry1 import TTtoSStt
from TIMBER.Tools.Common import CompileCpp

CompileCpp('TTmodules.cc')
test = TTtoSStt('raw_nano/StoAA_16.txt',16,1,1)

test.analysis1()
test.Snapshot()

