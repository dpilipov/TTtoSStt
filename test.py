from TTClass import TTClass
from TIMBER.Tools.Common import CompileCpp

CompileCpp('TTmodules.cc')
test = TTClass('raw_nano/StoAA_16.txt',16,1,1)

test.analysis1()
test.Snapshot()

