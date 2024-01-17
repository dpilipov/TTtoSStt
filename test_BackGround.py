from TTClass_BackGround import TTClass_BackGround
from TIMBER.Tools.Common import CompileCpp

CompileCpp('TTmodules.cc')
test_BackGround = TTClass_BackGround('raw_nano/StoAA_16.txt',16,1,1)
print("At 1")
test_BackGround.analysis1()
print("At 2")
test_BackGround.Snapshot()

