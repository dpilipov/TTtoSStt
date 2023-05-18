import ROOT
from TIMBER.Analyzer import Correction, CutGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from helpers import SplitUp
import TIMBER.Tools.AutoJME as AutoJME

AutoJME.AK8collection = 'Dijet'

class TTtoSStt:
    def __init__(self,inputfile,year,ijob,njobs):
        if inputfile.endswith('.txt'): 
            infiles = SplitUp(inputfile, njobs)[ijob-1]
	    # there is an issue with empty events TTrees, so make sure they don't make it through to the analyzer (mainly seen in V+Jets, esp at low HT)
	    invalidFiles = []
	    for iFile in infiles:
		#print('Adding {} to Analyzer'.format(iFile))
		f = ROOT.TFile.Open(iFile)
		if not f.Get('Events'):
		    print('\tWARNING: {} has no Events TTree - will not be added to analyzer'.format(iFile))
		    invalidFiles.append(iFile)
		    continue
		if not f.Get('Events').GetEntry():
		    print('\tWARNING: {} has an empty Events TTree - will not be added to analyzer'.format(iFile))
		    invalidFiles.append(iFile)
		f.Close()
	    inputFiles = [i for i in infiles if i not in invalidFiles]
	    if len(inputFiles) == 0:
		print("\n\t WARNING: None of the files given contain an Events TTree.")
	    self.a = analyzer(inputFiles)
        else:
            infiles = inputfile
            self.a = analyzer(infiles)

        if inputfile.endswith('.txt'):
            self.setname = inputfile.split('/')[-1].split('_')[0]
        else:
            self.setname = inputfile.split('/')[-1].split('_')[1]
        self.year = str(year)	# most of the time this class will be instantiated from other scripts with CLI args, so just convert to string for internal use
        self.ijob = ijob
        self.njobs = njobs
        self.config = OpenJSON('TTconfig.json')
	self.xsec = self.config['XSECS'][self.setname]
        self.trigs = {
        }

        if 'Data' in inputfile:		# SingleMuonDataX_year and DataX_year are possible data inputfile names
            self.a.isData = True
        else:
            self.a.isData = False

    def AddCutflowColumn(self, var, varName):
	'''
	for future reference:
	https://root-forum.cern.ch/t/rdataframe-define-column-of-same-constant-value/34851
	'''
	print('Adding cutflow information...\n\t{}\t{}'.format(varName, var))
	self.a.Define('{}'.format(varName),str(var))
	

    def getNweighted(self):
	if not self.a.isData:
	    return self.a.DataFrame.Sum("genWeight").GetValue()
	else:
	    return self.a.DataFrame.Count().GetValue()

    def analysis1(self):
        
