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
        self.year = str(year)   # most of the time this class will be instantiated from other scripts with CLI args, so just convert to string for internal use
        self.ijob = ijob
        self.njobs = njobs
        self.config = OpenJSON('TTconfig.json')
        self.xsec = self.config['XSECS'][self.setname]
        self.trigs = {
        }

        if 'Data' in inputfile:         # SingleMuonDataX_year and DataX_year are possible data inputfile names
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

    def Preselection(self):
        # Apply recommended flags
        flags = [
            'Flag_goodVertices',
            'Flag_globalSuperTightHalo2016Filter',
            'Flag_HBHENoiseFilter',
            'Flag_HBHENoiseIsoFilter',
            'Flag_EcalDeadCellTriggerPrimitiveFilter',
            'Flag_BadPFMuonFilter',
            'Flag_BadPFMuonDzFilter',
            'Flag_eeBadScFilter'
        ]
        if self.year == '17' or self.year == '18':
            flags.append('Flag_ecalBadCalibFilter')
        # string valid (existing in RDataFrame node) flags together w logical and
        MET_filters = self.a.GetFlagString(flags)
        self.a.Cut('flags', MET_filters)

        # Expecting at least 2 jets in the event
        self.a.Cut('njets','nFatJet >= 2')

        # Expecting at least 2 photons in the event
        self.a.Cut('nPhoton','nPhoton >= 2')

        # INFO: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#nanoAOD_Flags
        # drop any events whose dijets did not both pass tight jetId requirement
        self.a.Cut('jetId', 'Jet_jetId[0] > 1 && Jet_jetId[1] > 1')

        # just do a basic check for one of the jets and associated photons
        # The function PickJetAnd2Gamma() is a custom C++ code
        self.a.Define('JetAndDiphotonIdxs','PickJetAnd2Gamma(FatJet_pt, FatJet_eta, FatJet_phi, FatJet_msoftdrop, Photon_phi)')

        # cut out all events which didn't meet our criteria
        self.a.Cut('validEvents','JetAndDiphotonIdxs[0] > 0 && JetAndDiphotonIdxs[1] > 0 && JetAndDiphotonIdxs[2] > 0')

        # now create subcollections 
        self.a.ObjectFromCollection('JetCand','FatJet','JetAndDiphotonIdxs[0]',skip=['HLT_*'])
        self.a.ObjectFromCollection('PhoCand1','Photon','JetAndDiphotonIdxs[1]',skip=['HLT_*'])
        self.a.ObjectFromCollection('PhoCand2','Photon','JetAndDiphotonIdxs[2]',skip=['HLT_*'])

        # construct TLVectors
        self.a.Define('JetCand_vect','hardware::TLvector(JetCand_pt, JetCand_eta, JetCand_phi, JetCand_msoftdrop)')
        self.a.Define('PhoCand1_vect','hardware::TLvector(PhoCand1_pt, PhoCand1_eta, PhoCand1_phi, PhoCand1_mass)')
        self.a.Define('PhoCand2_vect','hardware::TLvector(PhoCand2_pt, PhoCand2_eta, PhoCand2_phi, PhoCand2_mass)')

        # now make invariant masses
        self.a.Define('mTprime','hardware::InvariantMass({JetCand_vect, PhoCand1_vect, PhoCand2_vect})')
        self.a.Define('mScalar','hardware::InvariantMass({PhoCand1_vect, PhoCand2_vect})')

        # now save it all out
        self.a.Snapshot(['mTprime','mScalar','JetAndDiphotonIdxs'],'{}_{}_snapshot.root'.format(self.setname,self.year),'Events',openOption='RECREATE',saveRunChain=True)

