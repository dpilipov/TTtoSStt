#include "ROOT/RVec.hxx"
#include "TIMBER/Framework/include/common.h"
#include <string>
#include <vector>
#include <random> // because TRandom wouldn't work in this case..

using namespace ROOT::VecOps;

/*****************************************************************
 *  Lepton veto functions 
 *  ---------------------
 *  Loops over all leptons in an event and returns true if 
 *  a lepton in the event meets the criteria for vetoing the 
 *  event. These functions should return a True if the criteria
 *  is met and False if the criteria is not met. The user is 
 *  in charge of correctly interpreting the results of the 
 *  function return, i.e., in an RDataFrame Filter() call.
 *****************************************************************/
bool TightMuVeto(int nMuon, RVec<bool> tightId, RVec<float> muonPt, RVec<float> muonRelIso, RVec<float> muonEta) {
    if (nMuon < 1) {return false;}	// don't veto event, there are no muons
    bool veto = false;
    for (int iMu = 0; iMu < muonPt.size(); iMu++) {
//DP EDIT
//	veto = (tightId[iMu] == 1) && (muonPt[iMu] > 30.) && (muonRelIso[iMu] < 0.15) && (std::abs(muonEta[iMu]) < 2.4); // veto if meets SL selection criteria
      veto = (tightId[iMu] == 1) && (muonPt[iMu] > 15.) && (muonRelIso[iMu] < 0.15) && (std::abs(muonEta[iMu]) < 3.0); 
	if (veto) {return veto;}
    }
    return veto;
} 

bool TightElVeto(int nElectron, RVec<bool> elIso, RVec<float> elPt, RVec<float> elEta) {
    if (nElectron < 1) {return false;}
    bool veto = false;
    for (int iEl = 0; iEl < elPt.size(); iEl++){
//DP EDIT
//	veto = (elIso[iEl] == 1) && (elPt[iEl] > 35.) && (std::abs(elEta[iEl])<2.5);
      veto = (elIso[iEl] == 1) && (elPt[iEl] > 15.) && (std::abs(elEta[iEl])<3.0);
	if (veto) {return veto;}
    }
    return veto;
}

bool GoodMuVeto(int nMuon, RVec<float> muonPt, RVec<bool> looseId, RVec<float> dxy, RVec<float> muonEta) {
    if (nMuon < 1) {return false;}
    bool veto = false;
    for (int iMu = 0; iMu < muonPt.size(); iMu++) {
// DP EDIT
//	veto = (muonPt[iMu] > 30.) && (looseId[iMu] == 1) && (std::abs(dxy[iMu]) < 0.02) && (std::abs(muonEta[iMu]) < 2.4);
      veto = (muonPt[iMu] > 15.) && (looseId[iMu] == 1) && (std::abs(dxy[iMu]) < 0.02) && (std::abs(muonEta[iMu]) < 3.0);
	if (veto) {return veto;} 
    }
    return veto;
}

bool GoodElVeto(int nElectron, RVec<float> elPt, RVec<bool> elIso, RVec<float> dxy, RVec<float> elEta) {
    if (nElectron < 1) {return false;}
    bool veto =false;
    for (int iEl = 0; iEl < elPt.size(); iEl++) {
// DP EDIT
//	veto = (elPt[iEl] > 35.) && (elIso[iEl] == 1) && (std::abs(dxy[iEl]) < 0.05) && (std::abs(elEta[iEl]) < 2.5);
      veto = (elPt[iEl] > 15.0) && (elIso[iEl] == 1) && (std::abs(dxy[iEl]) < 0.05) && (std::abs(elEta[iEl]) < 3.0);
	if (veto) {return veto;}
    }
    return veto;
}
// -------------------------------------------------------------------------------

double RAND() {
    const double from = 0.0;
    const double to = 1.0;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(from, to);
    return distr(generator);
}

float getSF(int jetCat, float pt, std::string _year, int _var) {
    if ((jetCat!=0) && (jetCat!=1)) {std::cerr<<"ERROR - INITIAL JET CATEGORY MUST BE 0 or 1. Current value is " << jetCat << std::endl;}
    // SF[_var][pt] (_var: 0=nom, 1=up, 2=down)
    // variations are described above, pt cats are [300, 400), [400, 480), [480, 600), [600, 1200) across all years
    // Scale factors located at:
    // https://indico.cern.ch/event/1152827/contributions/4840404/attachments/2428856/4162159/ParticleNet_SFs_ULNanoV9_JMAR_25April2022_PK.pdf
    // HP (tight WP)
    float SF2016APV_T[3][4] = {{1.10,1.06,1.04,1.00},{1.18,1.13,1.11,1.21},{1.03,1.00,0.98,0.91}};
    float SF2016_T[3][4]    = {{0.97,0.91,0.99,1.00},{1.07,0.96,1.05,1.09},{0.89,0.86,0.94,0.92}};
    float SF2017_T[3][4]    = {{1.12,0.96,1.00,0.93},{1.24,1.01,1.05,0.98},{1.02,0.92,0.95,0.87}};
    float SF2018_T[3][4]    = {{1.03,0.95,0.91,0.95},{1.12,1.00,0.95,1.02},{0.95,0.90,0.88,0.90}};
    // MP (loose WP)
    float SF2016APV_L[3][4] = {{1.23,1.07,1.04,1.06},{1.39,1.17,1.18,1.24},{1.09,1.02,0.99,0.96}};
    float SF2016_L[3][4]    = {{1.08,0.99,1.03,1.29},{1.19,1.05,1.10,1.54},{0.98,1.04,0.98,1.03}};
    float SF2017_L[3][4]    = {{1.11,1.01,1.05,1.00},{1.23,1.05,1.14,1.06},{1.03,0.97,1.01,0.96}};
    float SF2018_L[3][4]    = {{1.19,0.98,0.96,0.97},{1.31,1.02,1.00,1.02},{1.07,0.94,0.93,0.92}};
    // begin logic
    float SF;
    float SF1;
    float SF2;
    int ptCat;
    // get the pT category
// DP EDIT - the commented out line...
//    if ((pt >= 300) && (pt < 400))       { ptCat = 0; }
    if ((pt >= 350) && (pt < 400))       { ptCat = 0; }
    else if ((pt >= 400) && (pt < 480))  { ptCat = 1; }
    else if ((pt >= 480) && (pt < 600))  { ptCat = 2; }
    else if ((pt >= 600) && (pt < 1200)) { ptCat = 3; }
    else { return 1.0; }
    // get SF
    switch (jetCat) {
        case 0: { // jet is originally in fail
            if (_year=="2016APV") { SF = SF2016APV_L[_var][ptCat]; }
            else if (_year=="2016") { SF = SF2016_L[_var][ptCat]; }
            else if (_year=="2017") { SF = SF2017_L[_var][ptCat]; }
            else { SF = SF2018_L[_var][ptCat]; }
	    break;
	}
        case 1: { // jet originally in pass
            if (_year=="2016APV") { SF = SF2016APV_T[_var][ptCat]; }
            else if (_year=="2016") { SF = SF2016_T[_var][ptCat]; }
            else if (_year=="2017") { SF = SF2017_T[_var][ptCat]; }
            else { SF = SF2018_T[_var][ptCat]; }
	    break;
	}
    }
    if (SF<0.7) {std::cerr << "SF is " << SF << " for variation " << _var << " and ptCat " << ptCat << std::endl;}
    return SF;
}

int getDiPhotonCat(RVec<float> taggerVal, float taggerWP) {
    int wpCat;
    if ((taggerVal[0] >= taggerWP) && (taggerVal[1] >= taggerWP)) {
       wpCat = 2;
    } else if ((taggerVal[0] < taggerWP) && (taggerVal[1] < taggerWP)) {
       wpCat = 0;
    } else {
       wpCat = 1;
    }
//std::cout << "DiPhoton Cat " << taggerVal[0] << " " << taggerVal[1] << " " << wpCat << std::endl;
    return wpCat;
}

float getPhotonSF(RVec<float> pt, RVec<float> eta, RVec<float> taggerVal, float taggerWP, int ichoice) {
    int wpCat;
    int wpCat1 = int(taggerVal[0]);
    int wpCat2 = int(taggerVal[1]);
    int ptCat1=0;
    int ptCat2=0;
    int etaCat1=0;
    int etaCat2=0;
    float SF;
    float SF1;
    float SF2;
//    float SF2016APV_T[10][5] = {{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0}};
    float SF2016_T[10][5]    = {{1.045918345451355,1.0369843244552612,1.0328317880630493,1.022757649421692,1.0284552574157715},{0.9825396537780762,0.9904502034187317,0.9920424222946167,1.0052909851074219,1.0647886991500854},{1.041338562965393,1.015455961227417,1.0125983953475952,1.0856643915176392,1.0301003456115723},{1.0032051801681519,0.9972106218338013,0.9917582273483276,1.0233516693115234,1.0361111164093018},{0.9817578792572021,0.9759206771850586,0.9681440591812134,1.0055943727493286,1.0070126056671143},{0.9884488582611084,0.9816384315490723,0.969613254070282,0.9986013770103455,1.0166898965835571},{1.0097087621688843,0.9985975027084351,0.9944751262664795,1.025174856185913,1.026063084602356},{0.9778671860694885,1.007836937904358,0.9842767119407654,0.975944995880127,1.0263158082962036},{0.9760383367538452,0.9904109835624695,0.9826666712760925,1.0405954122543335,1.0272108316421509},{1.0439189672470093,1.0339462757110596,1.028493881225586,1.061662197113037,1.058432936668396}};
//    float SF2017_T[10][5]    = {{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0}};
//    float SF2018_T[10][5]    = {{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0}};
//        float SF2016APV_M[10][5] = {{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0}};
    float SF2016_M[10][5]    = {{1.0384024381637573,1.0325098037719727,1.0225000381469727,1.0295202732086182,1.0183823108673096},{0.9677419066429138,0.9853300452232361,0.9892473220825195,0.9869513511657715,1.028678297996521},{1.052901029586792,1.0067843198776245,1.0013889074325562,0.9956011772155762,1.0119940042495728},{1.0,0.9914634227752686,0.9890776872634888,1.0245699882507324,1.0260869264602661},{0.9548022747039795,0.9740419983863831,0.9670329689979553,1.0037267208099365,1.0162907838821411},{0.9831223487854004,0.9802225232124329,0.9694749712944031,0.9975185990333557,1.009987473487854},{0.9986187815666199,0.9926560521125793,0.9878345727920532,1.018587350845337,1.0270936489105225},{1.0260416269302368,1.0027472972869873,0.9707520604133606,1.0252366065979004,1.0160771608352661},{0.9731638431549072,0.9852941036224365,0.9844311475753784,1.0365407466888428,1.0182703733444214},{1.0428789854049683,1.0284605026245117,1.0124223232269287,1.0602705478668213,1.0550122261047363}};
//    float SF2017_M[10][5]    = {{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0}};
//    float SF2018_M[10][5]    = {{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0}};
//    float SF2016APV_L[10][5] = {{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0}};
    float SF2016_L[10][5]    = {{1.0037453174591064,1.0,0.9967776536941528,1.006355881690979,0.9851064085960388},{0.9605911374092102,0.969264566898346,0.9806034564971924,0.9892473220825195,1.0344061851501465},{1.035433053970337,0.9954954981803894, 0.9909604787826538, 1.0199296474456787,0.9351415038108826,},{0.9953106641769409,0.9913138151168823,0.9825897812843323,1.0178173780441284,0.9988725781440735},{0.989195704460144,0.9824753403663635,0.9769737124443054,1.0066889524459839,0.994356632232666},{0.9856114983558655,0.9835526347160339,0.9769483804702759,1.0033669471740723,1.0158730745315552},{0.9952940940856934,0.990217387676239,0.984749436378479,1.0066889524459839,1.0135135650634766},{0.9880319237709045,0.9966063499450684,0.9897843599319458,1.0023781061172485,1.0236612558364868},{0.9628252983093262,0.9790979027748108,0.982758641242981,1.0129870176315308,0.9967355728149414},{0.99622642993927,0.9934065937995911,0.9967880249023438,1.022459864616394,1.0201913118362427}};
//    float SF2017_L[10][5] = {{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0}};
//    float SF2018_L[10][5] = {{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0,1.0}};
    //
    if ((taggerVal[0] >= taggerWP) && (taggerVal[1] >= taggerWP)) {
       wpCat = 2;
    } else if ((taggerVal[0] < taggerWP) && (taggerVal[1] < taggerWP)) {
       wpCat = 0;
    } else {
       wpCat = 1;
    }
      if (eta[0] >= 2.0) {
    etaCat1 = 9;
  } else if ((eta[0] >= 1.566) && (eta[0] < 2.0)) {
    etaCat1 = 8;
  } else if ((eta[0] >= 1.444) && (eta[0] < 1.566)) {
    etaCat1 = 7;
  } else if ((eta[0] >= 0.8) && (eta[0] < 1.444)) {
    etaCat1 = 6;
  } else if ((eta[0] >= 0.0) && (eta[0] < 0.8)) {
    etaCat1 = 5;
  } else if ((eta[0] >= -0.8) && (eta[0] < 0.0)) {
    etaCat1 = 4;
  } else if ((eta[0] >= -1.444) && (eta[0] < -0.8)) {
    etaCat1 = 3;
  } else if ((eta[0] >= -1.566) && (eta[0] < -1.444)) {
    etaCat1 = 2;
  } else if ((eta[0] >= -2.0) && (eta[0] < -1.566)) {
    etaCat1 = 1;
  } else {
    etaCat1 = 0;
  }
  if (eta[1] >= 2.0) {
    etaCat2 = 9;
  } else if ((eta[1] >= 1.566) && (eta[1] < 2.0)) {
    etaCat2 = 8;
  } else if ((eta[1] >= 1.444) && (eta[1] < 1.566)) {
    etaCat2 = 7;
  } else if ((eta[1] >= 0.8) && (eta[1] < 1.444)) {
    etaCat2 = 6;
  } else if ((eta[1] >= 0.0) && (eta[1] < 0.8)) {
    etaCat2 = 5;
  } else if ((eta[1] >= -0.8) && (eta[1] < 0.0)) {
    etaCat2 = 4;
  } else if ((eta[1] >= -1.444) && (eta[1] < -0.8)) {
    etaCat2 = 3;
  } else if ((eta[1] >= -1.566) && (eta[1] < -1.444)) {
    etaCat2 = 2;
  } else if ((eta[1] >= -2.0) && (eta[1] < -1.566)) {
    etaCat2 = 1;
  } else {
    etaCat2 = 0;
  }
  if (pt[0] >= 120.0) {
    ptCat1 = 4;
  } else if ((pt[0] >= 80.0) && (pt[0] < 120.0)) {
    ptCat1 = 3;
  } else if ((pt[0] >= 50.0) && (pt[0] < 80.0)) {
    ptCat1 = 2;
  } else if ((pt[0] >= 35.0) && (pt[0] < 50.0)) {
    ptCat1 = 1;
  } else {
    ptCat1 = 0;
  }
  if (pt[1] >= 120.0) {
    ptCat2 = 4;
  } else if ((pt[1] >= 80.0) && (pt[1] < 120.0)) {
    ptCat2 = 3;
  } else if ((pt[1] >= 50.0) && (pt[1] < 80.0)) {
    ptCat2 = 2;
  } else if ((pt[1] >= 35.0) && (pt[1] < 50.0)) {
    ptCat2 = 1;
  } else {
    ptCat2 = 0;
  }


  switch (wpCat1) {
    case 0: 
      SF1 = 1.0;
    case 1: 
      SF1 = SF2016_L[etaCat1][ptCat1];
//      if (year=="2016APV") { SF1 = SF2016APV_L[etaCat1][ptCat1]; }
//      else if (year=="2016") { SF1 = SF2016_L[etaCat1][ptCat1]; }
//      else if (year=="2017") { SF1 = SF2017_L[etaCat1][ptCat1]; }
//      else { SF1 = SF2018_L[etaCat1][ptCat1]; }
    case 2:
      SF1 = SF2016_M[etaCat1][ptCat1];
//      if (year=="2016APV") { SF1 = SF2016APV_M[etaCat1][ptCat1]; }
//      else if (year=="2016") { SF1 = SF2016_M[etaCat1][ptCat1]; }
//      else if (year=="2017") { SF1 = SF2017_M[etaCat1][ptCat1]; }
//      else { SF1 = SF2018_M[etaCat1][ptCat1]; }
    case 3:
      SF1 = SF2016_T[etaCat1][ptCat1];
//      if (year=="2016APV") { SF1 = SF2016APV_T[etaCat1][ptCat1]; }
//      else if (year=="2016") { SF1 = SF2016_T[etaCat1][ptCat1]; }
//      else if (year=="2017") { SF1 = SF2017_T[etaCat1][ptCat1]; }
//      else { SF1 = SF2018_T[etaCat1][ptCat1]; }
  }
  switch (wpCat2) {
    case 0:
      SF2 = 1.0;
    case 1: 
      SF2 = SF2016_L[etaCat2][ptCat2];
//      if (year=="2016APV") { SF2 = SF2016APV_L[etaCat2][ptCat2]; }
//      else if (year=="2016") { SF2 = SF2016_L[etaCat2][ptCat2]; }
//      else if (year=="2017") { SF2 = SF2017_L[etaCat2][ptCat2]; }
//      else { SF2 = SF2018_L[etaCat2][ptCat2]; }
    case 2:
        SF2 = SF2016_M[etaCat2][ptCat2];
//      if (year=="2016APV") { SF2 = SF2016APV_M[etaCat2][ptCat2]; }
//      else if (year=="2016") { SF2 = SF2016_M[etaCat2][ptCat2]; }
//      else if (year=="2017") { SF2 = SF2017_M[etaCat2][ptCat2]; }
//      else { SF2 = SF2018_M[etaCat2][ptCat2]; }
    case 3:
        SF2 = SF2016_T[etaCat2][ptCat2];
//      if (year=="2016APV") { SF2 = SF2016APV_T[etaCat2][ptCat2]; }
//      else if (year=="2016") { SF2 = SF2016_T[etaCat2][ptCat2]; }
//      else if (year=="2017") { SF2 = SF2017_T[etaCat2][ptCat2]; }
//      else { SF2 = SF2018_T[etaCat2][ptCat2]; }
  }
  SF = SF1 * SF2;
//std::cout << "SF1 SF2 SF " << SF1 << " " <<SF2 << " " << SF << std::endl;
  if (ichoice == 1) SF = SF1;
  if (ichoice == 2) SF = SF2;
  return SF;
};

int getOriginalTopCat(float TvsQCD, float WP) {
    if (TvsQCD > WP) {return 1;}    // pass
    else {return 0;}                // fail
}

int getNewTopCat(float SF, int oldCat, float eff, double rand, bool invert) {
    int newCat = oldCat;    // to be (possibly) updated
    if (SF < 1) {
        // downgrade fraction (1-SF) of tagged -> untagged
        if ((oldCat == 1) && (rand < 1.-SF)) { 
	    newCat=0; 
	}
    }
    else {
        // upgrade fraction of untagged -> tagged
        if (oldCat == 0) {
            float num = 1.-SF;
            float den = 1.-(1./eff);
            float f = num/den;
            if (rand < f) {	
                newCat = 1;
            }
        }
    }
    //delete random;
    return newCat;
};

int updatePhotonTag(int photonCat, RVec<float> pt, RVec<float> eta, RVec<int> taggerVal, float taggerWP, float effP, float effF) {
    /* updates the tagger category for phi jets
 *      * https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#2a_Jet_by_jet_updating_of_the_b
 *           * Params:
 *                *   photonCat    = original photon category (0: failfail, 1: fail, 2: pass)
 *                     *   pt        = jet pt
 *                          *   taggerVal = particleNet tagger value
 *                              */
    float SF;// = getPhotonSF(pt, eta, taggerVal, taggerWP, 0);
    float SF1;// = getPhotonSF(pt, eta, taggerVal, taggerWP, 1);
    float SF2;// = getPhotonSF(pt, eta, taggerVal, taggerWP, 2);
    float SF_P;
    float SF_F;
    float eff_P = effP;
    float eff_F = effF;
    double rn = RAND();
    int newCat = photonCat;
    std::vector<float> TagVP(2);
    TagVP[1] = 1.0;
    TagVP[2]=1.0;
    SF = getPhotonSF(pt, eta, taggerVal, taggerWP, 0);
    SF1 = getPhotonSF(pt, eta, taggerVal, taggerWP, 1);
    SF2 = getPhotonSF(pt, eta, taggerVal, taggerWP, 2);
    SF_P = getPhotonSF(pt, eta, TagVP, taggerWP, 0);
    TagVP[2]=0.0;
    SF_F = getPhotonSF(pt, eta, TagVP, taggerWP, 0);    
    if ((SF_P < 1) && (SF_F < 1)) {
        if ( (newCat==2) && (rn < (1.-SF_P)) ) newCat=0;        // tight (2) -> untag (0)
        if ( (newCat==1) && (rn < (1.-SF_F)) ) newCat=0;        // loose (1) -> untag (0)
    }
    if ((SF_P > 1) && (SF_F > 1)) {
        float fF, fP;
        if (newCat==0) {
            float num = eff_P*(SF_P-1.);
            float den = 1.-(eff_F+eff_P);
            fP = num/den;
            if (rn < fP) newCat=2;      // untag (0) -> tight (2)
            else {
                rn = RAND();
                num = eff_F*(SF_F-1.);
                den = (1.-eff_F-eff_P)*(1.-fP);
                fF = num/den;
                if (rn < fF) newCat=1;  // loose (1) -> tight (2)
            }
        }
    }
    if ((SF_F < 1) && (SF_P > 1)) {
        if (newCat==0) {
            float num = eff_P*(SF_P-1.);
            float den = 1.-(eff_F+eff_P);
            float f = num/den;
            if (rn < f) newCat=2;           // untag (0) -> tight (2)
        }
        if (newCat==1) {
            if (rn < (1.-SF_F)) newCat=0;   // loose (1) -> untag (0)
        }
    }
    if ((SF_F > 1) && (SF_P < 1)) {
        if ((newCat==2) && (rn < (1.-SF_P))) newCat=0;  // tight (2) -> untag (0)
        if (newCat==0) {
            float num = eff_F*(SF_F-1.);
            float den = 1.-(eff_F+eff_P);
            float f = num/den;
            if (rn < f) newCat=1;       // untag (0) -> loose (1)
        }
    }
    return newCat;
}

RVec<int> PickTopWithSFs2(RVec<float> TvsQCD,
                  RVec<float> pt,
                  RVec<int> idxs,
                  float scoreCut,
                  float eff0,
                  float eff1,
                  std::string year,
                  int variation) {
    if (idxs.size() > 2) {
        std::cout << "PickTop - WARNING: You have input more than two indices. Only two jet indices accepted, assuming first two indices.\n";
    }

    std::vector<int> out(2);
    float WP = scoreCut;
    int idx0 = idxs[0];
    int idx1 = idxs[1];
    bool invertScore=false;
    bool isTop0=false;
    bool isTop1=false;
    int orig_score0 = getOriginalTopCat(TvsQCD[idx0], WP);
    int orig_score1 = getOriginalTopCat(TvsQCD[idx1], WP);
    float SF0 = getSF(orig_score0, pt[idx0], year, variation);
    float SF1 = getSF(orig_score1, pt[idx1], year, variation);
    if ((SF0==0) || (SF1==0)) {std::cerr<<"SF is 0\n";}
    double rand0 = RAND();
    double rand1 = RAND();
    int new_score0 = getNewTopCat(SF0, orig_score0, eff0, rand0, invertScore);
    int new_score1 = getNewTopCat(SF1, orig_score1, eff1, rand1, invertScore);
    isTop0 = (new_score0 == 1);
    isTop1 = (new_score1 == 1);
    if (invertScore) {
    if (isTop0 && isTop1) {
        out[0] = idx0;
        out[1] = idx1;
    } else if (isTop0) {
        out[0] = idx0;
        out[1] = idx1;
    } else if (isTop1) {
        out[0] = idx1;
        out[1] = idx0;
    } else {
        out[0] = -1;
        out[1] = -1;
    }}
    return out;
}

RVec<int> PickTopWithSFs(RVec<float> TvsQCD, 
                  RVec<float> HbbvsQCD, 
                  RVec<float> pt, 
                  RVec<int> idxs, 
                  float scoreCut, 
                  float eff0,
                  float eff1,
                  std::string year,
                  int variation,
                  bool invertScore=false) {
    if (idxs.size() > 2) {
        std::cout << "PickTop - WARNING: You have input more than two indices. Only two jet indices accepted, assuming first two indices.\n";
    }

    std::vector<int> out(2);
    float WP = scoreCut;
    int idx0 = idxs[0];
    int idx1 = idxs[1];

    // first, need to determine original tagger categories for each jet 
    int orig_score0 = getOriginalTopCat(TvsQCD[idx0], WP);
    int orig_score1 = getOriginalTopCat(TvsQCD[idx1], WP);
    // now, determine new tagger category
    float SF0 = getSF(orig_score0, pt[idx0], year, variation);
    float SF1 = getSF(orig_score1, pt[idx1], year, variation);
    if ((SF0==0) || (SF1==0)) {std::cerr<<"SF is 0\n";}
    double rand0 = RAND();
    double rand1 = RAND();
    int new_score0 = getNewTopCat(SF0, orig_score0, eff0, rand0, invertScore);
    int new_score1 = getNewTopCat(SF1, orig_score1, eff1, rand1, invertScore);

    // now perform top selection using the new tagger "scores" (aka categories)
    bool isTop0, isTop1;
    if (!invertScore) {
        isTop0 = (new_score0 == 1);
        isTop1 = (new_score1 == 1);
    }
    else {// higgs veto, use raw HbbvsQCD score.
        isTop0 = (new_score0 == 0) && (HbbvsQCD[idx0] < 0.2);
        isTop1 = (new_score1 == 0) && (HbbvsQCD[idx1] < 0.2);
    }
    // determine which is which 
    if (isTop0 && isTop1) {
	// DEBUG
	//std::cout << "PICKTOP (SIGNAL) - event has 2 top jets\n";
        // if both pass as Top, use the raw TvsQCD score to determine which is "real" top
        if (TvsQCD[idx0] > TvsQCD[idx1]) {
            out[0] = idx0;
            out[1] = idx1;
        }
        else {
            out[0] = idx1;
            out[1] = idx0;
        }
    }
    else if (isTop0) {
	// DEBUG
	//std::cout << "PICKTOP (SIGNAL) - event has 1 top jet\n";
        out[0] = idx0;
        out[1] = idx1;
    }
    else if (isTop1) {
	// DEBUG
	//std::cout << "PICKTOP (SIGNAL) - event has 1 top jet\n";
        out[0] = idx1;
        out[1] = idx0;
    }
    else {
	// DEBUG
	//std::cout << "PICKTOP (SIGNAL) - event has NO top jets\n";
        out[0] = -1;
        out[1] = -1;
    }
    return out;
}

RVec<int> PickDijets(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    int jet0Idx = -1;
    int jet1Idx = -1;
    for (int ijet = 0; ijet < pt.size(); ijet++) {
        if (jet1Idx == -1) {
// DP EDIT
//            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50) {
            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4) {
                if (jet0Idx == -1) {
                    jet0Idx = ijet;
                } else {
//DP EDIT
//                    if (abs(hardware::DeltaPhi(phi[jet0Idx], phi[ijet])) > M_PI/2) {
                    if (std::abs(eta[jet0Idx]-eta[ijet]) < 1.6) {
                        jet1Idx = ijet;
                        break;
                    }
                }
            }
        }       
    }
    return {jet0Idx,jet1Idx};
}

RVec<int> PickDiphotons(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    int ph0Idx = -1;
    int ph1Idx = -1;
    for (int iph = 0; iph < pt.size(); iph++) {
        if (ph1Idx == -1) {
            if (pt[iph] > 25) {
                if (ph0Idx == -1) {
                    ph0Idx = iph;
                } else {
                    ph1Idx = iph;
                    break;
                }
            }
        }
    }
    return {ph0Idx,ph1Idx};
}

RVec<float> PickDijetsV(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, RVec<float> Jet_btagCMVA) {
    int jet0Idx = -1;
    int jet1Idx = -1;
    float pt0 = -1.0;
    float pt1 = -1.0;
    float eta0 = -1.0;
    float eta1 = -1.0;
    float phi0 = -1.0;
    float phi1 = -1.0;
    float mass0 = -1.0;
    float mass1 = -1.0;
    float tag0 = -1.0;
    float tag1 = -1.0;
    float IdReg = -1.0; //DP: 0 for CR, 1 for SR
    for (int ijet = 0; ijet < pt.size(); ijet++) {
        if (jet1Idx == -1) {
// DP EDIT
// //            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50) {
            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4) {
                if (jet0Idx == -1) {
                    jet0Idx = ijet;
                } else {
//DP EDIT
////                    if (abs(hardware::DeltaPhi(phi[jet0Idx], phi[ijet])) > M_PI/2) {
                        jet1Idx = ijet;
                        break;
//                }
                }
            }
        }
    }
    if (jet0Idx>-1){
        pt0 = pt[jet0Idx];
        eta0 = eta[jet0Idx];
        phi0 = phi[jet0Idx];
        mass0 = mass[jet0Idx];
        tag0 = Jet_btagCMVA[jet0Idx];
    }
    if (jet1Idx>-1){
        pt1 = pt[jet1Idx];
        eta1 = eta[jet1Idx];
        phi1 = phi[jet1Idx];
        mass1 = mass[jet1Idx];
        tag1 = Jet_btagCMVA[jet1Idx];
    }
    if ((tag0>0.8) && (tag1>0.8)) {
       IdReg = 1.0;
    } else {
       if ((tag0>0.8) || (tag1>0.8)) IdReg = 0.0;
    } 
    return {IdReg,pt0,pt1,eta0,eta1,phi0,phi1,mass0,mass1,tag0,tag1};
}

RVec<int> PickDijetsV3(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, RVec<float> Jet_btagCMVA) {
    int jet0Idx = -1;
    int jet1Idx = -1;
    for (int ijet = 0; ijet < pt.size(); ijet++) {
        if (jet1Idx == -1) {
// DP EDIT
//            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50) {
            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && Jet_btagCMVA[ijet] > 0.8) {
                if (jet0Idx == -1) {
                    jet0Idx = ijet;
                } else {
//DP EDIT
////                    if (abs(hardware::DeltaPhi(phi[jet0Idx], phi[ijet])) > M_PI/2) {
                    if (std::abs(eta[jet0Idx]-eta[ijet]) < 1.6) {
                        jet1Idx = ijet;
                        break;
                    }
                }
            }
        }
    }
    return {jet0Idx,jet1Idx};
}

RVec<float> PickDijetsV3_ALL(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, RVec<float> Jet_btagCMVA) {
    int jet0Idx = -1;
    int jet1Idx = -1;
    float pt0 = -1.0;
    float pt1 = -1.0;
    float eta0 = -1.0;
    float eta1 = -1.0;
    float phi0 = -1.0;
    float phi1 = -1.0;
    float mass0 = -1.0;
    float mass1 = -1.0;
    float tag0 = -1.0;
    float tag1 = -1.0;
    float IdReg = -1.0; //DP: 0 for CR, 1 for SR
    for (int ijet = 0; ijet < pt.size(); ijet++) {
        if (jet1Idx == -1) {
// DP EDIT FIRST WORK ON SR:
//            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50) {
            if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && Jet_btagCMVA[ijet] > 0.8) {
                if (jet0Idx == -1) {
                    jet0Idx = ijet;
                } else {
//DP EDIT
//////                    if (abs(hardware::DeltaPhi(phi[jet0Idx], phi[ijet])) > M_PI/2) {
                    if (std::abs(eta[jet0Idx]-eta[ijet]) < 1.6 && Jet_btagCMVA[ijet] > 0.8) {
                        jet1Idx = ijet;
                        IdReg = 1.0;
                        break;
                    }
                }
            }
        }
    }
// DP EDIT NOW WORK ON CR:
    if (jet1Idx == -1 && jet0Idx > -1) {
         for (int ijet = 0; ijet < pt.size(); ijet++) {
            if (ijet != jet0Idx && jet1Idx == -1 && pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && Jet_btagCMVA[ijet] > 0.2 && Jet_btagCMVA[ijet] < 0.8 && std::abs(eta[jet0Idx]-eta[ijet]) < 1.6) {
                 jet1Idx = ijet;
                 IdReg = 0.0;
                 break;
            }
        }
    }
//DP EDIT: put everything out!
//    if (IdReg > -1) {
        pt0 = pt[jet0Idx];
        eta0 = eta[jet0Idx];
        phi0 = phi[jet0Idx];
        mass0 = mass[jet0Idx];
        tag0 = Jet_btagCMVA[jet0Idx];
        pt1 = pt[jet1Idx];
        eta1 = eta[jet1Idx];
        phi1 = phi[jet1Idx];
        mass1 = mass[jet1Idx];
        tag1 = Jet_btagCMVA[jet1Idx];
///    }
    return {IdReg, pt0,pt1,eta0,eta1,phi0,phi1,mass0,mass1,tag0,tag1};
}

RVec<int> PickDiphotons_X(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, float SmassTarget) {
    int ph0Idx = -1;
    int ph1Idx = -1;
    float SmassTry = 0.0;
    float Smass = -1.0;
    float dSmass = 1000000.0;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);
    for (int iph1 = 0; iph1 < pt.size(); iph1++) {
     Lvector1 = Lvector[iph1];
     for (int iph2 = iph1+1; iph2 < pt.size(); iph2++) {
       Lsum.SetCoordinates(0,0,0,0);
       Lsum = Lvector1+Lvector[iph2];
       SmassTry = Lsum.M();
       if (abs(SmassTry-SmassTarget) < dSmass) {
          dSmass = abs(SmassTry-SmassTarget);
          ph0Idx = iph1;
          ph1Idx = iph2;
          Smass = SmassTry;
       }
     }
    }
    return {ph0Idx,ph1Idx};
}

RVec<int> PickDiphotonsLeading(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    int ph1Idx = -1;
    int ph2Idx = -1;
    float Smass = -1.0;
    float dR = -1.0;
    float pt1 = -1.0;
    float pt2 = -1.0;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);
    for (int iph1 = 0; iph1 < pt.size(); iph1++) {
     if (pt[iph1]>pt1) {
       pt1 = pt[iph1];
       ph1Idx = iph1;
     }
    }
    for (int iph2 = 0; (iph2 < pt.size()); iph2++) {
     if ((iph2!=ph1Idx)&&(pt[iph2]>pt2)) {
       pt2 = pt[iph2];
       ph2Idx = iph2;
     }
    }
    if ((ph1Idx>-1)&&(ph2Idx>-1)) {
     Lvector1 = Lvector[ph1Idx];
     Lsum.SetCoordinates(0,0,0,0);
     Lvector2 = Lvector[ph2Idx];
     Lsum = Lvector1+Lvector2;
     Smass = Lsum.M();
     dR = hardware::DeltaR(Lvector1,Lvector2);
    }
    return {ph1Idx,ph2Idx};
}

RVec<int> PickDiphotonsLeadingOrdered(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    int ph1Idx = -1;
    int ph2Idx = -1;
    float Smass = -1.0;
    float dR = -1.0;
    float pt1 = -1.0;
    float pt2 = -1.0;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);
    if (pt.size()>1) {
       pt1 = pt[0];
       ph1Idx = 0;
       pt2 = pt[1];
       ph2Idx = 1;
    }
    if ((ph1Idx>-1)&&(ph2Idx>-1)) {
     Lvector1 = Lvector[ph1Idx];
     Lsum.SetCoordinates(0,0,0,0);
     Lvector2 = Lvector[ph2Idx];
     Lsum = Lvector1+Lvector2;
     Smass = Lsum.M();
     dR = hardware::DeltaR(Lvector1,Lvector2);
    }
    return {ph1Idx,ph2Idx};
}

RVec<int> FindMothersPdgId(RVec<int> id, RVec<int> idM) {
    RVec<int> motherId;
    int icurrM;
    for (int i = 0; i < idM.size(); i++) {
     icurrM = idM[i];
     if (icurrM>=0) {
      motherId[i] = id[icurrM];
     } else {
      motherId[i] = -1;
     }
    }
    return {motherId};
}

RVec<int> PickGENDiphotonsLeading(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, RVec<int> id, int IdUse, RVec<int> idM, int IdMUse) {
    int ph1Idx = -1;
    int ph2Idx = -1;
    float Smass = -1.0;
    float dR = -1.0;
    float pt1 = -1.0;
    float pt2 = -1.0;
    RVec<int> motherId=idM;
    int icurrM;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);

    for (int i = 0; i < idM.size(); i++) {
     icurrM = idM[i];
     if (icurrM>=0) {
      motherId[i] = id[icurrM];
     } else {
      motherId[i] = -1;
     }
    }

    for (int iph1 = 0; iph1 < pt.size(); iph1++) {
     if ((id[iph1]==IdUse)&&(motherId[iph1]==IdMUse)&&(pt[iph1]>pt1)) {
//     if ((id[iph1]==IdUse)&&(pt[iph1]>pt1)) {
       pt1 = pt[iph1];
       ph1Idx = iph1;
     }
    }
    for (int iph2 = 0; (iph2 < pt.size()); iph2++) {
     if ((iph2!=ph1Idx)&&(id[iph2]==IdUse)&&(idM[iph2]==idM[ph1Idx])) {
//     if ((id[iph2]==IdUse)&&(motherId[iph2]==IdMUse)&&(pt[iph2]>pt2)&&(idM[iph2]==idM[ph1Idx])) {
//     if ((id[iph2]==IdUse)&&(pt[iph2]>pt2)) {
       pt2 = pt[iph2];
       ph2Idx = iph2;
     }
    }
    if ((ph1Idx>-1)&&(ph2Idx>-1)) {
     Lvector1 = Lvector[ph1Idx];
     Lsum.SetCoordinates(0,0,0,0);
     Lvector2 = Lvector[ph2Idx];
     Lsum = Lvector1+Lvector2;
     Smass = Lsum.M();
     dR = hardware::DeltaR(Lvector1,Lvector2);
    }
    return {ph1Idx,ph2Idx};
}

float SmassCalc(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, float SmassTarget) {
    int ph0Idx = -1;
    int ph1Idx = -1;
    float SmassTry = 0.0;
    float Smass = -1.0;
    float dSmass = 1000000.0;
    float dR = -1.0;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);
    for (int iph1 = 0; iph1 < pt.size(); iph1++) {
     Lvector1 = Lvector[iph1];
     for (int iph2 = iph1+1; iph2 < pt.size(); iph2++) {
       Lsum.SetCoordinates(0,0,0,0);
       Lvector2 = Lvector[iph2];
       Lsum = Lvector1+Lvector2;
       SmassTry = Lsum.M();
       if (abs(SmassTry-SmassTarget) < dSmass) {
          dSmass = abs(SmassTry-SmassTarget);
          ph0Idx = iph1;
          ph1Idx = iph2;
          Smass = SmassTry;
          dR = hardware::DeltaR(Lvector1,Lvector2);
       }
     }
    }
    return Smass;
}

float TPmassCalcLeading(RVec<float> AA_vector, RVec<int> FJids, RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1=hardware::TLvector(AA_vector[0],AA_vector[2],AA_vector[4],AA_vector[6]);
    ROOT::Math::PtEtaPhiMVector Lvector2=hardware::TLvector(AA_vector[1],AA_vector[3],AA_vector[5],AA_vector[7]);
    ROOT::Math::PtEtaPhiMVector Lvector3;
    float TPmass=-1.0;
    if ((FJids[0]>-1)&&(AA_vector[0]>0.0)&&(AA_vector[1]>0.0)) {
       Lvector3=hardware::TLvector(pt[FJids[0]],eta[FJids[0]],phi[FJids[0]],mass[FJids[0]]);
       Lsum.SetCoordinates(0,0,0,0);
       Lsum = Lvector1+Lvector2+Lvector3;
       TPmass = Lsum.M();
    }
    return TPmass;
}

float TPmassCalcLeading80(RVec<float> AA_vector, RVec<int> FJids, RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, RVec<float> tag) {
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1=hardware::TLvector(AA_vector[0],AA_vector[2],AA_vector[4],AA_vector[6]);
    ROOT::Math::PtEtaPhiMVector Lvector2=hardware::TLvector(AA_vector[1],AA_vector[3],AA_vector[5],AA_vector[7]);
    ROOT::Math::PtEtaPhiMVector Lvector3;
    float TPmass=-1.0;
    int itop=-1;
    if ((FJids[0]>-1)&&(tag[FJids[0]]>0.8)) {
       itop=0;
    } else {
       if ((FJids[1]>-1)&&(tag[FJids[1]]>0.8)) {
         itop=1;
       }
    }
    if ((itop>-1)&&(AA_vector[0]>0.0)&&(AA_vector[1]>0.0)) {
       Lvector3=hardware::TLvector(pt[FJids[itop]],eta[FJids[itop]],phi[FJids[itop]],mass[FJids[itop]]);
       Lsum.SetCoordinates(0,0,0,0);
       Lsum = Lvector1+Lvector2+Lvector3;
       TPmass = Lsum.M();
    }
    return TPmass;
}

float TPmassCalcLeading70(RVec<float> AA_vector, RVec<int> FJids, RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, RVec<float> tag) {
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1=hardware::TLvector(AA_vector[0],AA_vector[2],AA_vector[4],AA_vector[6]);
    ROOT::Math::PtEtaPhiMVector Lvector2=hardware::TLvector(AA_vector[1],AA_vector[3],AA_vector[5],AA_vector[7]); 
    ROOT::Math::PtEtaPhiMVector Lvector3;
    float TPmass=-1.0;
    int itop=-1;
    if ((FJids[0]>-1)&&(tag[FJids[0]]>0.7)) {
       itop=0;
    } else {
       if ((FJids[1]>-1)&&(tag[FJids[1]]>0.7)) {
         itop=1;
       }
    }
    if ((itop>-1)&&(AA_vector[0]>0.0)&&(AA_vector[1]>0.0)) {
       Lvector3=hardware::TLvector(pt[FJids[itop]],eta[FJids[itop]],phi[FJids[itop]],mass[FJids[itop]]);
       Lsum.SetCoordinates(0,0,0,0);
       Lsum = Lvector1+Lvector2+Lvector3;
       TPmass = Lsum.M();
    }
    return TPmass;
}

float SmassCalcLeading(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    int ph1Idx = -1;
    int ph2Idx = -1;
    float Smass = -1.0;
    float dR = -1.0;
    float pt1 = -1.0;
    float pt2 = -1.0;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);
    for (int iph1 = 0; iph1 < pt.size(); iph1++) {
     if (pt[iph1]>pt1) {
       pt1 = pt[iph1];
       ph1Idx = iph1;
     }
    }
    for (int iph2 = 0; (iph2 < pt.size()); iph2++) {
     if ((iph2!=ph1Idx)&&(pt[iph2]>pt2)) {
       pt2 = pt[iph2];
       ph2Idx = iph2;
     }
    }
    if ((ph1Idx>-1)&&(ph2Idx>-1)) {
     Lvector1 = Lvector[ph1Idx];
     Lsum.SetCoordinates(0,0,0,0);
     Lvector2 = Lvector[ph2Idx];
     Lsum = Lvector1+Lvector2;
     Smass = Lsum.M();
     dR = hardware::DeltaR(Lvector1,Lvector2);
    }
    return Smass;
}

float SmassCalcLeadingOrdered(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    int ph1Idx = -1;
    int ph2Idx = -1;
    float Smass = -1.0;
    float dR = -1.0;
    float pt1 = -1.0;
    float pt2 = -1.0;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);
    if (pt.size()>1) {
       pt1 = pt[0];
       ph1Idx = 0;
       pt2 = pt[1];
       ph2Idx = 1;
    }
    if ((ph1Idx>-1)&&(ph2Idx>-1)) {
     Lvector1 = Lvector[ph1Idx];
     Lsum.SetCoordinates(0,0,0,0);
     Lvector2 = Lvector[ph2Idx];
     Lsum = Lvector1+Lvector2;
     Smass = Lsum.M();
     dR = hardware::DeltaR(Lvector1,Lvector2);
    }
    return Smass;
}


float GENmassCalcLeading(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, RVec<int> id, int IdUse, RVec<int> idM, int IdMUse) {
    int ph1Idx = -1;
    int ph2Idx = -1;
    float Smass = -1.0;
    float dR = -1.0;
    float pt1 = -1.0;
    float pt2 = -1.0;
    int icurrM= -1;
    RVec<int> motherId=idM;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);

    for (int i = 0; i < idM.size(); i++) {
     icurrM = idM[i];
     if (icurrM>=0) {
      motherId[i] = id[icurrM];
     } else {
      motherId[i] = -1;
     }
    }

    for (int iph1 = 0; iph1 < pt.size(); iph1++) {
      if ((id[iph1]==IdUse)&&(motherId[iph1]==IdMUse)&&(pt[iph1]>pt1)) {
        pt1 = pt[iph1];
        ph1Idx = iph1;
      }
    }
//std::cout << "IdUse, ph1Idx, pt1, idM[ph1Idx] " << IdUse << " ;" << ph1Idx << "; " << pt1 << " " << idM[ph1Idx] << std::endl;
    for (int iph2 = 0; (iph2 < pt.size()); iph2++) {
//std::cout << "iph2, id[iph2], idM[iph2], (id[iph2]==IdUse), (idM[iph2]==idM[ph1Idx]) " << iph2 << " " << id[iph2] << " " <<  idM[iph2] << " " <<  (id[iph2]==IdUse) << " " << (idM[iph2]==idM[ph1Idx]) << std::endl;
      if ((iph2!=ph1Idx)&&(id[iph2]==IdUse)&&(idM[iph2]==idM[ph1Idx])) {
//      if ((id[iph2]==IdUse)&&(motherId[iph2]==IdMUse)&&(pt[iph2]>pt2)&&(idM[iph2]==idM[ph1Idx])) {
        pt2 = pt[iph2];
        ph2Idx = iph2;
      }
    }
    if ((ph1Idx>-1)&&(ph2Idx>-1)) {
     Lvector1 = Lvector[ph1Idx];
     Lsum.SetCoordinates(0,0,0,0);
     Lvector2 = Lvector[ph2Idx];
     Lsum = Lvector1+Lvector2;
     Smass = Lsum.M();
     dR = hardware::DeltaR(Lvector1,Lvector2);
    }
//std::cout << "pt1, pt2, Smass " << pt1 << " " << ph1Idx << "; " << pt2 << " " << ph2Idx << " ;" << Smass << std::endl;
    return Smass;
}

float dRCalc(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, float SmassTarget) {
    int ph0Idx = -1;
    int ph1Idx = -1; 
    float SmassTry = 0.0;
    float Smass = -1.0;
    float dSmass = 1000000.0;
    float dR = -1.0;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);
    for (int iph1 = 0; iph1 < pt.size(); iph1++) {
     Lvector1 = Lvector[iph1];
     for (int iph2 = iph1+1; iph2 < pt.size(); iph2++) {
       Lsum.SetCoordinates(0,0,0,0);
       Lvector2 = Lvector[iph2];
       Lsum = Lvector1+Lvector2;
       SmassTry = Lsum.M();
       if (abs(SmassTry-SmassTarget) < dSmass) {
          dSmass = abs(SmassTry-SmassTarget);
          ph0Idx = iph1;
          ph1Idx = iph2;
          Smass = SmassTry;
          dR = hardware::DeltaR(Lvector1,Lvector2);
       }
     }
    }
    return dR;
}

float dRCalcLeading(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    int ph1Idx = -1;
    int ph2Idx = -1;
    float Smass = -1.0;
    float dR = -1.0;
    float pt1 = -1.0;
    float pt2 = -1.0;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);
    for (int iph1 = 0; iph1 < pt.size(); iph1++) {
     if (pt[iph1]>pt1) {
       pt1 = pt[iph1];
       ph1Idx = iph1;
     }
    }
    for (int iph2 = 0; (iph2 < pt.size()); iph2++) {
     if ((iph2!=ph1Idx)&&(pt[iph2]>pt2)) {
       pt2 = pt[iph2];
       ph2Idx = iph2;
     }
    }
    if ((ph1Idx>-1)&&(ph2Idx>-1)) {
     Lvector1 = Lvector[ph1Idx];
     Lsum.SetCoordinates(0,0,0,0);
     Lvector2 = Lvector[ph2Idx];
     Lsum = Lvector1+Lvector2;
     Smass = Lsum.M();
     dR = hardware::DeltaR(Lvector1,Lvector2);
    }
    return dR;
}

float dRCalcLeadingOrdered(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass) {
    int ph1Idx = -1;
    int ph2Idx = -1;
    float Smass = -1.0;
    float dR = -1.0;
    float pt1 = -1.0;
    float pt2 = -1.0;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);
    if (pt.size()>1) {
       pt1 = pt[0];
       ph1Idx = 0;
       pt2 = pt[1];
       ph2Idx = 1;
    }
    if ((ph1Idx>-1)&&(ph2Idx>-1)) {
     Lvector1 = Lvector[ph1Idx];
     Lsum.SetCoordinates(0,0,0,0);
     Lvector2 = Lvector[ph2Idx];
     Lsum = Lvector1+Lvector2;
     Smass = Lsum.M();
     dR = hardware::DeltaR(Lvector1,Lvector2);
    }
    return dR;
}


float GENdRCalcLeading(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, RVec<int> id, int IdUse, RVec<int> idM, int IdMUse) {
    int ph1Idx = -1;
    int ph2Idx = -1;
    float Smass = -1.0;
    float dR = -1.0; 
    float pt1 = -1.0;
    float pt2 = -1.0;
    RVec<int> motherId=idM;
    int icurrM;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1;
    ROOT::Math::PtEtaPhiMVector Lvector2;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(pt,eta,phi,mass);

    for (int i = 0; i < idM.size(); i++) {
     icurrM = idM[i];
     if (icurrM>=0) {
      motherId[i] = id[icurrM];
     } else {
      motherId[i] = -1;
     }
    }

    for (int iph1 = 0; iph1 < pt.size(); iph1++) {
     if ((id[iph1]==IdUse)&&(motherId[iph1]==IdMUse)&&(pt[iph1]>pt1)) {
       pt1 = pt[iph1];
       ph1Idx = iph1;
     }
    }

    for (int iph2 = 0; (iph2 < pt.size()); iph2++) {
      if ((iph2!=ph1Idx)&&(id[iph2]==IdUse)&&(idM[iph2]==idM[ph1Idx])) {
//     if ((id[iph2]==IdUse)&&(motherId[iph2]==IdMUse)&&(pt[iph2]>pt2)&&(idM[iph2]==idM[ph1Idx])) {
       pt2 = pt[iph2];
       ph2Idx = iph2;
     }
    }
    if ((ph1Idx>-1)&&(ph2Idx>-1)) {
     Lvector1 = Lvector[ph1Idx];
     Lsum.SetCoordinates(0,0,0,0);
     Lvector2 = Lvector[ph2Idx];
     Lsum = Lvector1+Lvector2;
     Smass = Lsum.M();
     dR = hardware::DeltaR(Lvector1,Lvector2);
    }
    return dR;
}

RVec<float> PickLeadingLepton(int nElectron, int nMuon, RVec<float> Electron_pt, RVec<float> Muon_pt, RVec<float> Electron_eta, RVec<float> Muon_eta, RVec<float> Electron_phi, RVec<float> Muon_phi, RVec<float> Electron_mass, RVec<float> Muon_mass) { 
    float lpt=-1.0;
    float leta=-5.0;
    float lphi=-5.0;
    float lmass=-1.0;
    if (nElectron+nMuon>0) {
       if (nElectron>0) {
         lpt = Electron_pt[0];
         leta = Electron_eta[0];
         lphi = Electron_phi[0];
         lmass = Electron_mass[0];
       }
       if (nMuon>0) {
         if (Muon_pt[0]>lpt) {
             lpt = Muon_pt[0];
             leta = Muon_eta[0];
             lphi = Muon_phi[0];
             lmass = Muon_mass[0];
         }
       }
    }
    return {lpt,leta,lphi,lmass};
}

RVec<float> PickLeadingQJets(RVec<float> Jet_pt, RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Jet_mass, RVec<float> Jet_hadronFlavour, int Flavour_Choice) {
    float qpt0=-1.0;
    float qpt1=-1.0;
    float qeta0=-5.0;
    float qeta1=-5.0;
    float qphi0=-5.0;
    float qphi1=-5.0;
    float qmass0=-1.0;
    float qmass1=-1.0;
    int id0=-1;
    int id1=-1;
    for (int i = 0; i < Jet_pt.size(); i++) { 
       if ((Jet_hadronFlavour[i]==Flavour_Choice)&&(Jet_pt[i]>qpt0)) {
            qpt0=Jet_pt[i];
            id0=i;
       }
    }
    if (id0>=0) {
       qeta0=Jet_eta[id0];
       qphi0=Jet_phi[id0];
       qmass0=Jet_mass[id0];
       for (int j = 0; j < Jet_pt.size(); j++) {
          if ((j!=id0)&&(Jet_hadronFlavour[j]==Flavour_Choice)&&(Jet_pt[j]>qpt1)) {
               qpt1=Jet_pt[j];
               id1=j;
          }
       }
       if (id1>=0) {
         qeta1=Jet_eta[id1];
         qphi1=Jet_phi[id1];
         qmass1=Jet_mass[id1];
       }
    }
    return {qpt0,qpt1,qeta0,qeta1,qphi0,qphi1,qmass0,qmass1};  
}

RVec<float> PickLeadingDiJets(RVec<float> Jet_pt, RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Jet_mass, RVec<float> Jet_hadronFlavour, int Flavour_Choice) {
    float qpt0=-1.0;
    float qpt1=-1.0;
    float qeta0=-5.0;
    float qeta1=-5.0;
    float qphi0=-5.0;
    float qphi1=-5.0;
    float qmass0=-1.0;
    float qmass1=-1.0;
    int id0=-1;
    int id1=-1;
    for (int i = 0; i < Jet_pt.size(); i++) {
       if ((Jet_hadronFlavour[i]!=Flavour_Choice)&&(Jet_pt[i]>qpt0)) {
            qpt0=Jet_pt[i];
            id0=i;
       }
    }
    if (id0>=0) {
       qeta0=Jet_eta[id0];
       qphi0=Jet_phi[id0];
       qmass0=Jet_mass[id0];
       for (int j = 0; j < Jet_pt.size(); j++) {
          if ((j!=id0)&&(Jet_hadronFlavour[j]!=Flavour_Choice)&&(Jet_pt[j]>qpt1)) {
               qpt1=Jet_pt[j];
               id1=j;
          }
       }
       if (id1>=0) {
         qeta1=Jet_eta[id1];
         qphi1=Jet_phi[id1];
         qmass1=Jet_mass[id1];
       }
    }
    return {qpt0,qpt1,qeta0,qeta1,qphi0,qphi1,qmass0,qmass1};
}

RVec<float> PickLeadingDiJetsDR(RVec<float> bjet_vec, float dRselection, RVec<float> Jet_pt, RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Jet_mass, RVec<float> Jet_hadronFlavour, int Flavour_Choice) {
    float qpt0=-1.0;
    float qpt1=-1.0;
    float qeta0=-5.0;
    float qeta1=-5.0;
    float qphi0=-5.0;
    float qphi1=-5.0;
    float qmass0=-1.0;
    float qmass1=-1.0;
    int id0=-1;
    int id1=-1;
    float b1pt=bjet_vec[0];
    float b2pt=bjet_vec[1];
    float b1eta=bjet_vec[2];
    float b2eta=bjet_vec[3];
    float b1phi=bjet_vec[4];
    float b2phi=bjet_vec[5];
    float b1mass=bjet_vec[6];
    float b2mass=bjet_vec[7];
    RVec<float> dR1=Jet_pt;
    RVec<float> dR2=Jet_pt;
    ROOT::Math::PtEtaPhiMVector Lsum;
    ROOT::Math::PtEtaPhiMVector Lvector1=hardware::TLvector(b1pt,b1eta,b1phi,b1mass);
    ROOT::Math::PtEtaPhiMVector Lvector2=hardware::TLvector(b2pt,b2eta,b2phi,b2mass);;
    RVec<ROOT::Math::PtEtaPhiMVector> Lvector=hardware::TLvector(Jet_pt,Jet_eta,Jet_phi,Jet_mass);
    for (int i = 0; i < Jet_pt.size(); i++) {
      dR1[i]=hardware::DeltaR(Lvector[i],Lvector1);
      dR2[i]=hardware::DeltaR(Lvector[i],Lvector2);
    }
    for (int i = 0; i < Jet_pt.size(); i++) {
     if ((dR1[i]>dRselection)&&(dR2[i]>dRselection)) {
       if ((Jet_hadronFlavour[i]!=Flavour_Choice)&&(Jet_pt[i]>qpt0)) {
            qpt0=Jet_pt[i];
            id0=i;
       }
     }
    }
    if (id0>=0) {
       qeta0=Jet_eta[id0];
       qphi0=Jet_phi[id0];
       qmass0=Jet_mass[id0];
       for (int j = 0; j < Jet_pt.size(); j++) {
         if ((dR1[j]>dRselection)&&(dR2[j]>dRselection)) {
          if ((j!=id0)&&(Jet_hadronFlavour[j]!=Flavour_Choice)&&(Jet_pt[j]>qpt1)) {
               qpt1=Jet_pt[j];
               id1=j;
          }
         }
       }
       if (id1>=0) {
         qeta1=Jet_eta[id1];
         qphi1=Jet_phi[id1];
         qmass1=Jet_mass[id1];
       }
    }
//if ((id0>=0)&&(id1>=0)) {
//std::cout << "final dijet choice " << id0 << " dR1 " << dR1[id0] << " " << id1 << " dR1 " << dR2[id1]  << std::endl;
//} else {
//std::cout << "nothing found"<< std::endl;
//}

    return {qpt0,qpt1,qeta0,qeta1,qphi0,qphi1,qmass0,qmass1};
}

RVec<float> PickLeadingDiPhotons(RVec<float> Photon_pt, RVec<float> Photon_eta, RVec<float> Photon_phi, RVec<float> Photon_mass, RVec<float> Photon_mvaID) {
    float qpt0=-1.0;
    float qpt1=-1.0;
    float qeta0=-5.0;
    float qeta1=-5.0;
    float qphi0=-5.0;
    float qphi1=-5.0;
    float qmass0=-1.0;
    float qmass1=-1.0;
    float qmvaID0=-5.0;
    float qmvaID1=-5.0;
    if (Photon_pt.size()>0) {
       qpt0=Photon_pt[0];
       qeta0=Photon_eta[0];
       qphi0=Photon_phi[0];
       qmass0=Photon_mass[0];
       qmvaID0=Photon_mvaID[0];
       if (Photon_pt.size()>1) {
         qpt1=Photon_pt[1];
         qeta1=Photon_eta[1];
         qphi1=Photon_phi[1];
         qmass1=Photon_mass[1];
         qmvaID1=Photon_mvaID[1];
       }
    }
    return {qpt0,qpt1,qeta0,qeta1,qphi0,qphi1,qmass0,qmass1,qmvaID0,qmvaID1};
}

RVec<int> PickLeadingDiPhotonsI(RVec<float> Photon_pt,RVec<int> Photon_cutBased) {
    int qcutBased0=-1;
    int qcutBased1=-1;
    if (Photon_pt.size()>0) {
       qcutBased0=Photon_cutBased[0];
       if (Photon_pt.size()>1) {
         qcutBased1=Photon_cutBased[1];
       }
    }
    return {qcutBased0,qcutBased1};
}

std::vector<float> ConvertToVecF(float f1, float f2) {
    std::vector<float> fout(2);
    fout[0]=f1;
    fout[1]=f2;
    return fout;
}

std::vector<int> PickTagF(RVec<float> tagScore, RVec<int> idxs, float scoreCut) {
    if (idxs.size()>2) {
        std::cout << "PickTagF -- WARNING: You have input more than two indices. Only two accepted. Assuming first two indices.";
    }
    std::vector<int> out(2);
    float WP = scoreCut;

    int idx0 = idxs[0];
    int idx1 = idxs[1];
    bool isTop0, isTop1;
    isTop0 = (WP < tagScore[idx0]);
    isTop1 = (WP < tagScore[idx1]);

    if (isTop0 && isTop1) {
        out[0] = idx0;
        out[1] = idx1;
    } else if (isTop0) {
        out[0] = idx0;
        out[1] = idx1;
    } else if (isTop1) {
        out[0] = idx1;
        out[1] = idx0;
    } else {
        out[0] = -1;
        out[1] = -1;
    }
    return out;
}

std::vector<int> PickTagI(RVec<int> tagScore, RVec<int> idxs, int scoreCut) {
    if (idxs.size()>2) {
        std::cout << "PickTagI -- WARNING: You have input more than two indices. Only two accepted. Assuming first two indices.";
    }
    std::vector<int> out(2);
    int IWP = scoreCut;
    int idx0 = idxs[0];
    int idx1 = idxs[1];
    bool isTop0, isTop1;

    isTop0 = (IWP < tagScore[idx0]);
    isTop1 = (IWP < tagScore[idx1]);
    if (isTop0 && isTop1) {
        out[0] = idx0;
        out[1] = idx1;
    } else if (isTop0) {
        out[0] = idx0;
        out[1] = idx1;
    } else if (isTop1) {
        out[0] = idx1;
        out[1] = idx0;
    } else {
        out[0] = -1;
        out[1] = -1;
    }
    return out;
}

std::vector<int> PickTop(RVec<float> mass, RVec<float> tagScore, RVec<int> idxs, std::pair<float,float> massCut, float scoreCut, bool invertScore=false) {
    if (idxs.size()>2) {
        std::cout << "PickTop -- WARNING: You have input more than two indices. Only two accepted. Assuming first two indices.";
    }
    std::vector<int> out(2);
    float WP = scoreCut;

    int idx0 = idxs[0];
    int idx1 = idxs[1];
    bool isTop0, isTop1;
    if (!invertScore) {
        isTop0 = (mass[idx0] > massCut.first) && (mass[idx0] < massCut.second) && (tagScore[idx0] > WP);
        isTop1 = (mass[idx1] > massCut.first) && (mass[idx1] < massCut.second) && (tagScore[idx1] > WP);
    } else {
	// if inverted, only accept jets meeting the loose score (>0.2) and not those meeting the tight cut (< WP)
        isTop0 = (tagScore[idx0] < WP) && (0.2 < tagScore[idx0]);
	isTop1 = (tagScore[idx1] < WP) && (0.2 < tagScore[idx1]);
        //isTop0 = (mass[idx0] > massCut.first) && (mass[idx0] < massCut.second) && (tagScore[idx0] < WP);
        //isTop1 = (mass[idx1] > massCut.first) && (mass[idx1] < massCut.second) && (tagScore[idx1] < WP);
    }
    
    if (isTop0 && isTop1) {
        if (tagScore[idx0] > tagScore[idx1]) {
            out[0] = idx0;
            out[1] = idx1;
        } else {
            out[0] = idx1;
            out[1] = idx0;
        }
    } else if (isTop0) {
        out[0] = idx0;
        out[1] = idx1;
    } else if (isTop1) {
        out[0] = idx1;
        out[1] = idx0;
    } else {
        out[0] = -1;
        out[1] = -1;
    }
    return out;
}

std::vector<int> PickTop_Simple(RVec<float> tagScore, RVec<int> idxs, float scoreCut, bool invertScore=false) {
    if (idxs.size()>2) {
        std::cout << "PickTop -- WARNING: You have input more than two indices. Only two accepted. Assuming first two indices.";
    }
    std::vector<int> out(2);
    float WP = scoreCut;

    int idx0 = idxs[0];
    int idx1 = idxs[1];
    bool isTop0, isTop1;
    if (!invertScore) {
//        isTop0 = (mass[idx0] > massCut.first) && (mass[idx0] < massCut.second) && (tagScore[idx0] > WP);
//        isTop1 = (mass[idx1] > massCut.first) && (mass[idx1] < massCut.second) && (tagScore[idx1] > WP);
        isTop0 = (tagScore[idx0] > WP);
        isTop1 = (tagScore[idx1] > WP);
    } else {
        // if inverted, only accept jets meeting the loose score (>0.2) and not those meeting the tight cut (< WP)
        isTop0 = (tagScore[idx0] < WP) && (0.2 < tagScore[idx0]);
        isTop1 = (tagScore[idx1] < WP) && (0.2 < tagScore[idx1]);
    }

    if (isTop0 && isTop1) {
        if (tagScore[idx0] > tagScore[idx1]) {
            out[0] = idx0;
            out[1] = idx1;
        } else {
            out[0] = idx1;
            out[1] = idx0;
        }
    } else if (isTop0) {
        out[0] = idx0;
        out[1] = idx1;
    } else if (isTop1) {
        out[0] = idx1;
        out[1] = idx0;
    } else {
        out[0] = -1;
        out[1] = -1;
    }
    return out;
}


std::vector<int> PickTopCRv2(RVec<float> mass, RVec<float> tagScore, RVec<float> HiggsScore, RVec<int> idxs, std::pair<float,float> massCut, float scoreCut, bool invertScore=false) {
    if (idxs.size() > 2) {
	std::cout << "PickTop -- WARNING: You have input more than two indices. Only two accepted. Assuming first two indices.";
    }
    std::vector<int> out(2);
    float WP = scoreCut;
    int idx0 = idxs[0];
    int idx1 = idxs[1];
    bool isTop0, isTop1;
    if (!invertScore) {	// signal region - apply mass window cut and top tag reqirement
        isTop0 = (mass[idx0] > massCut.first) && (mass[idx0] < massCut.second) && (tagScore[idx0] > WP);
        isTop1 = (mass[idx1] > massCut.first) && (mass[idx1] < massCut.second) && (tagScore[idx1] > WP);
    } else {	// control region - anti-top tag and Higgs veto on the top jet
	// same as CR_v1, but also require the top candidate jet to have Higgs tag < loose WP
	isTop0 = (tagScore[idx0] < WP) && (0.2 < tagScore[idx0]) && (HiggsScore[idx0] < 0.2);
        isTop1 = (tagScore[idx1] < WP) && (0.2 < tagScore[idx1]) && (HiggsScore[idx1] < 0.2);
    }
    if (isTop0 && isTop1) {	// if both jets meet top (anti)tag, then choose the one with higher top score as top jet
        if (tagScore[idx0] > tagScore[idx1]) {
            out[0] = idx0;
            out[1] = idx1;
        } else {
            out[0] = idx1;
            out[1] = idx0;
        }
    // otherwise, chose the one that passed tagging to be the top candidate
    } else if (isTop0) {
        out[0] = idx0;
        out[1] = idx1;
    } else if (isTop1) {
        out[0] = idx1;
        out[1] = idx0;
    } else {
        out[0] = -1;
        out[1] = -1;
    }
    return out;
}


bool MatchToGen(int pdgID, ROOT::Math::PtEtaPhiMVector jet,
                RVec<ROOT::Math::PtEtaPhiMVector> GenPart_vect,
                RVec<int> GenPart_pdgId) {
    bool found = false;
    for (int igp = 0; igp<GenPart_vect.size(); igp++) {
        if (abs(GenPart_pdgId[igp]) == pdgID) {
            if (hardware::DeltaR(jet,GenPart_vect[igp]) < 0.8) {
                found = true;
                break;
            }
        }
    }
    return found;
}

std::vector<int> PickTopGenMatch(RVec<ROOT::Math::PtEtaPhiMVector> Dijet_vect,
                                 RVec<ROOT::Math::PtEtaPhiMVector> GenPart_vect,
                                 RVec<int> GenPart_pdgId) {
    if (Dijet_vect.size()>2) {
        std::cout << "PickTopGenMatch -- WARNING: You have input more than two indices. Only two accepted. Assuming first two indices.";
    }
    int tIdx = -1;
    int hIdx = -1;
    for (int ijet = 0; ijet < 2; ijet++) {
        if (MatchToGen(6, Dijet_vect[ijet], GenPart_vect, GenPart_pdgId)) {
            tIdx = ijet;
        } else if (MatchToGen(25, Dijet_vect[ijet], GenPart_vect, GenPart_pdgId)) {
            hIdx = ijet;
        }
    }
    return {tIdx,hIdx};
}

RVec<int> PickJetAnd2Gamma(RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> mass, RVec<float> photon_phi) {
    int jetIdx = -1;
    int phoIdx1 = -1;
    int phoIdx2 = -1;
    for (int ijet=0; ijet < pt.size(); ijet++) {
// DP EDIT
//        if (pt[ijet] > 350 && std::abs(eta[ijet]) < 2.4 && mass[ijet] > 50) {
        if (pt[ijet] > 350 && std::abs(eta[ijet]) < 3.0 && mass[ijet] > 0) {
            jetIdx = ijet;
            break;      // we've found one candidate jet
        }
        // now check all the surrounding photons
        if (jetIdx < 0) {
            return {-1,-1,-1};
        }
        else {
            for (int ipho=0; ipho < photon_phi.size(); ipho++) {
                if (phoIdx2 == -1) {
                    if (abs(hardware::DeltaPhi(phi[jetIdx], photon_phi[ipho])) < M_PI/4) {
                        if (phoIdx1 == -1) {
                            phoIdx1 = ipho;
                        }
                        else {
                            if (abs(hardware::DeltaPhi(phi[jetIdx], photon_phi[ipho])) < M_PI/4) {
                                phoIdx2 = ipho;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    return {jetIdx, phoIdx1, phoIdx2};
}

RVec<int> PickDijetsNminus1(RVec<float> phi) {
// used for Nminus1 testing (only kinematic cut here will be deltaPhi)
    int jet0Idx = -1;
    int jet1Idx = -1;
    for (int ijet=0; ijet<phi.size(); ijet++) {
	if (jet1Idx == -1) {
	    if (jet0Idx == -1) {
		jet0Idx = ijet;
	    } else {
	    if (abs(hardware::DeltaPhi(phi[jet0Idx], phi[ijet])) > M_PI/2) {
		jet1Idx = ijet;
		break;
	    }
	    }
	}
    }
    return {jet0Idx,jet1Idx};
}
