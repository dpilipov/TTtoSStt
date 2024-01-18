#include<ROOT/RVec.hxx>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <vector>

/*
 * Class for applying ParticleNet scale factors on a jet-by-jet basis. 
 * For use in the T' -> t\phi analysis, with scale factors located at:
 * https://coli.web.cern.ch/coli/.cms/btv/boohft-calib/20220623_bb_TprimeB_useExpr_2016/4_fit/
 * 
 * Method: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#2a_Jet_by_jet_updating_of_the_b
 * There are a few things for consideration:
 *  1) You have to run this for each working point, so that jets that were reassigned for the first WP are used for the second.
 *     To do so, create a new column using TIMBER analyzer's Define() method and pass the new column to the updateTag() method.
 *  2) The efficiencies do not always obey e_hp < e_mp < 1, so this has to be accounted for in the equality check functions.
*/

using namespace ROOT::VecOps;

class PNetSaaSFHandler {
  private:
    float _wp;     // MP [0.8, 0.98], HP [0.98, 1.0]
    float _eff;    // efficiencies will be calculated via TIMBER then fed to constructor
    std::string _year;    // 2016APV, 2016, 2017, 2018
    int _var;             // 0: nominal, 1: up, 2: down, passed to constructor
    TRandom _rand;        // used for random number generation
    int _newTags[2]  = {0,0};      // number of photons in each new category [fail][pass]
    int _origTags[2] = {0,0};      // original num photons in each category [fail][pass]
    
    // SF[_var][pt]
    // variations are described above, pt cats are [400, 600), [600, 800), [800, +inf) across all years
    // HP (tight) [0.98, 1.0]   PASS
    float SF2016APV_T[3][3] = {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};
    float SF2016_T[3][3]    = {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};
    float SF2017_T[3][3]    = {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};
    float SF2018_T[3][3]    = {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};
    // MP (loose) [0.8, 0.98]  FAIL
    float SF2016APV_L[3][3] = {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};
    float SF2016_L[3][3]    = {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};
    float SF2017_L[3][3]    = {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};
    float SF2018_L[3][3]    = {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};

  public:
    PNetSaaSFHandler(float wp, float eff, std::string year, int var);  // default: wps={-0.995,-0.9}, effs={effl,effT}, var=0/1/2
    ~PNetSaaSFHandler();
    int getWPcat(float taggerVal);                                    // determine WP category: 0: fail, 1: loose, 2: tight
    float getSF(float pt, float taggerVal);                           // gets the proper SF based on jet's pt and score as well as internal variables _year, _var
    int updateTag(int photonCat, float pt, float taggerVal);	      // determines the jet's new tagger category
    int createTag(float taggerVal);				      // create new tagger category based on jet's original tagger value
    void printVals();						      // print the number of jets in each category

    // These functions should *NOT* be used, but are included for posterity
    // Essentially, instead of generating a vector of integers, just generate an integer
    // The reason is that we call these via TIMBER's Define() function, which takes in a C++ function or string and applies it to *EVERY* value in the column 
    // This means that we don't have to pass in the vectors themselves, just the column names
    RVec<int> updateTag(RVec<int> jetCats, RVec<float> pt, RVec<float> taggerVals);   // determines the jet's new tagger category 
    RVec<int> createTag(RVec<float> taggerVals);                                      // create vector of tagger categories based on jets' original tagger value.
};

void PNetSaaSFHandler::printVals() {
    // prints the number of original and new tagger values
    printf("Number of Original\n\tFail: %d\n\tPass: %d\n\tTotal: %d\n", _origTags[0], _origTags[1], _origTags[0]+_origTags[1]);
    printf("Number of New\n\tFail: %d\n\tPass: %d\n\tTotal: %d\n", _newTags[0], _newTags[1], _newTags[0]+_newTags[1]);
};

PNetSaaSFHandler::PNetSaaSFHandler(float wp, float eff, std::string year, int var) {
  _wp = wp;
  _eff = eff;
  _year = year;
  _var = var;
  // unique but repeatable random numbers. For repeated calls in the same event, random #s from Rndm() will be identical
  _rand = TRandom(1234);
};

PNetSaaSFHandler::~PNetSaaSFHandler() {
    // print out number of original and new tagger vals upon destruction
    printVals();
};

int PNetSaaSFHandler::getWPcat(float taggerVal) {
  // determine the WP category we're in, 0:fail, 1:loose, 2:tight
  int wpCat;
  if (taggerVal > _wp) { // pass
    wpCat = 1;
  }
  else {  // fail
    wpCat = 0;
  }
  return wpCat;
};

float PNetSaaSFHandler::getSF(float pt, float taggerVal) {
  /* getthe scale factor from the jet's year, score, and pt */
  float SF;
  int ptCat;
  int wpCat = getWPcat(taggerVal);
  // get the pT category
  if ((pt >= 400) && (pt < 600)) {
    ptCat = 0;
  }
  else if ((pt >= 600) && (pt < 800)) {
    ptCat = 1;
  }
  else if (pt > 800) {
    ptCat = 2;
  }
  else {
    // jet is outside of the pt range used in SF derivation, return no change
    return 1.0;
  }
  // get the SF
  switch (wpCat) {
    case 0:   // if jet is originally in fail, pass SF of 1.0 (no change)
      SF = 1.0;
    case 1:   // jet is in MP (loose)
      if (_year=="2016APV") { SF = SF2016APV_L[_var][ptCat]; }
      else if (_year=="2016") { SF = SF2016_L[_var][ptCat]; }
      else if (_year=="2017") { SF = SF2017_L[_var][ptCat]; }
      else { SF = SF2018_L[_var][ptCat]; }
    case 2:   // jet is in HP (tight)
      if (_year=="2016APV") { SF = SF2016APV_T[_var][ptCat]; }
      else if (_year=="2016") { SF = SF2016_T[_var][ptCat]; }
      else if (_year=="2017") { SF = SF2017_T[_var][ptCat]; }
      else { SF = SF2018_T[_var][ptCat]; }
  }
  return SF;
};


int PNetSaaSFHandler::createTag(float taggerVal) {
    /* Creates tagger categories for phi candidate jets
     * this MUST be called in TIMBER before running the rest of the script, as it places all jets into their respective categories for later use in updateTag()
     * This function is meant to be called after selecting the top and higgs in CR and SR (see THselection.py - getEfficiencies)
    */
    if (taggerVal > _wp) {
	_origTags[1]++;
	return 1;
    }
    else {
	_origTags[0]++;
	return 0;
    }
};

int PNetSaaSFHandler::updateTag(int photonCat, float pt, float taggerVal) {
    /* updates the tagger category for phi jets
     * https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#2a_Jet_by_jet_updating_of_the_b
     * Params:
     * 	 jetCat    = original jet category (0: fail, 1: loose, 2: pass)
     * 	 pt        = jet pt
     * 	 taggerVal = particleNet tagger value
    */ 
    float eff = _eff;
    float SF = getSF(pt, taggerVal);
    float SF_T = SF;
    float SF_L = SF;
    float eff_T = eff;
    float eff_L = eff;
    double rn = _rand.Rndm();
    int newCat = photonCat;	// grab the original tag category, will be updated
    // begin logica
    if ((SF_L < 1) && (SF_T < 1)) {
	if ( (newCat==2) && (rn < (1.-SF_T)) ) newCat=0;	// tight (2) -> untag (0)
	if ( (newCat==1) && (rn < (1.-SF_L)) ) newCat=0;	// loose (1) -> untag (0)
    }
    if ((SF_L > 1) && (SF_T > 1)) {
	float fL, fT;
	if (newCat==0) {
	    float num = eff_T*(SF_T-1.);
	    float den = 1.-(eff_L+eff_T);
	    fT = num/den;
	    if (rn < fT) newCat=2;	// untag (0) -> tight (2)
	    else {
		rn = _rand.Rndm();
		num = eff_L*(SF_L-1.);
		den = (1.-eff_L-eff_T)*(1.-fT);
		fL = num/den;
		if (rn < fL) newCat=1; 	// loose (1) -> tight (2)
	    }
	}
    }
    if ((SF_L < 1) && (SF_T > 1)) {
	if (newCat==0) {
	    float num = eff_T*(SF_T-1.);
	    float den = 1.-(eff_L+eff_T);
	    float f = num/den;
	    if (rn < f) newCat=2;	    // untag (0) -> tight (2)
	}
	if (newCat==1) {
	    if (rn < (1.-SF_L)) newCat=0;   // loose (1) -> untag (0)
	}
    }
    if ((SF_L > 1) && (SF_T < 1)) { 
	if ((newCat==2) && (rn < (1.-SF_T))) newCat=0;	// tight (2) -> untag (0)
	if (newCat==0) {
	    float num = eff_L*(SF_L-1.);
	    float den = 1.-(eff_L+eff_T);
	    float f = num/den;
	    if (rn < f) newCat=1;	// untag (0) -> loose (1)
	}
    }
    // update the new tag category array
    _newTags[newCat]++;
    // return the new tagger category value
    return newCat;
};

/* ----------------------------------------------------------------------------------
 *
 *    Do not use these functions, use the above versions 
 *
 * ----------------------------------------------------------------------------------
 */
RVec<int> PNetSaaSFHandler::createTag(RVec<float> taggerVals) {
  /* Creates tagger categories for phi candidate jets.
   * This MUST be called in TIMBER before running the rest of the script, as it places all jets into their respective categories for later use in updateTag()
   * example calling from TIMBER (after compiling class): analyzer.Define("ScaledPnetH","PNetXbbSFHandler.createTag(particleNetMD_HbbvsQCD);")
   * cat (int): 0-fail, 1-MP(loose), 2-HP(tight)
  */
  //std::cout << "in vec size: " << taggerVals.size() << "\n"; 
  //printf("Creating tag categories - 0: Fail, 1: Loose, 2: Tight\n");
  RVec<int> jetCats(taggerVals.size());
  for (size_t ijet=0; ijet<taggerVals.size(); ijet++) {   // loop over all jets
    int cat;
    float taggerVal = taggerVals[ijet];
    if (taggerVal > _wp) {   // 0.8 < tag < 0.98
      _origTags[1]++;
      jetCats[ijet] = 1;
    }
    else {    // tag < 0.8
      _origTags[0]++;
      jetCats[ijet] = 0;
    }
  }
  //printf("Finished creating tag categories. Initial values before btag reassignment:\n\tFail: %i\n\tLoose: %i\n\tTight: %i\n",_origTags[0],_origTags[1],_origTags[2]);
  return jetCats;
};

RVec<int> PNetSaaSFHandler::updateTag(RVec<int> jetCats, RVec<float> pt, RVec<float> taggerVals) {
  /* https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#2a_Jet_by_jet_updating_of_the_b 
   * params:
   *    jetCats = vector of ints representing the category into which each jet has been placed (0:fail, 1:loose, 2:tight)
   *        - should be passed the RVec created by createTag()
   *    pt = vector of floats holding each jet's pt
   * returns:
   *    cats = new vector of ints representing the jet categories after checking the four SF conditions
  */
  //printf("Updating tag categories - 0: Fail, 1: Loose, 2: Tight\n");
  RVec<int> cats(jetCats.size());
  float eff = _eff;
  float eff_T = eff;
  float eff_L = eff;
  for (size_t ijet=0; ijet<pt.size(); ijet++) {
    // get the SF for loose and tight using the jet's pt, tagger value. The getSF() function uses the internal year value and calculates the tagger WP
    // pt, taggerVals, and jetCats should be same length, so ijets should work for indexing
    float SF_L = getSF(pt[ijet], taggerVals[ijet]);
    float SF_T = getSF(pt[ijet], taggerVals[ijet]);
    int newCat = jetCats[ijet];   // grab the original tag category

    // generate the random number for the event 
    double rn = _rand.Rndm();
    //printf("\trn: %f",rn);

    /* BROKEN - results in very incorrect retag values across nom/up/down
    // check the four cases and modify the new category appropriately
    if ((SF_L < 1) && (SF_T < 1)) {
      if ( (newCat==2) && (rn < (1-SF_T)) ) newCat--;                         // demote from tight (2) to loose (1)
      if ( (newCat==1) && (rn < (1-SF_L)/(1-(eff_T/eff_L)*SF_T)) ) newCat--;  // demote from loose (1) to untag (0)
    }
    else if ((SF_L > 1) && (SF_T > 1)) {
      if ( (newCat==0) && (rn < (SF_L-1)/((1./eff_L)-1)) ) newCat++;          // promote from untag (0) to loose (1)
      if ( (newCat==1) && (rn < (SF_T-1)/((eff_L/eff_T)*SF_T-1)) ) newCat++;  // promote from loose (1) to tight (2)
    }
    else if ((SF_L < 1) && (SF_T > 1)) {
      if ( (newCat==0) && (rn < (eff_T*(SF_T-1))/(1-(eff_T+eff_L)) ) ) newCat++;   // promote from loose (1) to tight (2)
      if ( (newCat==1) && (rn < (1-SF_T)) ) newCat--;                              // demote from loose (1) to untag(0)
    }
    else if ((SF_L > 1) && (SF_T < 1)) {
      if ( (newCat==2) && (rn < (1-SF_T)) ) newCat--;                              // demote from tight (2) to loose (1)
      if ( (newCat==0) && (rn < (eff_L*(SF_L-1))/(1-(eff_T+eff_L))) ) newCat++;    // promote from untag (0) to loose (1)
    }
    */

    if ((SF_L < 1) && (SF_T < 1)) {
      if ( (newCat==2) && (rn < (1.-SF_T)) ) newCat--;                         // demote from tight (2) to loose (1)
      if ( (newCat==1) && (rn < (1.-SF_L)) ) newCat--;                         // demote from loose (1) to untag (0)
    }
    else if ((SF_L > 1) && (SF_T > 1)) {
      if (newCat==0) {
        if (rn < (eff_T*(SF_T-1.))/(1.-eff_L+eff_T)) newCat=2;                 // promote from untag (0)  to tight (2)
        if (rn < (eff_L*(SF_L-1.))/(1.-eff_L-eff_T)*(1.-(eff_T*(SF_T-1.))/(1.-eff_L+eff_T)) ) newCat++;   // promote from untag (0) to loose (1)
      }
    }
    else if ((SF_L < 1) && (SF_T > 1)) {
      if ( (newCat==0) && (rn < (eff_T*(SF_T-1.))/(1.-(eff_L+eff_T)) ) ) newCat=2;  // promote from untag (0) to tight (2)
      if ( (newCat==1) && (rn < 1.-SF_L) ) newCat--;                                // demote from loose (1) to untag (0)    
    }
    else if ((SF_L > 1) && (SF_T < 1)) {
      if ( (newCat==2) && (rn <  1.-SF_T) ) newCat=0;                               // demote from tight (2) to untag (0)
      if ( (newCat==0) && (rn < (eff_L*(SF_L-1.))/(1.-(eff_L+eff_T)) ) ) newCat++;  // promote from untag (0) to loose (1)
    }

    // append new category value to RVec
    cats[ijet] = newCat;
    // update the new tag category array
    _newTags[newCat]++;
  }
  // values in _newTags array MUST be updated during the 4 above subroutines based on their outcome
  //printf("Finished updating tag categories. New values after btag reassignment:\n\tFail: %i\n\tLoose: %i\n\tTight: %i\n",_newTags[0],_newTags[1],_newTags[2]);
  return cats;
};
