#define AddBackTestSelector_cxx
// The class definition in AddBackTestSelector.h has been generated automatically
#include "AddBackTestSelector.h"

int kValue_grif = 379; // default k-value for not piled-up hits
int kValue_zds = 59; // default k-value for not piled-up hits
int kValue_paces = 200; // default k-value for not piled-up hits //paces has different values for different detectors it seems 
//flags
bool PileUpRej	 = false; // If true, pile-up is rejected
bool AddbOn		 = false; // If true, addback spectra are sorted
bool SuppOn		 = true; // If true, suppressed spectra are sorted
bool AddbSuppOn	 = false; // If true, addback-suppressed spectra are sorted
// Ancillaries
bool BgoOn		 = true; // If true, spectra and matrices for BGos are sorted
bool ZdsOn		 = true; // If true, spectra and matrices for Zds are sorted
bool SceptarOn	 = false; // If true, spectra and matrices for Sceptar are sorted (ONLY IF ZDS IS ALSO ON!)
bool LabrOn		 = false; // If true, spectra and matrices for LaBr are sorted
bool LabrSuppOn	 = false; // If true, spectra and matrices for LaBr with suppression are sorted
bool PacesOn	 = false; // If true, spectra and matrices for Paces are sorted


////Global variables for the time coincidence windows

double ggtime   = 120.;//Time coincidence windows for GRIFFIN-GRIFFIN, assumed symmetrical
double gztime_l = -200.;//Time coincidence windows for GRIFFIN-ZDS, lower limit
double gztime_u = 0.;//Time coincidence windows for GRIFFIN-ZDS, upper limit
double gptime_l = -150.;//Time coincidence windows for GRIFFIN-PACES, lower limit
double gptime_u = -50.;//Time coincidence windows for GRIFFIN-PACES, upper limit

double randtime_l=1000.;//Random time coincidence window, common for all coincidences, lower limit
double randtime_u=2000.;////Random time coincidence window, common for all coincidences, upper limit

double gz_randtime_l=800.;//Random time coincidence window, common for all coincidences, lower limit
double gz_randtime_u=1800.;////Random time coincidence window, common for all coincidences, upper limit

double gp_randtime_l=200.;//Random time coincidence window, common for all coincidences, lower limit
double gp_randtime_u=1200.;////Random time coincidence window, common for all coincidences, upper limit

double MyCycleLength = 783.5;

// prompt coincidences
bool PromptCoincidence(TGriffinHit *h1, TGriffinHit *h2){  //Check if GRIFFIN hits are less then 120 ns apart.
   return std::fabs(h1->GetTime() - h2->GetTime()) < ggtime;
}
bool PromptCoincidence(TGriffinHit *h1, TZeroDegreeHit *h2){  // coincidence window GRIFFIN-ZDS: (-450,-325) ns
   return gztime_l < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < gztime_u;
}
bool PromptCoincidence(TGriffinHit *h1, TPacesHit *h2){  // coincidence window GRIFFIN-PACES: (-450,-325) ns
   return gptime_l < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < gptime_u;
}

// time random
bool TimeRandom(TGriffinHit *h1, TGriffinHit *h2){  //Check if hits are less then 1000 ns apart.
   return randtime_l < std::fabs(h1->GetTime() - h2->GetTime()) && std::fabs(h1->GetTime() - h2->GetTime()) < randtime_u;
}
bool TimeRandom(TGriffinHit *h1, TZeroDegreeHit *h2){  //Check if hits are less then 1000 ns apart.
   return gz_randtime_l < std::fabs(h1->GetTime() - h2->GetTime()) && std::fabs(h1->GetTime() - h2->GetTime()) < gz_randtime_u;
}
bool TimeRandom(TGriffinHit *h1, TPacesHit *h2){  //Check if hits are less then 1000 ns apart.
   return gp_randtime_l < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < gp_randtime_u;
}


bool IsOpposite(TGriffinHit* Hit1, TGriffinHit* Hit2)
{
   TVector3 Vec1 = TGriffin::GetPosition(Hit1->GetDetector(), Hit1->GetCrystal(), 90);
   TVector3 Vec2 = TGriffin::GetPosition(Hit2->GetDetector(), Hit2->GetCrystal(), 90);
   return TMath::RadToDeg()*(Vec1.Angle(Vec2)) > 179.5;
} 

void AddBackTestSelector::CreateHistograms() {
   // get cycle length from ODB - this is in nanoseconds! (Actually it's not in ns, I got it in mu s!)
   //fCycleLength = fPpg->OdbCycleLength();
   //std::cout<<"got ODB cycle length "<<fCycleLength<<" mu s, "<<fCycleLength/1e6<<" s"<<std::endl;
	
////////////////////////////////////////////////////////////////////////////////////////////////////

	int NBins = 12000;
	double FirstBin = 0.25; double LastBin = 6000.25;
	//singles plots
	fH1["gMult"] = new TH1D("gMult","singles multiplicity", 100, -0.5, 99.5);
	fH1["aMult"] = new TH1D("aMult","addback multiplicity", 100, -0.5, 99.5);
	
	fH1["aFrags"] = new TH1D("aFrags","addback fragments", 100, -0.5, 99.5);
	
	fH2["gEA"] = new TH2D("gEA","bare #gamma singles", 64, 0.5, 64.5, NBins, FirstBin, LastBin);
	fH2["aEA"] = new TH2D("aEA","addback #gamma singles", 64, 0.5, 64.5, NBins, FirstBin, LastBin);
	
	fH2["aF1EA"] = new TH2D("aF1EA","addback #gamma singles One Frag", 64, 0.5, 64.5, NBins, FirstBin, LastBin);
	fH2["aF2EA"] = new TH2D("aF2EA","addback #gamma singles Two Frag", 64, 0.5, 64.5, NBins, FirstBin, LastBin);
	fH2["aF3EA"] = new TH2D("aF3EA","addback #gamma singles Three Frag", 64, 0.5, 64.5, NBins, FirstBin, LastBin);
	fH2["aF4EA"] = new TH2D("aF4EA","addback #gamma singles Four Frag", 64, 0.5, 64.5, NBins, FirstBin, LastBin);
	
	//fH2["Det11_EnVsCharge"] = new TH2D("Det11_EnVsCharge","GetArrayNumber()==11 energy versus charge", 2000, 0.5, 2000.5, 2000, 0.5, 2000.5);
/////////////////////////////////////////////////////////////////////////// Multiplicity Histograms
///////////////////////////////////////////////////////////////////////////  GRIFFIN ENERGY Spectra
   
   //Send histograms to Output list to be added and written.
   for(auto it : fH1) {
      GetOutputList()->Add(it.second);
	  // it->Sumw2();//If this works, it will propagate errors correctly
   }
   for(auto it : fH2) {
      GetOutputList()->Add(it.second);
	  // it.Sumw2();//If this works, it will propagate errors correctly
   }
   for(auto it : fHSparse) {
      GetOutputList()->Add(it.second);
	  // it.Sumw2();//If this works, it will propagate errors correctly
   }
}


void AddBackTestSelector::FillHistograms() {

///////////////////////////////////////   //GRIFFIN No addback with suppression



	fH1.at("gMult")->Fill( fGrif->GetMultiplicity() );
	fH1.at("aMult")->Fill( fGrif->GetAddbackMultiplicity() );
	
	
	for(auto g = 0; g < fGrif->GetMultiplicity(); g++){
		auto grif = fGrif->GetGriffinHit(g);
		if(PileUpRej) if(grif->GetKValue() != kValue_grif) continue;
		fH2.at("gEA")->Fill( grif->GetArrayNumber(), grif->GetEnergy() );
		
		//if( grif->GetArrayNumber() == 11 ) fH2.at("Det11_EnVsCharge")->Fill( grif->GetEnergy(), grif->GetCharge() );
	}

	
	//for(auto g = 0; g < fGrif->GetAddbackMultiplicity(); ++g) {
	for(auto g = 0; g < fGrif->GetSuppressedAddbackMultiplicity(fBgo); ++g) {
		auto grif = fGrif->GetSuppressedAddbackHit(g);
		//auto grif = fGrif->GetAddbackHit(g);
		fH1.at("aFrags")->Fill( fGrif->GetNAddbackFrags(g) );
		if(PileUpRej) if(grif->GetKValue() != kValue_grif) continue;
		fH2.at("aEA")->Fill( grif->GetArrayNumber(), grif->GetEnergy() );
		if( fGrif->GetNAddbackFrags(g) == 1 ) fH2.at("aF1EA")->Fill( grif->GetArrayNumber(), grif->GetEnergy() );
		if( fGrif->GetNAddbackFrags(g) == 2 ) fH2.at("aF2EA")->Fill( grif->GetArrayNumber(), grif->GetEnergy() );
		if( fGrif->GetNAddbackFrags(g) == 3 ) fH2.at("aF3EA")->Fill( grif->GetArrayNumber(), grif->GetEnergy() );
		if( fGrif->GetNAddbackFrags(g) == 4 ) fH2.at("aF4EA")->Fill( grif->GetArrayNumber(), grif->GetEnergy() );
	}
 }



