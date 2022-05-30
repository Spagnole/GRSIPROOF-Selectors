#define MyGGGCubeSelector_cxx
// The class definition in MyGGGCubeSelector.h has been generated automatically
#include "MyGGGCubeSelector.h"

int kValue_grif = 379; // default k-value for not piled-up hits

void MyGGGCubeSelector::CreateHistograms() {
	// define Trees
	fTree["ggg_tree"] = new TTree("ggg_tree", "tripple coincident events");
		// set branch addresses for output tree (these can be different for different trees)
	// we save the entry number as well to detect possible double counting
	// four coincident gammas a,b,c,d will be saved as a,b,c, a,b,d, a,c,d, and b,c,d
	fTree["ggg_tree"]->Branch("gE", &gE,"gE[3]/F");
	//fTree["ggg_tree"]->Branch("gTD", &gTD,"gTD[3]/F");
	fTree["ggg_tree"]->Branch("prompt_time", &prompt_time,"prompt_time[3]/O");
	//fTree["ggg_tree"]->Branch("gamma_mult", &gamma_mult);
	fTree["ggg_tree"]->Branch("CycleTime", &CycleTime, "CycleTime/F");
	//fTree["ggg_tree"]->Branch("BeamOff", &BeamOff);
	//We can also create histograms at the same time, maybe for some diagnostics or simple checks
	fH1["gE"] = new TH1D("gE","#total projection triple coincidence gamma-energy", 6000, 0.5, 6000.5); 
	fH1["ggTD"] = new TH1D("ggTD","#gamma-#gamma time differences", 10000, -5000, 5000); 



	//Send histograms to Output list to be added and written.
	for(auto it : fH1) {
		GetOutputList()->Add(it.second);
	}
	for(auto it : fH2) {
		GetOutputList()->Add(it.second);
	}
	for(auto it : fTree) {
		GetOutputList()->Add(it.second);
	}
}

double MyCycleLength = 783.5;

bool PromptCoincidence(TGriffinHit* g, TZeroDegreeHit* z){
	//Check if hits are less then 300 ns apart.
	return std::fabs(g->GetTime() - z->GetTime()) < 300.;
}

bool PromptCoincidence(TGriffinHit* g, TSceptarHit* s){
	//Check if hits are less then 300 ns apart.
	return std::fabs(g->GetTime() - s->GetTime()) < 300.;
}

bool PromptCoincidence(TGriffinHit* g1, TGriffinHit* g2){
	//Check if hits are less then 500 ns apart.
	//return std::fabs(g1->GetTime() - g2->GetTime()) < 300.;
	return -100. < (g2->GetTime() - g1->GetTime()) && (g2->GetTime() - g1->GetTime()) < 300.;
}

bool TimeRandom(TGriffinHit *h1, TGriffinHit *h2){  //Check if hits are less then 1000 ns apart.
   //return 500. < std::fabs(h1->GetTime() - h2->GetTime()) && std::fabs(h1->GetTime() - h2->GetTime()) < 1500.;
   return 500. < (h2->GetTime() - h1->GetTime()) && (h2->GetTime() - h1->GetTime()) < 1500.;
}

void MyGGGCubeSelector::FillHistograms() {
	// we could check multiplicities here and skip events where we do not have at least
	// three suppressed addback energies and a beta, but we want to fill some general 
	// histograms without these cuts.
	//fSuppressedAddback.resize(3, 0.);
	
	if( fGrif->GetSuppressedMultiplicity(fGriffinBgo) > 2 && fGrif->GetSuppressedMultiplicity(fGriffinBgo) < 8){	
		for(auto g1 = 0; g1 < fGrif->GetSuppressedMultiplicity(fGriffinBgo); ++g1) {
			auto grif1 = fGrif->GetSuppressedHit(g1);
			gE[0] = 0.; gE[1] = 0.; gE[2] = 0.;
			gE[0] = grif1->GetEnergy();
			gamma_mult = fGrif->GetSuppressedMultiplicity(fGriffinBgo);
			if(grif1->GetKValue() != kValue_grif) continue;
			for( auto g2 = g1+1; g2 < fGrif->GetSuppressedMultiplicity(fGriffinBgo); ++g2){
				if( g1 == g2 ) continue;
				auto grif2 = fGrif->GetSuppressedHit(g2);
				if(grif2->GetKValue() != kValue_grif) continue;
				gE[1] = grif2->GetEnergy();
				for( auto g3 = g2+1; g3 < fGrif->GetSuppressedMultiplicity(fGriffinBgo); ++g3){
					if( g1 == g3 || g2 == g3 ) continue;
					auto grif3 = fGrif->GetSuppressedHit(g3);
					if(grif3->GetKValue() != kValue_grif) continue;
					gE[2] = grif3->GetEnergy();
					prompt_time[0] = 0; prompt_time[1] = 0; prompt_time[2] = 0;
					CycleTime = fmod(grif1->GetTimeStampNs()/1E9, MyCycleLength);
					//gTD[0] = grif1->GetTime() - grif2->GetTime();
					//gTD[1] = grif1->GetTime() - grif3->GetTime();
					//gTD[2] = grif2->GetTime() - grif3->GetTime();
					
					
					if( PromptCoincidence(grif1,grif2) ) prompt_time[0] = 1;
						else if( TimeRandom(grif1,grif2) ) prompt_time[0] = 0;
							else continue;
					if( PromptCoincidence(grif1,grif3) ) prompt_time[1] = 1;
						else if( TimeRandom(grif1,grif3) ) prompt_time[1] = 0;
							else continue;
					if( PromptCoincidence(grif2,grif3) ) prompt_time[2] = 1;
						else if( TimeRandom(grif2,grif3) ) prompt_time[2] = 0;
							else continue;
					
					//CycleTime = fmod(grif1->GetTimeStampNs()/1E9, MyCycleLength)+0.5; 		
					if( gE[0] > 0. && gE[1] > 0. && gE[2] > 0.){
						fTree.at("ggg_tree")->Fill();
						fH1.at("gE")->Fill( gE[0] );	fH1.at("gE")->Fill( gE[1] );	fH1.at("gE")->Fill( gE[2] );
						fH1.at("ggTD")->Fill( grif2->GetTime() - grif1->GetTime() );
						fH1.at("ggTD")->Fill( grif3->GetTime() - grif1->GetTime() );
						fH1.at("ggTD")->Fill( grif3->GetTime() - grif2->GetTime() );
					}
				}
			
			}
		}
	}
}
