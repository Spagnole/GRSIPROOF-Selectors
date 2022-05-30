#define GammaSelector_cxx
// The class definition in GammaSelector.h has been generated automatically
#include "GammaSelector.h"

////Global variables for the time coincidence windows

double ggtime=120.;//Time coincidence windows for GRIFFIN-GRIFFIN, assumed symmetrical
double gstime_l=-210.;//Time coincidence windows for GRIFFIN-SCEPTAR, lower limit
double gstime_u=10.;//Time coincidence windows for GRIFFIN-SCEPTAR, upper limit
double gztime_l=-500.;//Time coincidence windows for GRIFFIN-ZDS, lower limit
double gztime_u=200.;//Time coincidence windows for GRIFFIN-ZDS, upper limit
double gptime_l=-1000.;//Time coincidence windows for GRIFFIN-PACES, lower limit
double gptime_u=300.;//Time coincidence windows for GRIFFIN-PACES, upper limit

double randtime_l=1000.;//Random time coincidence window, common for all coincidences, lower limit
double randtime_u=2000.;////Random time coincidence window, common for all coincidences, upper limit

// prompt coincidences
bool PromptCoincidence(TGriffinHit *h1, TGriffinHit *h2){  //Check if GRIFFIN hits are less then 120 ns apart.
   return std::fabs(h1->GetTime() - h2->GetTime()) < ggtime;
}
bool PromptCoincidence(TGriffinHit *h1, TSceptarHit *h2){  // coincidence window GRIFFIN-SCEPTAR: (-150,-20) ns
   //return -150. < h1->GetTimeStampNs() - h2->GetTimeStampNs() && h1->GetTimeStampNs() - h2->GetTimeStampNs() < -20.;
   return gstime_l < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < gstime_u;
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
bool TimeRandom(TGriffinHit *h1, TSceptarHit *h2){  //Check if hits are less then 1000 ns apart.
   //return 1000. < std::fabs(h1->GetTimeStampNs() - h2->GetTimeStampNs()) && std::fabs(h1->GetTimeStampNs() - h2->GetTimeStampNs()) < 2000.;
   return randtime_l < std::fabs(h1->GetTime() - h2->GetTime()) && std::fabs(h1->GetTime() - h2->GetTime()) < randtime_u;
}
bool TimeRandom(TGriffinHit *h1, TZeroDegreeHit *h2){  //Check if hits are less then 1000 ns apart.
   return randtime_l < std::fabs(h1->GetTime() - h2->GetTime()) && std::fabs(h1->GetTime() - h2->GetTime()) < randtime_u;
}
bool TimeRandom(TGriffinHit *h1, TPacesHit *h2){  //Check if hits are less then 1000 ns apart.
   return randtime_l < std::fabs(h1->GetTime() - h2->GetTime()) && std::fabs(h1->GetTime() - h2->GetTime()) < randtime_u;
}


void GammaSelector::CreateHistograms() {
   // get cycle length from ODB - this is in nanoseconds! (Actually it's not in ns, I got it in mu s!)
   fCycleLength = fPpg->OdbCycleLength();
   std::cout<<"got ODB cycle length "<<fCycleLength<<" mu s, "<<fCycleLength/1e6<<" s"<<std::endl;



/////////////////////////////////////////////////////////////////////////// Multiplicity Histograms
/*
  fH1["gMult"] = new TH1D("gMult","GRIFFIN Multiplicity",20, 0, 20);
  fH1["aMult"] = new TH1D("aMult","Addback Multiplicity",20, 0, 20);
  fH1["sMult"] = new TH1D("sMult","SCEPTAR Multiplicity",20, 0, 20);
  fH1["zMult"] = new TH1D("zMult","ZDS Multiplicity",20, 0, 20);
  fH1["pMult"] = new TH1D("pMult","PACES Multiplicity",20, 0, 20);
  fH1["lMult"] = new TH1D("lMult","LaBr3(Ce) Multiplicity",20, 0, 20);
  fH1["tMult"] = new TH1D("tMult","TAC Multiplicity",20, 0, 20);
*/
///////////////////////////////////////////////////////////////////////////  GRIFFIN ENERGY Spectra

   fH2["gEA"] = new TH2D("gEA","#gamma singles vs. crystal #", 64, 0., 64., 16384, 0., 8192.); //griffin energy (all channels)
   fH1["gE"] = new TH1D("gE","#gamma singles", 16384, 0., 8192.);
   fH1["gEnonSuppr"] = new TH1D("gEnonSuppr","#gamma singles, non-suppressed", 16384, 0., 8192.);
   fH2["gEC"] = new TH2D("gEC","#gamma-energy vs cycle time", fCycleLength/1e4, 0, fCycleLength/1e6, 4000, 0, 4000); //griffin energy (all channels)
   fH2["gzEC"] = new TH2D("gzEC","ZDS-tagged #gamma-energy vs cycle time", fCycleLength/1e4, 0, fCycleLength/1e6, 4000, 0, 2000); //griffin energy (all channels)
   fH2["ggE"] = new TH2D("ggE","#gamma-#gamma matrix, time-random subtracted", 4096, 0, 4096, 4096, 0, 4096); //griffin energy (all channels)
   fH2["ggEBg"] = new TH2D("ggEBg","#gamma-#gamma matrix, time random BG", 4096, 0, 4096, 4096, 0, 4096); //griffin energy (all channels)
   fH1["ggT"] = new TH1D("ggT", "GRIFFIN-GRIFFIN timing", 4000, -2000., 2000.);
   fH2["gT"] = new TH2D("gT", "#gamma energy vs. time", 3600, 0., 3600., 4000, 0., 4000.);

   // Cycle-gated gamma-gamma matrices
   fH2["ggzEC_bg"] = new TH2D("ggzEC_bg", "#gamma-#gamma matrix, background tape cycle", 4096, 0, 4096, 4096, 0, 4096);
   fH2["ggzEC_beam_on"] = new TH2D("ggzEC_beam_on", "#gamma-#gamma matrix, beam on tape cycle", 4096, 0, 4096, 4096, 0, 4096);
   fH2["ggzEC_beam_off"] = new TH2D("ggzEC_beam_off", "#gamma-#gamma matrix, beam off tape cycle", 4096, 0, 4096, 4096, 0, 4096);

/*
   fH2["aEA"] = new TH2D("aEA","Addback #gamma vs crystal #", 64, 0., 64., 16384, 0., 8192.);
   fH1["aE"] = new TH1D("aE","Addback #gamma", 16384, 0., 8192.);
   fH2["aEC"] = new TH2D("aEC","Addback #gamma-energy vs cycle time", fCycleLength/1e4, 0, fCycleLength/1e6, 8192, 0, 8192); //griffin energy (all channels)
   fH2["aaE"] = new TH2D("aaE","Addback #gamma-#gamma matrix", 8192, 0., 8192., 8192, 0., 8192.); //griff$
   fH2["aaEBg"] = new TH2D("aaEBg","Addback #gamma-#gamma matrix, time random BG", 8192, 0, 8192, 8192, 0, 8192);
   fH1["aaT"] = new TH1D("aaT", "Addback-Addback timing", 4000, -2000., 2000.);
*/
///////////////////////////////////////////////////////////////////////////  Beta ENERGY Spectra

//   fH2["sE"] = new TH2D("sE","SCEPTAR Energy vs. Detector #", 20 , 0, 20, 4096, 0, 16384);
//   fH2["sEC"] = new TH2D("sEC","SCEPTAR Energy vs cycle time", fCycleLength/1e4, 0, fCycleLength/1e6, 4096, 0, 16384);
   fH1["zE"] = new TH1D("zE","ZDS Energy", 4096, 0, 4096);
   fH2["zEC"] = new TH2D("zEC","ZDS Energy vs cycle time", fCycleLength/1e4, 0, fCycleLength/1e6, 4096, 0, 16384);

///////////////////////////////////////////////////////////////////////////  Beta Tagged GRIFFIN Timing Spectra

//   fH1["gsT"] = new TH1D("gsT", "Griffin-SCEPTAR timing", 4000, -2000, 2000);
//   fH2["gsTE"] = new TH2D("gsTE", "Griffin-SCEPTAR timing vs. SCEPTAR Energy", 4000, -2000, 2000, 4096,0, 16384);
//   fH2["gsTD"] = new TH2D("gsTD", "Griffin-SCEPTAR timing vs. SCEPTAR Detector #", 20, 0, 20, 4000, -2000, 2000);
   fH1["gzT"] = new TH1D("gzT", "Griffin-ZDS timing", 4000, -2000, 2000);
   fH2["gzTE"] = new TH2D("gzTE", "Griffin-ZDS timing vs. ZDS Energy", 4000, -2000, 2000, 4096,0, 16384);

   fH2["gzTD2"] = new TH2D("gzTD2", "Griffin-ZDS timing vs. HPGe Detector #", 64, 0, 64, 4000, -2000, 2000);
//   fH2["gsTD2"] = new TH2D("gsTD2", "Griffin-SCEPTAR timing vs. HPGe Detector #", 64., 0, 64., 4000, -2000, 2000);
/*
   fH1["asT"] = new TH1D("asT", "Griffin-Addback-SCEPTAR timing", 4000, -2000, 2000);
   fH2["asTE"] = new TH2D("asTE", "Griffin-Addback-SCEPTAR timing vs. SCEPTAR Energy", 4000, -2000, 2000, 4096,0, 16384);
   fH2["asTD"] = new TH2D("asTD", "Griffin-Addback-SCEPTAR timing vs. SCEPTAR Detector #", 20, 0, 20, 4000, -2000, 2000);
   fH1["azT"] = new TH1D("azT", "Griffin-Addback-ZDS timing", 4000, -2000, 2000);
   fH2["azTE"] = new TH2D("azTE", "Griffin-Addback-ZDS timing vs. ZDS Energy", 4000, -2000, 2000, 4096,0, 16384);
*/
///////////////////////////////////////////////////////////////////////////  Beta Tagged GRIFFIN Energy Spectra

//   fH1["gsE"] = new TH1D("gsE","#gamma singles, SCEPTAR Prompt Time Gated", 16384, 0., 8192.);
//   fH2["gsEA"] = new TH2D("gsEA","#gamma singles vs. crystal #, SCEPTAR Prompt Time Gated", 64, 0., 64., 16384, 0., 8192.);
   fH1["gzE"] = new TH1D("gzE","#gamma singles, ZDS Prompt Time Gated", 16384, 0., 8192.);
   fH2["gzEA"] = new TH2D("gzEA","#gamma singles vs. crystal #, ZDS Prompt Time Gated", 64, 0., 64., 16384, 0., 8192.);
//   fH1["gbE"] = new TH1D("gbE","#gamma singles, any Beta Prompt Time Gated", 16384, 0., 8192.);
//   fH2["gbEA"] = new TH2D("gbEA","#gamma singles vs. crystal #, any Beta Prompt Time Gated", 64, 0., 64., 16384, 0., 8192.);
/*
   fH1["gsE_Bg"] = new TH1D("gsE_Bg","#gamma singles, SCEPTAR Random Time Gated", 16384, 0., 8192.);
   fH2["gsEA_Bg"] = new TH2D("gsEA_Bg","#gamma singles vs. crystal #, SCEPTAR Random Time Gated", 64, 0., 64., 16384, 0., 8192.);
   fH1["gzE_Bg"] = new TH1D("gzE_Bg","#gamma singles, ZDS Random Time Gated", 16384, 0., 8192.);
   fH2["gzEA_Bg"] = new TH2D("gzEA_Bg","#gamma singles vs. crystal #, ZDS Random Time Gated", 64, 0., 64., 16384, 0., 8192.);
   fH1["gbE_Bg"] = new TH1D("gbE_Bg","#gamma singles, any #beta Random Time Gated", 16384, 0., 8192.);
   fH2["gbEA_Bg"] = new TH2D("gbEA_Bg","#gamma singles vs. crystal #, any #beta Random Time Gated", 64, 0., 64., 16384, 0., 8192.);

   fH1["asE"] = new TH1D("asE","Addback #gamma, SCEPTAR Prompt Time Gated", 16384, 0., 8192.);
   fH2["asEA"] = new TH2D("asEA","Addback #gamma vs. crystal #, SCEPTAR Prompt Time Gated", 64, 0., 64., 16384, 0., 8192.);
   fH1["azE"] = new TH1D("azE","#Addback gamma, ZDS Prompt Time Gated", 16384, 0., 8192.);
   fH2["azEA"] = new TH2D("azEA","Addback #gamma vs. crystal #, ZDS Prompt Time Gated", 64, 0., 64., 16384, 0., 8192.);
   fH1["abE"] = new TH1D("abE","Addback #gamma, any Beta Prompt Time Gated", 16384, 0., 8192.);
   fH2["abEA"] = new TH2D("abEA","Addback #gamma vs. crystal #, any Beta Prompt Time Gated", 64, 0., 64., 16384, 0., 8192.);

   fH1["asE_Bg"] = new TH1D("asE_Bg","Addback #gamma, SCEPTAR Random Time Gated", 16384, 0., 8192.);
   fH2["asEA_Bg"] = new TH2D("asEA_Bg","Addback #gamma vs. crystal #, SCEPTAR Random Time Gated", 64, 0., 64., 16384, 0., 8192.);
   fH1["azE_Bg"] = new TH1D("azE_Bg","Addback #gamma, ZDS Random Time Gated", 16384, 0., 8192.);
   fH2["azEA_Bg"] = new TH2D("azEA_Bg","Addback #gamma vs. crystal #, ZDS Random Time Gated", 64, 0., 64., 16384, 0., 8192.);
   fH1["abE_Bg"] = new TH1D("abE_Bg","Addback #gamma, any Beta Random Time Gated", 16384, 0., 8192.);
   fH2["abEA_Bg"] = new TH2D("abEA_Bg","Addback #gamma vs. crystal #, any Beta Random Time Gated", 64, 0., 64., 16384, 0., 8192.);
*/
///////////////////////////////////////////////////////////////////////////  Beta Tagged GRIFFIN-GRIFFIN Energy Spectra

//  fH2["ggsE"] = new TH2D("ggsE","#gamma-#gamma matrix, SCEPTAR Prompt Time Gated", 8192, 0, 8192, 8192, 0, 8192);
//   fH2["ggsEBg"] = new TH2D("ggsEBg","#gamma-#gamma matrix, SCEPTAR time random BG", 8192, 0, 8192, 8192, 0, 8192);
   fH2["ggzE"] = new TH2D("ggzE","#gamma-#gamma matrix, ZDS Prompt Time Gated, time-random subtracted", 4096, 0, 4096, 4096, 0, 4096);
   fH2["ggzEBg"] = new TH2D("ggzEBg","#gamma-#gamma matrix, ZDS time random BG", 4096, 0, 4096, 4096, 0, 4096);
//   fH2["ggbE"] = new TH2D("ggbE","#gamma-#gamma matrix, any #beta Prompt Time Gated", 8192, 0, 8192, 8192, 0, 8192);
//   fH2["ggbEBg"] = new TH2D("ggbEBg","#gamma-#gamma matrix, and #beta time random BG", 8192, 0, 8192, 8192, 0, 8192);
/*
   fH2["aasE"] = new TH2D("aasE","Addback #gamma-#gamma matrix, SCEPTAR Prompt Time Gated", 8192, 0, 8192, 8192, 0, 8192);
   fH2["aasEBg"] = new TH2D("aasEBg","Addback #gamma-#gamma matrix, SCEPTAR time random BG", 8192, 0, 8192, 8192, 0, 8192);
   fH2["aazE"] = new TH2D("aazE","Addback #gamma-#gamma matrix, ZDS Prompt Time Gated", 8192, 0, 8192, 8192, 0, 8192);
   fH2["aazEBg"] = new TH2D("aazEBg","Addback #gamma-#gamma matrix, ZDS time random BG", 8192, 0, 8192, 8192, 0, 8192);
   fH2["aabE"] = new TH2D("aabE","Addback #gamma-#gamma matrix, any #beta Prompt Time Gated", 8192, 0, 8192, 8192, 0, 8192);
   fH2["aabEBg"] = new TH2D("aabEBg","Addback #gamma-#gamma matrix, and #beta time random BG", 8192, 0, 8192, 8192, 0, 8192);
*/
   ///////////////////////////////////////////////////////////////////////////  PACES ENERGY Spectra

   fH1["pE"] = new TH1D("pE","PACES Energy", 8192, 0, 8192);
   fH2["pEA"] = new TH2D("pEA","PACES Energy vs. Detector #", 6 , 0, 6, 4096, 0, 4096);
   fH2["pEC"] = new TH2D("pEC","PACES Energy vs cycle time", fCycleLength/1e4, 0, fCycleLength/1e6, 4096, 0, 4096);
   
///////////////////////////////////////////////////////////////////////////  PACES-GRIFFIN Energy Spectra

   fH1["gpT"] = new TH1D("gpT", "PACES-Griffin timing", 4000, -2000, 2000);
   fH2["gpE"] = new TH2D("gpE","PACES-#gamma matrix, time-random subtracted", 4096, 0, 4096, 4096, 0, 4096); 
   //Daniel changed this 13/7/21 when troubleshooting the sort code. GRSIProof keeps indicating a memory leak associated with gpE
   //fH2["gpE"] = new TH2D("gpE","PACES-#gamma matrix, time-random subtracted", 10,0,10,10,0,10);
   fH2["gpEBg"] = new TH2D("gpE","PACES-#gamma matrix", 4096, 0, 4096, 4096, 0, 4096);
      
   for(int i = 1; i < 6; ++i) {//Lazy loop to create histograms for individual PACES crystals
      fH2[Form("gpE_%d",i)] = new TH2D(Form("gpE_%d",i),Form("#gamma-CE matrix;E_{#gamma} [keV];E_{#gamma} [keV] crystal %d, time-random subtracted",i), 4000, 0., 4000., 2000, 0., 2000.);
      fH2[Form("gpEBG_%d",i)] = new TH2D(Form("gpEBG_%d",i),Form("#gamma-CE matrix BG;E_{#gamma} [keV];E_{#gamma} [keV] crystal %d",i), 4000, 0., 4000., 2000, 0., 2000.);
   }
   

///////////////////////////////////////////////////////////////////////////  LaBr3(Ce) + TAC ENERGY Spectra

   fH1["lE"] = new TH1D("lE","LaBr3(Ce) Energy", 4096, 0, 4096);
   fH1["lEnonSuppr"] = new TH1D("lEnonSuppr","LaBr3(Ce) Energy no suppressed", 4096, 0, 4096);
   fH2["lED"] = new TH2D("lED","LaBr3(Ce) Energy vs. Detector #", 8, 0, 8, 4096, 0, 4096);
   fH2["lEC"] = new TH2D("lEC","LaBr3(Ce) Energy vs cycle time", fCycleLength/1e4, 0, fCycleLength/1e6, 4096, 0, 4096);
   fH1["lzE"] = new TH1D("lzE","LaBr3(Ce) tagged with ZDS and its TAC;E_{#gamma} [keV];Counts/keV", 4096, 0., 4096); //LaBr singles energy
   fH2["lltE"] = new TH2D("lltE","LaBr3(Ce) #gamma-#gamma matrix with corresponding TAC tag ;START E_{#gamma} [keV];STOP E_{#gamma} [keV]", 2000, 0., 2000., 2000, 0., 2000.); //LaBr-LaBr with the right TAC coincidence
   fH1["lTAC"] = new TH1D("lTAC","LaBr TAC Time [ps];Counts/10 ps", 5000, 0, 50000); //Not very useful, they need to be gainmatched
   fH1["zTAC"] = new TH1D("zTAC","ZDS TAC Time [ps];Counts/10 ps", 5000, 0, 50000); 
//   fH2["tE"] = new TH2D("tE","TAC Time [ps];Counts/10 ps vs Detector #", 8, 0,8, 5000, 0, 50000); //Not very useful, they need to be gainmatched
//   fH1["tzE"] = new TH1D("tzE","ZDS TAC Time [ps];Counts/10 ps", 5000, 0, 50000); //Not very useful, they need to be gainmatched

   
   
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


void GammaSelector::FillHistograms() {

///////////////////////////////////////   //GRIFFIN No addback with suppression

   for(auto g = 0; g < fGrif->GetSuppressedMultiplicity(fBgo); ++g) {
      auto grif = fGrif->GetSuppressedHit(g);
      fH1.at("gE")->Fill(grif->GetEnergy());
      fH2.at("gEA")->Fill(grif->GetArrayNumber(), grif->GetEnergy()); 

      fH2.at("gT")->Fill(grif->GetTime()/1e9, grif->GetEnergy());

      fH2.at("gEC")->Fill(fmod(grif->GetTimeStampNs()/1e3, fCycleLength)/1e6, grif->GetEnergy());

      //GRIFFIN-GRIFFIN
      for(auto g2 = 0; g2 < fGrif->GetSuppressedMultiplicity(fBgo); ++g2) {

	if(g==g2) continue;

         auto grif2 = fGrif->GetSuppressedHit(g2);

          fH1.at("ggT")->Fill(grif->GetTime()-grif2->GetTime());

         if(PromptCoincidence(grif, grif2)) {
            fH2.at("ggE")->Fill(grif->GetEnergy(), grif2->GetEnergy());
         } else if(TimeRandom(grif, grif2)) {
            fH2.at("ggEBg")->Fill(grif->GetEnergy(), grif2->GetEnergy(), 2*ggtime/(2*(randtime_u-randtime_l)));//The factors 2* come from using abs() values
            fH2.at("ggE")->Fill(grif->GetEnergy(), grif2->GetEnergy(), -2*ggtime/(2*(randtime_u-randtime_l)));//To subtract time random from the true coincidences
         }
      }

      //GRIFFIN-ZDS
      for(auto z = 0; z < fZds->GetMultiplicity(); ++z) {
         auto zds = fZds->GetZeroDegreeHit(z);

       fH1.at("gzT")->Fill(grif->GetTime() - zds->GetTime());
       fH2.at("gzTE")->Fill(grif->GetTime() - zds->GetTime(), zds->GetEnergy());   

       fH2.at("gzTD2")->Fill(grif->GetArrayNumber(), grif->GetTime() - zds->GetTime());   

         if(PromptCoincidence(grif, zds)) {

         fH1.at("gzE")->Fill(grif->GetEnergy());

         fH2.at("gzEC")->Fill(fmod(grif->GetTimeStampNs()/1e3, fCycleLength)/1e6, grif->GetEnergy());

//       fH1.at("gbE")->Fill(grif->GetEnergy());

            //GRIFFIN-GRIFFIN-ZDS
            for(auto g2 = 0; g2 < fGrif->GetSuppressedMultiplicity(fBgo); ++g2) {

		if(g==g2) continue;

               auto grif2 = fGrif->GetSuppressedHit(g2);
               if(PromptCoincidence(grif, grif2)) {
                  fH2.at("ggzE")->Fill(grif->GetEnergy(), grif2->GetEnergy());
                  
                  // bg cycle time=(1.50, 1.70) 200ms
                  if ((fmod(grif->GetTimeStampNs()/1e3, fCycleLength)/1e6 > 1.50) && (fmod(grif->GetTimeStampNs()/1e3, fCycleLength)/1e6 < 1.70)) {
                      fH2.at("ggzEC_bg")->Fill(grif->GetEnergy(), grif2->GetEnergy());
		  }
                  // beam on cycle time=(1.70, 7.70) 6sec
                  if ((fmod(grif->GetTimeStampNs()/1e3, fCycleLength)/1e6 > 1.70) && (fmod(grif->GetTimeStampNs()/1e3, fCycleLength)/1e6 < 7.70)) {
                      fH2.at("ggzEC_beam_on")->Fill(grif->GetEnergy(), grif2->GetEnergy());
		  }
                  // beam off cycle time=(7.70, 8.70) 1sec
                  if ((fmod(grif->GetTimeStampNs()/1e3, fCycleLength)/1e6 > 7.70) && (fmod(grif->GetTimeStampNs()/1e3, fCycleLength)/1e6 < 8.70)) {
                      fH2.at("ggzEC_beam_off")->Fill(grif->GetEnergy(), grif2->GetEnergy());
		  }

//                  fH2.at("ggbE")->Fill(grif->GetEnergy(), grif2->GetEnergy());
               } 
	       else if(TimeRandom(grif, grif2)) {
                  fH2.at("ggzEBg")->Fill(grif->GetEnergy(), grif2->GetEnergy(), (gztime_u-gztime_l)/(2*(randtime_u-randtime_l)));
                  fH2.at("ggzE")->Fill(grif->GetEnergy(), grif2->GetEnergy(), -(gztime_u-gztime_l)/(2*(randtime_u-randtime_l)));

//                  fH2.at("ggbEBg")->Fill(grif->GetEnergy(), grif2->GetEnergy());
               }
            }
         }
         else if(TimeRandom(grif, zds)) {

         fH1.at("gzE")->Fill(grif->GetEnergy(), -(gztime_u-gztime_l)/(2*(randtime_u-randtime_l)));
         fH2.at("gzEC")->Fill(fmod(grif->GetTimeStampNs()/1e3, fCycleLength)/1e6, grif->GetEnergy(), -(gztime_u-gztime_l)/(2*(randtime_u-randtime_l)));
         }
         
      }

      //GRIFFIN-SCEPTAR
	/*  
      for(auto s = 0; s < fScep->GetMultiplicity(); ++s) {
         auto scep = fScep->GetSceptarHit(s);

       fH2.at("sE")->Fill(scep->GetDetector(), scep->GetEnergy());
       fH1.at("gsT")->Fill(grif->GetTime() - scep->GetTime());
       fH2.at("gsTE")->Fill(grif->GetTime() - scep->GetTime(), scep->GetEnergy());      
       fH2.at("gsTD")->Fill(grif->GetTime() - scep->GetTime(), scep->GetDetector());

       fH2.at("gsTD2")->Fill(grif->GetArrayNumber(), grif->GetTime() - scep->GetTime()); 

         if(PromptCoincidence(grif, scep)) {

       fH1.at("gsE")->Fill(grif->GetEnergy());

       fH1.at("gbE")->Fill(grif->GetEnergy());

            //GRIFFIN-GRIFFIN-SCEPTAR
            for(auto g2 = 0; g2 < fGrif->GetSuppressedMultiplicity(fBgo); ++g2) {

		if(g==g2) continue;

               auto grif2 = fGrif->GetSuppressedHit(g2);
               if(PromptCoincidence(grif, grif2)) {
                  fH2.at("ggsE")->Fill(grif->GetEnergy(), grif2->GetEnergy());

                  fH2.at("ggbE")->Fill(grif->GetEnergy(), grif2->GetEnergy());
               } 
	       else if(TimeRandom(grif, grif2)) {
                  fH2.at("ggsEBg")->Fill(grif->GetEnergy(), grif2->GetEnergy());

                  fH2.at("ggbEBg")->Fill(grif->GetEnergy(), grif2->GetEnergy());
               }
            }
         }
      }
	*/
  
       //GRIFFIN-PACES
	for(auto p = 0; p < fPaces->GetMultiplicity(); ++p) {
         auto paces = fPaces->GetPacesHit(p);
		 
       fH1.at("gpT")->Fill(grif->GetTime() - paces->GetTime());

         if(PromptCoincidence(grif, paces)) {

			fH2.at("gpE")->Fill(grif->GetEnergy(), paces->GetEnergy());
			fH2.at(Form("gpE_%d",paces->GetDetector()))->Fill(grif->GetEnergy(), paces->GetEnergy());
        
         }
		 
         else if(TimeRandom(grif, paces)) {

			fH2.at("gpEBg")->Fill(grif->GetEnergy(), paces->GetEnergy(), (gptime_u-gptime_l)/(2*(randtime_u-randtime_l)));
			fH2.at(Form("gpEBG_%d",paces->GetDetector()))->Fill(grif->GetEnergy(), paces->GetEnergy(), (gptime_u-gptime_l)/(2*(randtime_u-randtime_l)));
			fH2.at("gpE")->Fill(grif->GetEnergy(), paces->GetEnergy(), -(gptime_u-gptime_l)/(2*(randtime_u-randtime_l)));
			fH2.at(Form("gpE_%d",paces->GetDetector()))->Fill(grif->GetEnergy(), paces->GetEnergy(), -(gptime_u-gptime_l)/(2*(randtime_u-randtime_l)));
        
         }
      }
  }//First GRIFFIN loop

////////////////////////////////// GRIFFIN no addback, no suppression
   for(auto g = 0; g < fGrif->GetMultiplicity(); ++g) {
      auto grif = fGrif->GetGriffinHit(g);
      fH1.at("gEnonSuppr")->Fill(grif->GetEnergy());
   }
   
   
////////////////////////////////// ZDS singles

   for(auto z = 0; z < fZds->GetMultiplicity(); ++z) {
      auto zds = fZds->GetZeroDegreeHit(z);
      fH1.at("zE")->Fill(zds->GetEnergy());
      fH2.at("zEC")->Fill(fmod(zds->GetTimeStampNs()/1e3, fCycleLength)/1e6, zds->GetEnergy());
   }


////////////////////////////////// PACES singles
   for(auto p = 0; p < fPaces->GetMultiplicity(); ++p) {
      auto paces = fPaces->GetPacesHit(p);
      fH1.at("pE")->Fill(paces->GetEnergy());
      fH2.at("pEA")->Fill(paces->GetDetector(), paces->GetEnergy()); 
      fH2.at("pEC")->Fill(fmod(paces->GetTimeStampNs()/1e3, fCycleLength)/1e6, paces->GetEnergy());
   }
////////////////////////////////// LaBr suppressed
  for(auto l = 0; l < fLaBr->GetSuppressedMultiplicity(fLaBrBgo); l++){
     auto labr = fLaBr->GetSuppressedHit(l);
     fH1.at("lE")->Fill(labr->GetEnergy());
     fH2.at("lED")->Fill(labr->GetEnergy(), labr->GetDetector());
	 
	 if(fTAC->GetMultiplicity()==1&&fZds->GetMultiplicity()==1&&fLaBr->GetMultiplicity()==1&&fTAC->GetHit(0)->GetDetector()==8){
		fH1.at("lzE")->Fill(labr->GetEnergy());
	 }
  }
   
/////////////////////////////////// LaBr x2 suppressed with TAC coincidence

    if(fTAC->GetMultiplicity()==1&&fLaBr->GetSuppressedMultiplicity(fLaBrBgo)==2) {	
        auto labr1 = fLaBr->GetSuppressedHit(0);		
        auto labr2 = fLaBr->GetSuppressedHit(1);
		
	auto l1 = labr1->GetDetector();
	auto l2 = labr2->GetDetector();
	auto tac = fTAC->GetHit(0)->GetDetector();
	
	if(tac<8) {fH1.at("lTAC")->Fill(fTAC->GetHit(0)->GetEnergy());}
	else if(tac==8) {fH1.at("zTAC")->Fill(fTAC->GetHit(0)->GetEnergy());}

	if(l1<l2&&l1==tac) {
		fH2.at("lltE")->Fill(labr1->GetEnergy(), labr2->GetEnergy());
	}else if(l2<l1&&l2==tac) {
		fH2.at("lltE")->Fill(labr2->GetEnergy(), labr1->GetEnergy());			
	}
    }

   
////////////////////////////////// LaBr no suppression
  for(auto l = 0; l < fLaBr->GetMultiplicity(); l++){
     auto labr = fLaBr->GetLaBrHit(l);
     fH1.at("lEnonSuppr")->Fill(labr->GetEnergy());
   }


}


