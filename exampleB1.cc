//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: exampleB1.cc 86065 2014-11-07 08:51:15Z gcosmo $
//
/// \file exampleB1.cc
/// \brief Main program of the B1 example

#include "B1DetectorConstruction.hh"
#include "B1DetectorConstructionNEMA.hh"
#include "B1DetectorConstructionSTEAK.hh"
#include "B1ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
//#include "QBBC.hh"
#include "B1PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"
#include "G4StepLimiter.hh"
#include "G4UserLimits.hh"
#include "G4StepLimiterPhysics.hh"

#include <stdio.h>
#include <stdlib.h>
#include "SteppingVerbose.hh"

#include "G4GlobalConfig.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
using namespace std;


int main(int argc,char** argv)
{
	
	G4bool VisFlag=false;
	G4bool QuickFlag=false;

	// Detect interactive mode (if no arguments) and define UI session
	G4UIExecutive* ui = 0;
	
	G4double sphereDistY=0, CenterSphere=0, DetConf=0*mm, AbsorberThickness=1*mm, TBRvalue=1;
	G4int FilterFlag=1, SourceChoice=1, IsotopeChoice=1, AbsorberMaterial=1, QuickFlagCommandLine=0, SphereSelect=6, nThreadIn=-1;
	
	G4String MacroName ="";
	G4String FileNameLabel="";

	G4double DTmis=600; //in s
	G4double AttSorg[10]={ //in Bq
		350, //PSr
		2.53e3, //ExtSr Rm
		26.5e3, //Y - according to rosa's data for Z scan (9/1/18)
		23.24e3, //Co60
		16.04e3, //Na22
		33.89e3, //Ba133
		37.84e3, //Cs137
		99,
		99,
		27.110e3, //Na22 nuda PG
	};
	
	G4bool NoOfPrimToGenChangeFlag=false;
	/*
	 2.53e3 - ExtSrRome
	 23.24e3 - Co60
	 16.04e3 - Na22
	 33.89e3 - Ba133
	 37.84e3 - Cs137
	 1 - Pointlike Sr
	 2 - Extended Sr
	 3 - ExtY
	 4 - Co60
	 5 - Na22
	 6 - Ba133
	 7 - Cs137
	 8 - FlatEle 0-3MeV
	 9 - FlatGamma 0-1 MeV
	 10 - Na22 nuda
	 */
	
	G4double PixelThickness=2.5; //um
	G4int NoOfPrimToGen=99, Verbose=0;

	for(int i=1;i<argc;i++)
		if(argv[i][0] =='-')
		{
			G4String option(argv[i]);
			G4cout<<"option: "<<i<<" "<<option<<G4endl;
			if(option.compare("-DetConf")==0)
			{
				DetConf=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-AbsT")==0)
			{
				AbsorberThickness=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-AbsMat")==0)
			{
				AbsorberMaterial=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-Center")==0)
			{
				CenterSphere=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-Sphere")==0)
			{
				SphereSelect=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-NThreads")==0)
			{
				nThreadIn=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-TBR")==0)
			{
				TBRvalue=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-Source")==0)
			{
				SourceChoice=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-SphereDistY")==0)
			{
				sphereDistY=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-Vis")==0)
			{
				VisFlag=stoi (argv[++i], NULL);;
			}
			else if(option.compare("-Verbose")==0)
			{
				Verbose=stoi (argv[++i], NULL);;
			}
			else if(option.compare("-Isotope")==0)
			{
				IsotopeChoice=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-PixT")==0)
			{
				PixelThickness=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-NPrim")==0)
			{
				NoOfPrimToGen=stoi (argv[++i], NULL);;
			}
			else if(option.compare("-Quick")==0)
			{
				QuickFlagCommandLine=strtod (argv[++i], NULL);;
			}else if(option.compare("-Label")==0)
			{
				FileNameLabel= argv[++i];;
			}
			
		}
		else
		{
			MacroName = argv[i]; //se ho trovato una macro (senza il "-" davanti) significa che NON voglio l'interattivo
			VisFlag=false;
			QuickFlag=false;
		}
	
	if (VisFlag) QuickFlag=true;
	
	if (NoOfPrimToGen!=99) NoOfPrimToGenChangeFlag=true;
	
//	if (SourceChoice==3) DTmis=400;
//	if (SourceChoice==10) DTmis=1000;
//
//	if (!NoOfPrimToGenChangeFlag && (SourceChoice==2 || SourceChoice==3 || (SourceChoice>=4 && SourceChoice<=7) || SourceChoice==10)) { // If still 99 it means I did not choose a precise value via command line, so let's compute it! -   To be fixed: what to do in case of PSr/Y
//		NoOfPrimToGen=DTmis*AttSorg[(int)SourceChoice-1];
//	}
	G4cout<<"\n############## \nI WILL GENERATE n= "<<NoOfPrimToGen<<" primaries \n##############"<<G4endl;
	if (QuickFlagCommandLine) QuickFlag=true;

	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4long seed = time(NULL);
	if (VisFlag) seed=12345; //If vis was requested same always the same seed to have reproducibility
	G4Random::setTheSeed(seed);

//	G4int SphereSelect=SourceChoice;
	
	G4String MaterialiAssorbitore[3]= {"Cu","Al","ABS"};
	
//	G4String FileNamePrim="Primaries";
	G4String OutFileName="NEMA";
	G4String FileNameCommonPart;
//
//	FileNameCommonPart.append("_X"+ std::to_string((G4int)sphereDistY));
//	FileNameCommonPart.append("_Z"+ std::to_string((G4int)(100*CenterSphere)));
//
//	if (DetConf>=0) FileNameCommonPart.append("_AbsDz" + std::to_string((G4int)(1000*AbsorberThickness))+"_AbsHole" + std::to_string((G4int)(100*DetConf)) +"_AbsMat" + MaterialiAssorbitore[AbsorberMaterial-1]);
//	else FileNameCommonPart.append("_NoAbs");
//
	//	FileNameCommonPart.append("_Fil" + std::to_string((G4int)FilterFlag));
	
	if (SphereSelect>0) {//phantom NEMA
		FileNameCommonPart.append("_Nema_s" + std::to_string((G4int)SphereSelect));
	}else{//bistecca
		if (SphereSelect<-6)
		FileNameCommonPart.append("_Steak_Triangle");
		else 		FileNameCommonPart.append("_Steak_s"+ std::to_string(-(G4int)SphereSelect));

	}
	
	if (sphereDistY!=0) FileNameCommonPart.append("_Dist"+ std::to_string((G4int)sphereDistY));
	if (DetConf!=0) FileNameCommonPart.append("_Det"+ std::to_string((G4int)DetConf));

	if (SourceChoice==2) {//fondo
		FileNameCommonPart.append("_F");
	}else{//segnale
		FileNameCommonPart.append("_S");
	}
		//
	if (IsotopeChoice==1) FileNameCommonPart.append("_Ho");
	if (IsotopeChoice==2) FileNameCommonPart.append("_Tc");
	if (IsotopeChoice==3) FileNameCommonPart.append("_I");
	if (IsotopeChoice==4) FileNameCommonPart.append("_Cu");
	if (IsotopeChoice==5) FileNameCommonPart.append("_HoS");


	if (VisFlag) FileNameCommonPart.append("TEST"); //if it was a TEST run under vis
	if (QuickFlagCommandLine) FileNameCommonPart.append("_Quick"); //if it was a TEST run under vis

	if (FileNameLabel!="") FileNameCommonPart.append("_" + FileNameLabel);
	
	if (NoOfPrimToGenChangeFlag) FileNameCommonPart.append("_N"+std::to_string((G4int)NoOfPrimToGen)); //if it was a TEST run under vis
	
//	FileNamePrim.append(FileNameCommonPart+".dat");
	OutFileName.append(FileNameCommonPart);
	
//	std::ofstream primFile(FileNamePrim, std::ios::out);
	

	// Construct the default run manager
	//
#ifdef G4MULTITHREADED
	  G4MTRunManager* runManager = new G4MTRunManager;
//	runManager->SetNumberOfThreads( G4Threading::G4GetNumberOfCores() );
	G4int numOfThreads=(Verbose>0||VisFlag)?1:nThreadIn!=-1? nThreadIn:G4Threading::G4GetNumberOfCores() -2;
	runManager->SetNumberOfThreads(numOfThreads);
#else
	G4RunManager* runManager = new G4RunManager;
#endif
	////	G4VSteppingVerbose::SetInstance(new SteppingVerbose); //to use my SteppingVerbose

	// Set mandatory initialization classes
	// Detector construction
	if (SphereSelect>0) { // geometria NEMA-like
	runManager->SetUserInitialization(new B1DetectorConstructionNEMA(sphereDistY, CenterSphere, DetConf, AbsorberThickness, AbsorberMaterial, FilterFlag, SphereSelect, IsotopeChoice, QuickFlag, PixelThickness));
	} else { // geometria bistecca-like
		runManager->SetUserInitialization(new B1DetectorConstructionSTEAK(sphereDistY, CenterSphere, DetConf, AbsorberThickness, AbsorberMaterial, FilterFlag, SphereSelect, IsotopeChoice, QuickFlag, PixelThickness));
	}
		
	// Physics list
	//G4VModularPhysicsList* physicsList = new QBBC;
	//physicsList->SetVerboseLevel(1);
	
	//  runManager->SetUserInitialization(new B1PhysicsList);
	
	B1PhysicsList* physicsList=new B1PhysicsList;
	physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	runManager->SetUserInitialization(physicsList);
	
	// User action initialization
	//	runManager->SetUserInitialization(new B1ActionInitialization(sphereDistY, CenterSphere, CollHoleDiam, FilterFlag, primFile, TBRvalue,SphereSelect, SphereSelect));
	runManager->SetUserInitialization(new B1ActionInitialization(sphereDistY, CenterSphere, DetConf, FilterFlag, TBRvalue, SphereSelect, IsotopeChoice, OutFileName));
	
	// Initialize visualization
	//
	G4VisManager* visManager = new G4VisExecutive;
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();
	
	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	
	runManager->Initialize();

	if ( VisFlag ) { //Prepare for vis
		ui = new G4UIExecutive(argc, argv);
	}

	if ( ! ui ) {
		// batch mode
		if (MacroName!="") {
			G4String command = "/control/execute ";
			UImanager->ApplyCommand(command+MacroName);
		} else {
			UImanager->ApplyCommand("/gps/pos/type Volume");
			UImanager->ApplyCommand("/gps/pos/centre 0. 0. 0. cm");
			UImanager->ApplyCommand("/gps/pos/shape Para");
			
			if (SphereSelect>0) { //NEMA like
				if (SourceChoice==2) { //fondo
					UImanager->ApplyCommand("/gps/pos/halfx 30 cm");
					UImanager->ApplyCommand("/gps/pos/halfy 30 cm");
					UImanager->ApplyCommand("/gps/pos/halfz 30 cm");
					UImanager->ApplyCommand("/gps/pos/confine Phantom");
				} else { //segnale
					UImanager->ApplyCommand("/gps/pos/halfx 4 cm");
					UImanager->ApplyCommand("/gps/pos/halfy 14 cm");
					UImanager->ApplyCommand("/gps/pos/halfz 2 cm");
						UImanager->ApplyCommand("/gps/pos/confine Sphere");
				}
				
			}
			
			else { //STEAK like
				
				if (SourceChoice==2) { //fondo
					UImanager->ApplyCommand("/gps/pos/halfx 15 cm");
					UImanager->ApplyCommand("/gps/pos/halfy 15 cm");
					UImanager->ApplyCommand("/gps/pos/halfz 2 cm");
					UImanager->ApplyCommand("/gps/pos/confine Phantom");
				} else { //segnale
					UImanager->ApplyCommand("/gps/pos/halfx 4 cm");
					UImanager->ApplyCommand("/gps/pos/halfy 7 cm");
					UImanager->ApplyCommand("/gps/pos/halfz 2 cm");
						UImanager->ApplyCommand("/gps/pos/confine Sphere");
				}
				
			}

			if (IsotopeChoice<=4) {
			UImanager->ApplyCommand("/gps/ene/type/Mono");
			UImanager->ApplyCommand("/gps/ene/mono 0 MeV");
			} else {
				G4cout<<("AAAA FOTONE")<<G4endl;
				UImanager->ApplyCommand("/gps/ene/type/Mono");
				UImanager->ApplyCommand("/gps/ene/mono 80 keV");
				UImanager->ApplyCommand("/gps/ang/type iso");
			}
			UImanager->ApplyCommand("/tracking/verbose " + std::to_string(Verbose));
			UImanager->ApplyCommand("/run/beamOn " + std::to_string(NoOfPrimToGen));
		}
		
	}
	else {
		UImanager->ApplyCommand("/gps/pos/type Volume");
		UImanager->ApplyCommand("/gps/pos/centre 0. 0. 0. cm");
		UImanager->ApplyCommand("/gps/pos/shape Para");
		
		if (SphereSelect>0) { //NEMA like
			if (SourceChoice==2) { //fondo
				UImanager->ApplyCommand("/gps/pos/halfx 30 cm");
				UImanager->ApplyCommand("/gps/pos/halfy 30 cm");
				UImanager->ApplyCommand("/gps/pos/halfz 30 cm");
				UImanager->ApplyCommand("/gps/pos/confine Phantom");
			} else { //segnale
				UImanager->ApplyCommand("/gps/pos/halfx 4 cm");
				UImanager->ApplyCommand("/gps/pos/halfy 14 cm");
				UImanager->ApplyCommand("/gps/pos/halfz 2 cm");
					UImanager->ApplyCommand("/gps/pos/confine Sphere");
			}
			
		}
		
		else { //STEAK like
			
			if (SourceChoice==2) { //fondo
				UImanager->ApplyCommand("/gps/pos/halfx 15 cm");
				UImanager->ApplyCommand("/gps/pos/halfy 15 cm");
				UImanager->ApplyCommand("/gps/pos/halfz 2 cm");
				UImanager->ApplyCommand("/gps/pos/confine Phantom");
			} else { //segnale
				UImanager->ApplyCommand("/gps/pos/halfx 4 cm");
				UImanager->ApplyCommand("/gps/pos/halfy 7 cm");
				UImanager->ApplyCommand("/gps/pos/halfz 2 cm");
					UImanager->ApplyCommand("/gps/pos/confine Sphere");
			}
			
		}

		if (IsotopeChoice<=4) {
		UImanager->ApplyCommand("/gps/ene/type/Mono");
		UImanager->ApplyCommand("/gps/ene/mono 0 MeV");
		} else {
			UImanager->ApplyCommand("/gps/ene/type/Mono");
			UImanager->ApplyCommand("/gps/ene/mono 80 keV");
			UImanager->ApplyCommand("/gps/ang/type iso");
		}
		UImanager->ApplyCommand("/tracking/verbose " + std::to_string(Verbose));
		UImanager->ApplyCommand("/run/beamOn " + std::to_string(NoOfPrimToGen));
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete ui;
	}
	
#ifdef G4MULTITHREADED
	// Manually merge root files at the end (and delete temporary root files) since the automatic G4 way fails to merge vectors
	G4cout<<"hadd -f "<< OutFileName<<".root ";
	G4String tempFileList="";
	for (int iThread=0; iThread<numOfThreads; iThread++) {
		tempFileList.append(OutFileName+ "_t"+to_string(iThread) + ".root ");
	}
	G4String comandoHadd="hadd -f  " + OutFileName + ".root " + tempFileList;
	G4String comandoRm="rm " + tempFileList;

	system(comandoHadd);
	system(comandoRm);

#endif
	

	delete visManager;
	delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
