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
// $Id: B1RunAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1SteppingAction.hh"


#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "B1Analysis.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction(G4String FileName)
: G4UserRunAction(),
fEdep("Edep", 0.),
fEdep2("Edep2", 0.),
fEdkin("Edkin", 0.)
, fFileName(FileName)

{
	// Register accumulable to the accumulable manager
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->RegisterAccumulable(fEdep);
	accumulableManager->RegisterAccumulable(fEdep2);
	accumulableManager->RegisterAccumulable(fEdkin);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{
	delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run* run)
{
	// inform the runManager to save random number seed
	G4RunManager::GetRunManager()->SetRandomNumberStore(false);
	
	// reset accumulable to their initial values
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->Reset();
	
	
	CreateHistogram();
	
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	
	nbEventInRun = run->GetNumberOfEventToBeProcessed();
//	analysisManager->FillNtupleIColumn(0,42, nbEventInRun);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
//	G4int nofEvents = run->GetNumberOfEvent();
//	if (nofEvents == 0) return;
//	
//	// Merge accumulable
//	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
//	accumulableManager->Merge();
//	
//	// Compute dose = total energy deposit in a run and its variance
//	//
//	G4double edep  = fEdep.GetValue();
//	G4double edep2 = fEdep2.GetValue();
//	
//	//G4double edkin  = fEdkin.GetValue();
//	
//	G4double rms = edep2 - edep*edep/nofEvents;
//	if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;
//	
//	const B1DetectorConstruction* detectorConstruction
//	= static_cast<const B1DetectorConstruction*>
//	(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//	G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
//	G4double dose = edep/mass;
//	G4double rmsDose = rms/mass;
//	
//	// Run conditions
//	//  note: There is no primary generator action object for "master"
//	//        run manager for multi-threaded mode.
//	const B1PrimaryGeneratorAction* generatorAction
//	= static_cast<const B1PrimaryGeneratorAction*>
//	(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
//	G4String runCondition;
//	
////	if (generatorAction)
////	{
////		const G4GeneralParticleSource* particleGun = generatorAction->GetParticleGun();
////		runCondition += particleGun->GetParticleDefinition()->GetParticleName();
////		runCondition += " of ";
////		G4double particleEnergy = particleGun->GetParticleEnergy();
////		runCondition += G4BestUnit(particleEnergy,"Energy");
////	}
//	
//	// Print
//	//
//	if (IsMaster()) {
//		G4cout
//		<< G4endl
//		<< "--------------------End of Global Run-----------------------";
//	}
//	else {
//		G4cout
//		<< G4endl
//		<< "--------------------End of Local Run------------------------";
//	}
//	
//	G4cout
//	<< G4endl
//	<< " The run consists of " << nofEvents << " "<< runCondition
//	<< G4endl
//	<< G4endl
//	// << " N. of hits in scoring volume :" << fnofHits
//	// << G4endl
//	<< " Cumulated dose per run, in scoring volume : "
//	<< G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
//	<< G4endl
//	<< "------------------------------------------------------------"
//	<< G4endl
//	<< G4endl;
//	
//	///////////////
//	// Write Histo
//	//
	WriteHistogram();
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
	fEdep  += edep;
	fEdep2 += edep*edep;
}

void B1RunAction::AddEdkin(G4double edkin)
{
	fEdkin  += edkin;
}


void B1RunAction::CreateHistogram()
{
	// Book histograms, ntuple
	//	G4cout << "##### Create analysis manager " << "  " << this << G4endl;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	
	G4cout << "Using " << analysisManager->GetType() << " analysis manager" << G4endl;
	
	// Create directories
	analysisManager->SetVerboseLevel(1);

	analysisManager->OpenFile(fFileName);

	// Creating ntuple
	
	analysisManager->CreateNtuple("B1", "physics");
	analysisManager->CreateNtuple("Source", "SourceNtuple");
//
//	analysisManager->CreateNtupleDColumn(0,"Eabs");                           //0
//	analysisManager->CreateNtupleDColumn(0,"EabsComp", fRunEAbsComp); //1
//	analysisManager->CreateNtupleDColumn(0,"PreFilterTrackN");                  //2b
//	analysisManager->CreateNtupleDColumn(0,"PreFilterPart", fRunPreFilterPart); //3b
//	analysisManager->CreateNtupleDColumn(0,"PreFilterEn", fRunPreFilterEn); //4b
	analysisManager->CreateNtupleDColumn(0,"PreCmosTrackN");                  //2
	analysisManager->CreateNtupleDColumn(0,"PreCmosPart", fRunPart); //3
	analysisManager->CreateNtupleDColumn(0,"PreCmosEn", fRunEnPre); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosX", fRunPreCmosX); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosY", fRunPreCmosY); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosZ", fRunPreCmosZ); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosDirX", fRunPreCmosDirX); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosDirY", fRunPreCmosDirY); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosDirZ", fRunPreCmosDirZ); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosVX", fRunPreCmosVX); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosVY", fRunPreCmosVY); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosVZ", fRunPreCmosVZ); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosEnPrim", fRunEnPrePrim); //4
	analysisManager->CreateNtupleDColumn(0,"PreCmosEventNum", fRunPreCmosTrackID); //4
	analysisManager->CreateNtupleIColumn(0,"PreCmosOrigReg", fRunPreCmosOrigReg); //4
	analysisManager->CreateNtupleIColumn(0,"PreCmosPrestepReg", fRunPreCmosPrestepReg); //4
	analysisManager->CreateNtupleIColumn(0,"PreCmosCreatorProcess", fRunPreCmosCreatorProcess); //4
//	analysisManager->CreateNtupleDColumn(0,"InCmosTrackN");                   //5
//	analysisManager->CreateNtupleDColumn(0,"InCmosPart", fRunPartCmos); //6
//	analysisManager->CreateNtupleDColumn(0,"InCmosEn", fRunEnCmos); //7
//	analysisManager->CreateNtupleDColumn(0,"InCmosEnPrim", fRunEnCmosPrim); //7b
//	analysisManager->CreateNtupleFColumn(0,"InCmosTime", fRunEnCmosTime); //7c
//	analysisManager->CreateNtupleDColumn(0,"InCmosX", fRunXCmos); //8
//	analysisManager->CreateNtupleDColumn(0,"InCmosY", fRunYCmos); //9
//	analysisManager->CreateNtupleDColumn(0,"InCmosZ", fRunZCmos); //10
//	analysisManager->CreateNtupleDColumn(0,"PixelID", fRunPixNo); //11
//	analysisManager->CreateNtupleDColumn(0,"PixXPos", fRunPixXpos); //13
//	analysisManager->CreateNtupleDColumn(0,"PixYPos", fRunPixYpos); //14
	analysisManager->CreateNtupleDColumn(0,"SourceX");                           //14
	analysisManager->CreateNtupleDColumn(0,"SourceY");                           //15
	analysisManager->CreateNtupleDColumn(0,"SourceZ");                           //16

	analysisManager->CreateNtupleDColumn(0,"SourceCosX", fRunCosX); //17
	analysisManager->CreateNtupleDColumn(0,"SourceCosY", fRunCosY); //18
	analysisManager->CreateNtupleDColumn(0,"SourceCosZ", fRunCosZ); //19

	analysisManager->CreateNtupleDColumn(0,"SourceEne", fRunEnGen); //20
	analysisManager->CreateNtupleDColumn(0,"SourcePart", fRunPartGen); //20
//	analysisManager->CreateNtupleDColumn(0,"SourceIsotope", fRunIsotopeGen); //21
//	analysisManager->CreateNtupleIColumn(0,"Nev");							//22
	analysisManager->CreateNtupleSColumn(0,"SourceReg");                           //17

	analysisManager->CreateNtupleDColumn(1,"AllX");                           //0
	analysisManager->CreateNtupleDColumn(1,"AllY");                           //1
	analysisManager->CreateNtupleDColumn(1,"AllZ");                           //2
	analysisManager->CreateNtupleSColumn(1,"AllReg");                           //17

	analysisManager->CreateNtupleDColumn(1,"AllCosX", fRunCosX);                           //3
	analysisManager->CreateNtupleDColumn(1,"AllCosY", fRunCosY);                           //4
	analysisManager->CreateNtupleDColumn(1,"AllCosZ", fRunCosZ);                           //5
	
	analysisManager->CreateNtupleDColumn(1,"AllEne", fRunEnGen);                           //6
	analysisManager->CreateNtupleDColumn(1,"AllPart", fRunPartGen);                           //6
	analysisManager->CreateNtupleDColumn(1,"AllIsotope", fRunIsotopeGen);                           //7
	analysisManager->CreateNtupleDColumn(1,"ExitX", fRunXExit);                           //8
	analysisManager->CreateNtupleDColumn(1,"ExitY", fRunYExit);                           //9
	analysisManager->CreateNtupleDColumn(1,"ExitZ", fRunZExit);                           //10
	analysisManager->CreateNtupleDColumn(1,"ExitCosX", fRunCosXExit);                           //11
	analysisManager->CreateNtupleDColumn(1,"ExitCosY", fRunCosYExit);                           //12
	analysisManager->CreateNtupleDColumn(1,"ExitCosZ", fRunCosZExit);                           //13
	
	analysisManager->CreateNtupleDColumn(1,"ExitEne", fRunEnExit);                           //14
	analysisManager->CreateNtupleDColumn(1,"ExitPart", fRunPartExit);                           //15
	analysisManager->CreateNtupleDColumn(1,"ExitParentID", fRunParentIDExit);                           //16
	analysisManager->CreateNtupleIColumn(1,"ExitProcess", fRunExitProcess); //17
	analysisManager->CreateNtupleDColumn(1,"ExitTrackN"); //18

	analysisManager->FinishNtuple(0);
	analysisManager->FinishNtuple(1);
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void B1RunAction::WriteHistogram()
{
	
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	
	// save histograms
	//
	analysisManager->Write();
	analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

