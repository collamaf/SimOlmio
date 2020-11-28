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
// $Id: B1PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4RunManager.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include "G4IonTable.hh"
#include "G4ChargedGeantino.hh"

#include "B1RunAction.hh"
#include "B1Analysis.hh"


#include "G4Event.hh"

#include <iostream>
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using std::ofstream;
using std::ios;
using std::endl;


B1PrimaryGeneratorAction::B1PrimaryGeneratorAction(B1EventAction* eventAction, G4double TBR, G4int SourceSelect)
: G4VUserPrimaryGeneratorAction(),
fGPS(0) ,
evtPrimAction(eventAction), fTBR(TBR), fSourceSelect(SourceSelect)

{
	G4int n_particle = 1;
	fGPS  = new G4GeneralParticleSource();

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("geantino");

	fGPS->SetParticleDefinition(particle);
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
	delete fGPS;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries (G4Event* anEvent)
{

	
	switch (fSourceSelect) {
		case 3:
			fSourceZ=39;
			fSourceA=90;
			break;

		case 4:
			fSourceZ=27;
			fSourceA=60;
			break;

		case 5:
		case 10:
			fSourceZ=11;
			fSourceA=22;
			break;
			
		case 6:
			fSourceZ=56;
			fSourceA=133;
			break;

		case 7:
			fSourceZ=55;
			fSourceA=137;
			break;
			
		default: //Sr
			fSourceZ=67;
			fSourceA=166;
			break;
	}
	
	G4double ionCharge   = 0.*eplus;
	G4double excitEnergy = 0.*keV;

	G4ParticleDefinition* sourceION
	= G4IonTable::GetIonTable()->GetIon(fSourceZ,fSourceA,excitEnergy);

	fGPS->SetParticleDefinition(sourceION);
	fGPS->SetParticleCharge(ionCharge);

	fGPS->GeneratePrimaryVertex(anEvent);

	
	
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

