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
// $Id: B3StackingAction.cc 66536 2012-12-19 14:32:36Z ihrivnac $
// 
/// \file B3StackingAction.cc
/// \brief Implementation of the B3StackingAction class

#include "B1StackingAction.hh"
#include "B1RunAction.hh"
#include "B1EventAction.hh"

#include "G4Track.hh"
#include "G4NeutrinoE.hh"
#include "G4VProcess.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1StackingAction::B1StackingAction(B1RunAction* runAction, B1EventAction* EventAction)
:runStackAction(runAction), fEventAction(EventAction)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1StackingAction::~B1StackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
B1StackingAction::ClassifyNewTrack(const G4Track* track)
{
	G4int debug=0;
	
	if (fabs(track->GetDynamicParticle() ->GetPDGcode())==12) return fKill; //kill neutrinos
	
	const G4VProcess* creator=track->GetCreatorProcess();
	std::string CreatorProcname="undefined";
	if(creator) CreatorProcname=creator->GetProcessName();
	
	
	
	if (debug) G4cout<<"CMOSDEBUG PROVA STACKING creata nuova traccia tipo= "<< track->GetDynamicParticle() ->GetPDGcode()<<", MotherIsotope Val= "<< runStackAction->GetMotherIsotope()
		<<G4endl;
	
	fEventAction->ResetPassCounterSource(); //collamaf: at each new track we reset the pass counter
	
	fEventAction->ResetPassCounterCmos(); //collamaf: at each new track we reset the pass counter
	
	
	if (track->GetTrackID()==1) {
//		G4VPhysicalVolume* sourceVol = track->GetOriginTouchableHandle()->GetVolume();

		if (debug) G4cout<<"OLMIO qui prendo la posizione = "<< track->GetDynamicParticle() ->GetPDGcode()<<", POS= "<< track->GetPosition()<<" TrackStatus= "<<track->GetTrackStatus()<<G4endl;
		
		//
//		G4String volName="undefined";
//		if (track->GetVolume()->GetName()) {
//			volName=track->GetVolume()->GetName();
//		}
		fEventAction->SetSourceX(track->GetPosition().x());
		fEventAction->SetSourceY(track->GetPosition().y());
		fEventAction->SetSourceZ(track->GetPosition().z());
//		fEventAction->SetSourceRegion( sourceVol->GetName());
		//		fEventAction->SetSourceRegion("ciao");
		//		G4cout<<"CMOSDEBUG PROVA STACKI"<<G4endl;
	}
//
	if (track->GetTrackID()==4) {
	fEventAction->SetSourceRegion( track->GetVolume()->GetName());
		if (debug) G4cout<<"OLMIO e qui prendo il volume = "<< track->GetVolume()->GetName()<<G4endl;

	}
	
	if (CreatorProcname=="RadioactiveDecayBase" && track->GetDynamicParticle() ->GetPDGcode()<9e8) {
		runStackAction->SetMotherIsotope(track->GetParentID()-1);
		(runStackAction->SetMotherEnergy(track->GetKineticEnergy()/CLHEP::keV));
		(runStackAction->SetMotherPart(track->GetDynamicParticle()->GetPDGcode()));
		(runStackAction->SetMotherTime(track->GetGlobalTime()/CLHEP::ns));
		(runStackAction->GetRunEnGen()).push_back(track->GetKineticEnergy()/CLHEP::keV);
		(runStackAction->GetRunPartGen()).push_back(track->GetDynamicParticle() ->GetPDGcode());
		(runStackAction->GetRunIsotopeGen()).push_back(track->GetParentID()-1);
		(runStackAction->GetRunCosX()).push_back(track->GetMomentumDirection().x());
		(runStackAction->GetRunCosY()).push_back(track->GetMomentumDirection().y());
		(runStackAction->GetRunCosZ()).push_back(track->GetMomentumDirection().z());
	}
	
	
	return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
