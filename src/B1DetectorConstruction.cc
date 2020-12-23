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
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class
///
///
///

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4StepLimiter.hh"
#include "G4UserLimits.hh"
#include "G4Region.hh"
#include "G4EllipticalTube.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction(G4double x0, G4double CenterSphere, G4double CollHoleDiam, G4double CollThickness, G4int CollMaterial, G4int FilterFlag, G4double SphereSelect, G4int IsotopeChoice, G4bool QuickFlag, G4double PxThickness)
: G4VUserDetectorConstruction(),
fScoringVolume(0), fX0Scan(x0), fCenterSphere(CenterSphere), fCollHoleDiam(CollHoleDiam), fCollThickness(CollThickness), fCollMaterial(CollMaterial), fFilterFlag(FilterFlag), fSphereSelect(SphereSelect), fIsotopeChoice(IsotopeChoice), fQuickFlag(QuickFlag), fPixelThickness(PxThickness)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();
	
	// Option to switch on/off checking of volumes overlaps
	G4bool checkOverlaps = false;
	
	//###################################################################
	//###################################################
	// MATERIALS DEFINITION
	//##########################
	
	G4double z, a, density;
	G4String name, symbol;
	G4int ncomponents, natoms;
	
	// Some usefuk elements
	a = 1.01*g/mole;
	G4Element* elH = new G4Element (name="Hydrogen", symbol="H", z=1.,a );
	a = 12.01*g/mole;
	G4Element* elC = new G4Element (name="Carbon", symbol="C", z=6.,a );
	a = 16.00*g/mole;
	G4Element* elO = new G4Element (name="Oxygen", symbol="O", z=8.,a );
	a = 14.00*g/mole;
	G4Element* elN = new G4Element (name="Nitrogen", symbol="N", z=7.,a );
	G4Element* elSi = new G4Element("Silicon" ,"Si" , 14., 28.09*g/mole);

	G4double densityAlu = 2.600*g/cm3;
	G4NistManager::Instance()->BuildMaterialWithNewDensity("MyAlu","G4_Al",densityAlu);
	G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");


	//###################################################
	// Glass resin for sensor filter
	//##########################
	density = 4.000*g/cm3; //4 for MT9V011, 2.43 for MT9V115
	if (fIsotopeChoice==2) density=2.43*g/cm3;
	G4Material* Resin = new G4Material (name="Resin", density, ncomponents=3);
	Resin->AddElement (elH, natoms=30);
	Resin->AddElement (elC, natoms=20);
	Resin->AddElement (elO, natoms=2);
	
	//###################################################
	// AGAR AGAR Source - AgarAgar should be C14 H24 O9
	//##########################
	G4double Agardensity = 1.030*g/cm3;
	G4Material* AgarAgar = new G4Material (name="AgarAgar", Agardensity, ncomponents=3);
	AgarAgar->AddElement (elH, natoms=24);
	AgarAgar->AddElement (elC, natoms=14);
	AgarAgar->AddElement (elO, natoms=9);
	
	//Define PMMA (C502H8)
	// NIST reference
 G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3);
 PMMA -> AddElement(elC, 5);
 PMMA -> AddElement(elO, 2);
 PMMA -> AddElement(elH, 8);
	

	
	//###################################################################
	//###################################################
	// MATERIAL'S ASSIGNMENT
	//##########################
	
//	G4Material* SourceExtY_mat = AgarAgar;
//	G4Material* ABSaround_mat = ABS;
//	G4Material* ABSbehind_mat = ABS;
//	G4Material* SourceExtSR_mat = nist->FindOrBuildMaterial("MyAlu");
//	G4Material* Resin_mat = Resin;
//	G4Material* shapeColl_mat = nist->FindOrBuildMaterial("G4_Cu");
//	G4Material* shapeDummy_mat = world_mat;
//	G4Material* pix_mat = nist->FindOrBuildMaterial("G4_Si");
//	G4Material* Cmos_mat = pix_mat;
//	G4Material* carrier_mat = FR4;
//	//	carrier_mat=world_mat; //to remove carrier behind CMOS
//	if (fCollMaterial==2) 	shapeColl_mat=nist->FindOrBuildMaterial("MyAlu");
//	else if (fCollMaterial==3) shapeColl_mat=ABS;
//	G4Material* GammaSource_mat = polycarbonate;
//	G4Material* SourcePSR_mat=nist->FindOrBuildMaterial("MyPlastic");
//	if (fSphereSelect==8 || fSphereSelect==9) SourceExtSR_mat=world_mat;
//
//	G4Material* Na22nudeSource_mat = nist->FindOrBuildMaterial("MyAlu");

	G4Material* PhantomShell_mat = PMMA;
	G4Material* PhantomInt_mat =nist->FindOrBuildMaterial("G4_WATER");
	G4Material* SphereShell_mat=nist->FindOrBuildMaterial("G4_GLASS_PLATE");
	G4Material* SphereInt_mat =nist->FindOrBuildMaterial("G4_WATER");

	G4double phantom_sizeX = 30*cm;
	G4double phantom_sizeY = 22*cm;
	G4double phantom_sizeZ  = 18*cm;
	G4double steak_thickness= 3*cm;
	
	G4double phantom_thickness = 2.5*mm;
	G4double sphere_thickness = 5*mm;
	
	G4int whichSphere=abs(G4int(fSphereSelect))-1;

	G4double sphere_diameter[6]= {
		10*mm, 13*mm, 17*mm, 22*mm, 28*mm, 37*mm
	};
	G4double distSphereSurface=70*mm;
	G4double distGammaCamera=25*cm;
	
	G4double sphereCenterZ=0.5*phantom_sizeZ-distSphereSurface-0.5*sphere_diameter[whichSphere];
	if (fCenterSphere==0) sphereCenterZ=0*cm; //sfera centrata
	G4ThreeVector spherePos= G4ThreeVector(0,0,sphereCenterZ);
	G4ThreeVector gammaCameraPos= G4ThreeVector(0,phantom_sizeY*0.5+distGammaCamera,0);
	
	G4double gammaCamera_sizeXZ=100*cm;
	G4double gammaCamera_sizeY=10*um;

	
	
	//###################################################################
	//###################################################
	// DEFINITIONS OF VOLUMES
	//##########################
	//###################################################################
	
	//###################################################
	// WORLD
	//##########################
	
	G4double world_sizeXY = 1.*m;
	G4double world_sizeZ  = 1.*m;
	
	G4Box* solidWorld =
	new G4Box("World",                       //its name
						0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
	
	G4LogicalVolume* logicWorld =
	new G4LogicalVolume(solidWorld,          //its solid
											world_mat,           //its material
											"World");            //its name
	
	G4VPhysicalVolume* physWorld =
	new G4PVPlacement(0,                     //no rotation
										G4ThreeVector(),       //at (0,0,0)
										logicWorld,            //its logical volume
										"World",               //its name
										0,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	
	//###################################################
	// Phantom
	//##########################
	

	
	G4VSolid* solidPhantomLarge;
	G4VSolid* solidPhantomInt;
	G4VSolid* solidSphereLarge;
	G4VSolid* solidSphereInt;
	G4RotationMatrix* phantomRot= new G4RotationMatrix();

	if (fSphereSelect>0){ //fantoccio simil NEMA
		solidPhantomLarge =	new G4EllipticalTube("PhantomLarge",                       //its name
						0.5*phantom_sizeX+phantom_thickness, 0.5*phantom_sizeY+phantom_thickness, 0.5*phantom_sizeZ+phantom_thickness);     //its size
		
		solidPhantomInt =
		new G4EllipticalTube("PhantomInt",                       //its name
							0.5*phantom_sizeX, 0.5*phantom_sizeY, 0.5*phantom_sizeZ);     //its size
		
		solidSphereLarge= new G4Orb("SphereLarge", 0.5*sphere_diameter[whichSphere]+sphere_thickness);
		
		solidSphereInt= new G4Orb("SphereInt", 0.5*sphere_diameter[whichSphere]);
		
	} else { //fantoccio simil bistecca
		solidPhantomLarge =	new G4Tubs("PhantomLarge",0,0.5*phantom_sizeX+phantom_thickness, 0.5*steak_thickness+phantom_thickness,0*deg,360*deg);     //its size
		
		solidPhantomInt =		new G4Tubs("PhantomInt",0,0.5*phantom_sizeX,0.5*steak_thickness,0*deg,360*deg);     //its size
		
		solidSphereLarge= new G4Tubs("SphereLarge",0,0.5*sphere_diameter[whichSphere]+sphere_thickness, 0.5*steak_thickness+phantom_thickness,0,360*deg);
		
		solidSphereInt= new G4Tubs("SphereInt",0,0.5*sphere_diameter[whichSphere],0.5*steak_thickness,0,360*deg);
		
		phantomRot->rotateX(90*deg);
	}

	G4SubtractionSolid* solidSphereShell=new G4SubtractionSolid("SphereLarge-SphereInt", solidSphereLarge, solidSphereInt);

	
	G4SubtractionSolid* solidPhantomShell=new G4SubtractionSolid("PhantomLarge-PhantomInt", solidPhantomLarge, solidPhantomInt);
	
	G4SubtractionSolid* solidPhantom=new G4SubtractionSolid("PhantomInt-SphereLarge", solidPhantomInt, solidSphereLarge, new G4RotationMatrix, spherePos);

	
	
	G4LogicalVolume* logicPhantomShell =
	new G4LogicalVolume(solidPhantom,          //its solid
											PhantomShell_mat,           //its material
											"PhantomShell");            //its name
	
	G4VPhysicalVolume* physPhantomShell =
	new G4PVPlacement(phantomRot,                     //no rotation
										G4ThreeVector(),       //at (0,0,0)
										logicPhantomShell,            //its logical volume
										"PhantomShell",               //its name
										logicWorld,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	
	
	G4LogicalVolume* logicPhantom =
	new G4LogicalVolume(solidPhantom,          //its solid
											PhantomInt_mat,           //its material
											"Phantom");            //its name
	
	G4VPhysicalVolume* physPhantomInt =
	new G4PVPlacement(phantomRot,                     //no rotation
										G4ThreeVector(),       //at (0,0,0)
										logicPhantom,            //its logical volume
										"Phantom",               //its name
										logicWorld,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	
	//###################################################
	// Sphere
	//##########################
	

	
	
	G4LogicalVolume* logicSphereShell =
	new G4LogicalVolume(solidSphereShell,          //its solid
											SphereShell_mat,           //its material
											"SphereShell");            //its name
	
	G4VPhysicalVolume* physSphereShell =
	new G4PVPlacement(phantomRot,                     //no rotation
										spherePos,       //at (0,0,0)
										logicSphereShell,            //its logical volume
										"SphereShell",               //its name
										logicWorld,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	G4LogicalVolume* logicSphere =
	new G4LogicalVolume(solidSphereInt,          //its solid
											SphereInt_mat,           //its material
											"Sphere");            //its name
	
	G4VPhysicalVolume* physSphere =
	new G4PVPlacement(phantomRot,                     //no rotation
										spherePos,       //at (0,0,0)
										logicSphere,            //its logical volume
										"Sphere",               //its name
										logicWorld,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	
	
	
	G4Box* solidGammaCamera =
	new G4Box("GammmaCamera",                       //its name
						0.5*gammaCamera_sizeXZ, 0.5*gammaCamera_sizeY, 0.5*gammaCamera_sizeXZ);     //its size
	
	G4LogicalVolume* logicGammaCamera =
	new G4LogicalVolume(solidGammaCamera,          //its solid
											world_mat,           //its material
											"GammaCamera");            //its name
	
	G4VPhysicalVolume* physGammaCamera =
	new G4PVPlacement(0,                     //no rotation
										gammaCameraPos,       //at (0,0,0)
										logicGammaCamera,            //its logical volume
										"GammaCamera",               //its name
										logicWorld,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking

	
	
	
	// Set scoring volume
	//Pixelated CMOS
//	fScoringVolume = logicPhantom;
	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
