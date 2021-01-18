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
// $Id: B1DetectorConstructionSTEAK.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstructionSTEAK.cc
/// \brief Implementation of the B1DetectorConstructionSTEAK class
///
///
///

#include "B1DetectorConstructionSTEAK.hh"

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
#include "G4ExtrudedSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstructionSTEAK::B1DetectorConstructionSTEAK(G4double sphereDistY, G4double CenterSphere, G4double DetConf, G4double CollThickness, G4int CollMaterial, G4int FilterFlag, G4double SphereSelect, G4int IsotopeChoice, G4bool QuickFlag, G4double PxThickness)
: G4VUserDetectorConstruction(),
fScoringVolume(0), fsphereDistY(sphereDistY), fCenterSphere(CenterSphere), fDetConf(DetConf), fCollThickness(CollThickness), fCollMaterial(CollMaterial), fFilterFlag(FilterFlag), fSphereSelect(SphereSelect), fIsotopeChoice(IsotopeChoice), fQuickFlag(QuickFlag), fPixelThickness(PxThickness)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstructionSTEAK::~B1DetectorConstructionSTEAK()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstructionSTEAK::Construct()
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

	G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

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
	

	G4Material* PhantomShell_mat = PMMA;
	G4Material* PhantomInt_mat =nist->FindOrBuildMaterial("G4_WATER");
	G4Material* SphereShell_mat=nist->FindOrBuildMaterial("G4_GLASS_PLATE");
	G4Material* SphereInt_mat =nist->FindOrBuildMaterial("G4_WATER");
	G4Material* Collimator_mat =nist->FindOrBuildMaterial("G4_Pb");

	G4double phantom_sizeX = 30*cm;
	G4double phantom_sizeY = 22*cm;
	G4double phantom_sizeZ  = 18*cm;
	G4double steak_thickness= 3*cm;
	
	G4double vesselExt_sizeX=25*cm;
	G4double vesselExt_sizeY=25*cm;
	G4double vesselExt_sizeZ=3*cm;

	G4double vesselExt_thickness=1000*um;
	G4double vesselInt_thickness=500*um;

	G4double triangle_vertices=3*cm;
	
//	if (fDetConf!=0) phantom_sizeX=phantom_sizeY;

	G4double phantom_thickness = 2.5*mm;
	G4double sphere_thickness = 5*mm;
	
	G4int whichSphere=abs(G4int(fSphereSelect))-1;


	G4double sphere_diameter[6]= {
		10*mm, 13*mm, 17*mm, 22*mm, 28*mm, 37*mm
	};
	
	
	G4double distSphereSurface=70*mm;
	G4double distGammaCamera=25*cm;
	G4double detLinearDistance=7*cm;
	G4double distGammaCameraCollimator=5*cm;
	
	G4double sphereCenterZ=0.5*phantom_sizeZ-distSphereSurface-0.5*sphere_diameter[whichSphere];
	if (fCenterSphere==0) sphereCenterZ=0*cm; //sfera centrata

	G4double sphereCenterY=0*cm;


	G4double gammaCamera_sizeXY=100*cm;
	G4double gammaCamera_sizeZ=10*um;
	
	G4double collimator_sizeXY=100*cm;
	G4double collimator_sizeZ=40.64*mm;
	
	G4double collimatorHole_Diameter=2.94*mm;
	G4double collimatorHole_SeptalThickness=1.13*mm*0.5;

	G4double safeDistance=5*mm;

//	G4double collimatorPosZ=1;
	
	G4double collimatorPosZ=vesselExt_sizeZ*0.5+distGammaCamera-collimator_sizeZ*0.5;
	
	G4ThreeVector spherePos= G4ThreeVector(0,sphereCenterY,sphereCenterZ);
	G4ThreeVector collimatorPos= G4ThreeVector(0,0, collimatorPosZ);
	G4ThreeVector gammaCameraPos= G4ThreeVector(0,0, vesselExt_sizeZ*0.5+distGammaCamera);
	G4ThreeVector trianglePos= G4ThreeVector(0, -triangle_vertices*0.5, 0);


	//###################################################################
	//###################################################
	// DEFINITIONS OF VOLUMES
	//##########################
	//###################################################################
	
	//###################################################
	// WORLD
	//##########################
	
	G4double world_sizeXY = 2.*m;
	G4double world_sizeZ  = 2.*m;
	
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
	
	G4double booleanDistance=0*um;
	
	G4VSolid* solidOuterVesselLarge;
	G4VSolid* solidOuterVesselInt;
	G4VSolid* solidInnerVesselLarge;
	G4VSolid* solidInnerVesselInt;
	G4VSolid* solidSphereLarge;
	G4VSolid* solidSphereInt;

	solidOuterVesselLarge=new G4Box("OuterVesselLarge",
															(vesselExt_sizeX+vesselExt_thickness)*0.5,
															(vesselExt_sizeY+vesselExt_thickness)*0.5,
															vesselExt_sizeZ*0.5
															);

		
	
	solidOuterVesselInt=new G4Box("OuterVesselInt",
															(vesselExt_sizeX)*0.5,
															(vesselExt_sizeY)*0.5,
															vesselExt_sizeZ*0.5+booleanDistance
															);
	
	
	std::vector<G4TwoVector> verticiExt;
	std::vector<G4TwoVector> verticiInt;
	std::vector<G4ExtrudedSolid::ZSection> pianiExt;
	std::vector<G4ExtrudedSolid::ZSection> pianiInt;

	if (whichSphere==6) { //triangolo
	verticiExt.push_back({0,2*triangle_vertices+vesselInt_thickness});
	verticiExt.push_back({triangle_vertices+vesselInt_thickness,-(triangle_vertices+vesselInt_thickness)});
	verticiExt.push_back({-(triangle_vertices+vesselInt_thickness),-(triangle_vertices+vesselInt_thickness)});
	

	verticiInt.push_back({0,2*triangle_vertices});
	verticiInt.push_back({triangle_vertices,-triangle_vertices});
	verticiInt.push_back({-triangle_vertices,-triangle_vertices});

	pianiExt.push_back(G4ExtrudedSolid::ZSection(-0.5*vesselExt_sizeZ-booleanDistance,G4TwoVector(0,0),1));
	pianiExt.push_back(G4ExtrudedSolid::ZSection(0.5*vesselExt_sizeZ+booleanDistance,G4TwoVector(0,0),1));

	
	pianiInt.push_back(G4ExtrudedSolid::ZSection(-0.5*vesselExt_sizeZ-booleanDistance,G4TwoVector(0,0),1));
	pianiInt.push_back(G4ExtrudedSolid::ZSection(0.5*vesselExt_sizeZ+booleanDistance,G4TwoVector(0,0),1));
	
	
	
	solidInnerVesselLarge= new G4ExtrudedSolid("InnerVesselExt",
																						 verticiExt,
																						 pianiExt
																						 );

	solidInnerVesselInt= new G4ExtrudedSolid("InnerVesselInt",
																						 verticiInt,
																						 pianiInt
																						 );
	} else if (whichSphere==7) { //Sauron
		double triangleSizeExt=triangle_vertices+vesselInt_thickness;
		double lenghtFactor=2.5;
		verticiExt.push_back({0,lenghtFactor*triangleSizeExt});
		verticiExt.push_back({triangleSizeExt,0});
		verticiExt.push_back({0,-lenghtFactor*triangleSizeExt});
		verticiExt.push_back({-triangleSizeExt,0});
		
		verticiInt.push_back({0,lenghtFactor*triangle_vertices});
		verticiInt.push_back({triangle_vertices,0});
		verticiInt.push_back({0,-lenghtFactor*triangle_vertices});
		verticiInt.push_back({-triangle_vertices,0});

		pianiExt.push_back(G4ExtrudedSolid::ZSection(-0.5*vesselExt_sizeZ-booleanDistance,G4TwoVector(0,0),1));
		pianiExt.push_back(G4ExtrudedSolid::ZSection(0.5*vesselExt_sizeZ+booleanDistance,G4TwoVector(0,0),1));

		
		pianiInt.push_back(G4ExtrudedSolid::ZSection(-0.5*vesselExt_sizeZ-booleanDistance,G4TwoVector(0,0),1));
		pianiInt.push_back(G4ExtrudedSolid::ZSection(0.5*vesselExt_sizeZ+booleanDistance,G4TwoVector(0,0),1));
		
		
		
		solidInnerVesselLarge= new G4ExtrudedSolid("InnerVesselExt",
																							 verticiExt,
																							 pianiExt
																							 );

		solidInnerVesselInt= new G4ExtrudedSolid("InnerVesselInt",
																							 verticiInt,
																							 pianiInt
																							 );
	
	} else if (whichSphere==8) { //Papillon
		double triangleSizeExt=triangle_vertices+vesselInt_thickness;
		double lenghtFactor=2;

		verticiExt.push_back({triangleSizeExt,lenghtFactor*triangleSizeExt}); //A
		verticiExt.push_back({-triangleSizeExt,lenghtFactor*triangleSizeExt}); //D
		verticiExt.push_back({-vesselInt_thickness,0}); //O
		verticiExt.push_back({-triangleSizeExt,-lenghtFactor*triangleSizeExt}); //C
		verticiExt.push_back({triangleSizeExt,-lenghtFactor*triangleSizeExt}); //B
		verticiExt.push_back({vesselInt_thickness,0}); //O

		verticiInt.push_back({triangle_vertices,lenghtFactor*triangle_vertices}); //A
		verticiInt.push_back({-triangle_vertices,lenghtFactor*triangle_vertices}); //D
		verticiInt.push_back({-0.0001,0}); //O
		verticiInt.push_back({-triangle_vertices,-lenghtFactor*triangle_vertices}); //C
		verticiInt.push_back({triangle_vertices,-lenghtFactor*triangle_vertices}); //B
		verticiInt.push_back({0.0001,0}); //O

		pianiExt.push_back(G4ExtrudedSolid::ZSection(-0.5*vesselExt_sizeZ-booleanDistance,G4TwoVector(0,0),1));
		pianiExt.push_back(G4ExtrudedSolid::ZSection(0.5*vesselExt_sizeZ+booleanDistance,G4TwoVector(0,0),1));

		
		pianiInt.push_back(G4ExtrudedSolid::ZSection(-0.5*vesselExt_sizeZ-booleanDistance,G4TwoVector(0,0),1));
		pianiInt.push_back(G4ExtrudedSolid::ZSection(0.5*vesselExt_sizeZ+booleanDistance,G4TwoVector(0,0),1));
		
		
		
		solidInnerVesselLarge= new G4ExtrudedSolid("InnerVesselExt",
																							 verticiExt,
																							 pianiExt
																							 );

		solidInnerVesselInt= new G4ExtrudedSolid("InnerVesselInt",
																							 verticiInt,
																							 pianiInt
																							 );
	}	else  {
	solidInnerVesselLarge= new G4Tubs("InnerVesselExt",0,
																		0.5*sphere_diameter[whichSphere]+sphere_thickness,
																		0.5*steak_thickness,
																		0,360*deg);

	solidInnerVesselInt= new G4Tubs("InnerVesselInt",0,
																	0.5*sphere_diameter[whichSphere],
																	0.5*steak_thickness+booleanDistance,
																	0,360*deg);
	}
	
	G4SubtractionSolid* solidOuterVesselShell=new G4SubtractionSolid("OuterVesselLarge-OuterVesselInt", solidOuterVesselLarge, solidOuterVesselInt);


	
	G4SubtractionSolid* solidInnerVesselShell=new G4SubtractionSolid("InnerVesselExt-InnerVesselInt", solidInnerVesselLarge, solidInnerVesselInt);
	
	G4SubtractionSolid* solidOuterVessel=new G4SubtractionSolid("OuterVesselInt-InnerVesselExt", solidOuterVesselInt, solidInnerVesselLarge);
	
//	G4SubtractionSolid* solidOuterVessel=new G4SubtractionSolid("OuterVesselInt-InnerVesselExt", solidOuterVesselInt, solidInnerVesselShell);

	G4LogicalVolume* logicOuterVesselShell =
	new G4LogicalVolume(solidOuterVesselShell,          //its solid
											PhantomShell_mat,           //its material
											"PhantomShell");            //its name
	
	new G4PVPlacement(0,                     //no rotation
										G4ThreeVector(),       //at (0,0,0)
										logicOuterVesselShell,            //its logical volume
										"PhantomShell",               //its name
										logicWorld,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking
	


	G4LogicalVolume* logicOuterVessel =
	new G4LogicalVolume(solidOuterVessel,          //its solid
											PhantomInt_mat,           //its material
											"Phantom");            //its name

	new G4PVPlacement(0,                     //no rotation
										G4ThreeVector(),       //at (0,0,0)
										logicOuterVessel,            //its logical volume
										"Phantom",               //its name
										logicWorld,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking


	G4LogicalVolume* logicInnerVesselShell =
	new G4LogicalVolume(solidInnerVesselShell,          //its solid
											SphereShell_mat,           //its material
											"SphereShell");            //its name

	new G4PVPlacement(0,                     //no rotation
										G4ThreeVector(),       //at (0,0,0)
										logicInnerVesselShell,            //its logical volume
										"SphereShell",               //its name
										logicWorld,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking

	G4LogicalVolume* logiInnerVessel =
	new G4LogicalVolume(solidInnerVesselInt,          //its solid
											SphereInt_mat,           //its material
											"Sphere");            //its name

	new G4PVPlacement(0,                     //no rotation
										G4ThreeVector(),       //at (0,0,0)
										logiInnerVessel,            //its logical volume
										"Sphere",               //its name
										logicWorld,                     //its mother  volume
										false,                 //no boolean operation
										0,                     //copy number
										checkOverlaps);        //overlaps checking


	G4Box* solidGammaCamera =
	new G4Box("GammmaCamera",                       //its name
						0.5*gammaCamera_sizeXY, 0.5*gammaCamera_sizeXY, 0.5*gammaCamera_sizeZ);     //its size

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
//
//
//	G4Box* solidCollimatorAbs =
//	new G4Box("CollimatorAbs",                       //its name
//						0.5*collimator_sizeXY, 0.5*collimator_sizeXY, 0.5*collimator_sizeZ);     //its size

//	G4LogicalVolume* logicCollimatorAbs =
//	new G4LogicalVolume(solidCollimatorAbs,          //its solid
//											Collimator_mat,           //its material
//											"CollimatorAbs");            //its name
//
//	new G4PVPlacement(0,                     //no rotation
//										collimatorPos,       //at (0,0,0)
//										logicCollimatorAbs,            //its logical volume
//										"CollimatorAbs",               //its name
//										logicWorld,                     //its mother  volume
//										false,                 //no boolean operation
//										0,                     //copy number
//										checkOverlaps);        //overlaps checking

	
	G4double collimatorHoleAext=(collimatorHole_Diameter+collimatorHole_SeptalThickness)/2.;
	G4double collimatorHoleBext=(collimatorHole_Diameter+collimatorHole_SeptalThickness)*sqrt(3)/3;
	G4double collimatorHoleCext=(collimatorHole_Diameter+collimatorHole_SeptalThickness)/2.*sqrt(3)/3;

	
	std::vector<G4TwoVector> verticiHoleExt;

	verticiHoleExt.push_back({collimatorHoleCext,-collimatorHoleAext}); //vertex F
	verticiHoleExt.push_back({-collimatorHoleCext,-collimatorHoleAext}); //vertex E
	verticiHoleExt.push_back({-collimatorHoleBext,0}); //vertex D
	verticiHoleExt.push_back({-collimatorHoleCext,collimatorHoleAext}); //vertex C
	verticiHoleExt.push_back({collimatorHoleCext,collimatorHoleAext}); //vertex B
	verticiHoleExt.push_back({collimatorHoleBext,0}); //vertex A


	std::vector<G4ExtrudedSolid::ZSection> pianiHoleExt;
	
	pianiHoleExt.push_back(G4ExtrudedSolid::ZSection(-0.5*collimator_sizeZ,G4TwoVector(0,0),1));
	pianiHoleExt.push_back(G4ExtrudedSolid::ZSection(0.5*collimator_sizeZ,G4TwoVector(0,0),1));

	
	G4VSolid* solidCollimatorHoleExt= new G4ExtrudedSolid("CollimatorHoleExt",
																										 verticiHoleExt,
																										 pianiHoleExt
																						 );
	
	G4double collimatorHoleA=collimatorHole_Diameter/2.;
	G4double collimatorHoleB=collimatorHole_Diameter*sqrt(3)/3;
	G4double collimatorHoleC=collimatorHole_Diameter/2.*sqrt(3)/3;

	
	std::vector<G4TwoVector> verticiHole;

	verticiHole.push_back({collimatorHoleC,-collimatorHoleA}); //vertex F
	verticiHole.push_back({-collimatorHoleC,-collimatorHoleA}); //vertex E
	verticiHole.push_back({-collimatorHoleB,0}); //vertex D
	verticiHole.push_back({-collimatorHoleC,collimatorHoleA}); //vertex C
	verticiHole.push_back({collimatorHoleC,collimatorHoleA}); //vertex B
	verticiHole.push_back({collimatorHoleB,0}); //vertex A


	std::vector<G4ExtrudedSolid::ZSection> pianiHole;
	
	pianiHole.push_back(G4ExtrudedSolid::ZSection(-0.5*collimator_sizeZ-100*um,G4TwoVector(0,0),1));
	pianiHole.push_back(G4ExtrudedSolid::ZSection(0.5*collimator_sizeZ+100*um,G4TwoVector(0,0),1));

	
	G4VSolid* solidCollimatorHole= new G4ExtrudedSolid("CollimatorHole",
																										 verticiHole,
																										 pianiHole
																						 );

	G4SubtractionSolid* solidCollimator=new G4SubtractionSolid("CollimatorHoleExt-CollimatorHole", solidCollimatorHoleExt, solidCollimatorHole);
	
	G4int nHoleColumn=100;
	G4int nHoleRow=100;

	G4LogicalVolume* logicCollimator =
	new G4LogicalVolume(solidCollimator,          //its solid
											Collimator_mat,           //its material
											"Collimator");            //its name

	for (int iHoleColumn=0; iHoleColumn<nHoleColumn; iHoleColumn++) {
		for (int iHoleRow=0; iHoleRow<nHoleRow; iHoleRow++) {
			
			//Verso Alto Destra
			new G4PVPlacement(0,                     //no rotation
												G4ThreeVector(iHoleRow*2.*collimatorHoleAext*sqrt(3)/2.,
																			iHoleColumn*(collimatorHole_Diameter+collimatorHole_SeptalThickness)+iHoleRow*2.*collimatorHoleAext/2.,
																			collimatorPosZ),       //at (0,0,0)
												logicCollimator,            //its logical volume
												"Collimator",               //its name
												logicWorld,                     //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
			
			//Verso Alto Sinistra
			new G4PVPlacement(0,                     //no rotation
												G4ThreeVector(-(iHoleRow+1)*2.*collimatorHoleAext*sqrt(3)/2.,
																			iHoleColumn*(collimatorHole_Diameter+collimatorHole_SeptalThickness)-(iHoleRow+1)*2.*collimatorHoleAext/2.,
																			collimatorPosZ),       //at (0,0,0)
												logicCollimator,            //its logical volume
												"Collimator",               //its name
												logicWorld,                     //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
			
			
			//Verso Basso Destra
			new G4PVPlacement(0,                     //no rotation
												G4ThreeVector(iHoleRow*2.*collimatorHoleAext*sqrt(3)/2.,
																			-(iHoleColumn+1)*(collimatorHole_Diameter+collimatorHole_SeptalThickness)+iHoleRow*2.*collimatorHoleAext/2.,
																			collimatorPosZ),       //at (0,0,0)
												logicCollimator,            //its logical volume
												"Collimator",               //its name
												logicWorld,                     //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
			
			//Verso Basso Sinistra
			new G4PVPlacement(0,                     //no rotation
												G4ThreeVector(-(iHoleRow+1)*2.*collimatorHoleAext*sqrt(3)/2.,
																			-(iHoleColumn+1)*(collimatorHole_Diameter+collimatorHole_SeptalThickness)-(iHoleRow+1)*2.*collimatorHoleAext/2.,
																			collimatorPosZ),       //at (0,0,0)
												logicCollimator,            //its logical volume
												"Collimator",               //its name
												logicWorld,                     //its mother  volume
												false,                 //no boolean operation
												0,                     //copy number
												checkOverlaps);        //overlaps checking
			
			
			
		}
	}
	// Set scoring volume
	//Pixelated CMOS
//	fScoringVolume = logicPhantom;
	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
