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
/// \file B4/B4a/src/DetectorConstruction.cc
/// \brief Implementation of the B4::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "FiberSD.hh"
#include "G4SDManager.hh"

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

namespace B4
{

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

// Construct: 定义材料并构建几何
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();
  return DefineVolumes();
}

// DefineMaterials: 注册需要的材料
void DetectorConstruction::DefineMaterials()
{
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_Cu");
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");

  // 额外示例材料（与原模板兼容）
  G4double a; G4double z; G4double density;
  new G4Material("liquidArgon", z=18., a=39.95*g/mole, density=1.390*g/cm3);

  // 真空 (Galactic)
  new G4Material("Galactic", z=1., a=1.01*g/mole, density=universe_mean_density,
                 kStateGas, 2.73*kelvin, 3.e-18*pascal);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

// DefineVolumes: 完整几何（World / Calorimeter / Towers / Rods / Fibers）
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters according to DREAM specifications

  G4int nofTowersX = 4;        // 4 towers in X direction
  G4int nofTowersY = 4;        // 4 towers in Y direction
  G4int nofRodsPerTowerX = 16; // 16 rods per tower in X
  G4int nofRodsPerTowerY = 16; // 16 rods per tower in Y

  
  G4double rodLength = 2.0*m;           // Copper rod length
  G4double rodCrossSection = 4.0*mm;    // 4mm x 4mm cross section
  G4double centerHoleDiameter = 2.5*mm; // Center aperture diameter
  G4double fiberDiameter = 0.8*mm;      // Fiber diameter
  
  G4double towerSpacing = 1.0*cm;       // Spacing between towers
  G4double rodSpacing = 0.1*mm;         // Spacing between rods

  // Calculate overall dimensions
  G4double towerSizeXY = nofRodsPerTowerX * rodCrossSection + 
                        (nofRodsPerTowerX - 1) * rodSpacing;
  G4double detectorSizeXY = nofTowersX * towerSizeXY + 
                           (nofTowersX - 1) * towerSpacing;
  G4double worldSizeXY = 1.5 * detectorSizeXY;
  G4double worldSizeZ = 1.5 * rodLength;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto copperMaterial = G4Material::GetMaterial("G4_Cu");
  auto scintillatorMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto quartzMaterial = G4Material::GetMaterial("G4_SILICON_DIOXIDE");
  

  if ( ! defaultMaterial || ! copperMaterial || ! scintillatorMaterial || ! quartzMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name

  auto worldPV
    = new G4PVPlacement(
                 nullptr,  // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 nullptr,         // its mother  volume
                 false,           // no boolean operation
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps

  G4RotationMatrix* detectorRotation = new G4RotationMatrix();
  detectorRotation->rotateY(2.0*degree);   // 水平旋转 +2°
  detectorRotation->rotateX(0.7*degree);   // 俯仰 +0.7°

  //
  // Detector envelope
  //
  auto detectorS
    = new G4Box("Detector",     // its name
                 detectorSizeXY/2, detectorSizeXY/2, rodLength/2); // its size

  auto detectorLV
    = new G4LogicalVolume(
                 detectorS,     // its solid
                 defaultMaterial,  // its material
                 "Detector");   // its name

  new G4PVPlacement(
                 detectorRotation,  // rotation
                 G4ThreeVector(),          // at (0,0,0)
                 detectorLV,                  // its logical volume
                 "Detector",            // its name
                 worldLV,                  // its mother  volume
                 false,                    // no boolean operation
                 0,                        // copy number
                 fCheckOverlaps);          // checking overlaps

  //
  // Create towers (4x4 grid)
  //
  for (G4int iTowerX = 0; iTowerX < nofTowersX; iTowerX++) {
    for (G4int iTowerY = 0; iTowerY < nofTowersY; iTowerY++) {
      
      // Calculate tower position
      G4double towerPosX = (iTowerX - (nofTowersX-1)/2.0) * (towerSizeXY + towerSpacing);
      G4double towerPosY = (iTowerY - (nofTowersY-1)/2.0) * (towerSizeXY + towerSpacing);
      
      // Create tower logical volume
      auto towerS
        = new G4Box("Tower",           // its name
                     towerSizeXY/2, towerSizeXY/2, rodLength/2); // its size

      auto towerLV
        = new G4LogicalVolume(
                     towerS,           // its solid
                     defaultMaterial,  // its material
                     "Tower");         // its name

      new G4PVPlacement(
                     nullptr,  // no rotation
                     G4ThreeVector(towerPosX, towerPosY, 0), // position
                     towerLV,          // its logical volume
                     "Tower",          // its name
                     detectorLV,       // its mother  volume
                     false,            // no boolean operation
                     iTowerX * nofTowersY + iTowerY, // copy number
                     fCheckOverlaps);  // checking overlaps

      //
      // Create copper rods within tower (16x16 grid)
      //
      for (G4int iRodX = 0; iRodX < nofRodsPerTowerX; iRodX++) {
        for (G4int iRodY = 0; iRodY < nofRodsPerTowerY; iRodY++) {
          
          // Calculate rod position within tower
          G4double rodPosX = (iRodX - (nofRodsPerTowerX-1)/2.0) * (rodCrossSection + rodSpacing);
          G4double rodPosY = (iRodY - (nofRodsPerTowerY-1)/2.0) * (rodCrossSection + rodSpacing);
          
          // Create copper rod
          auto rodS
            = new G4Box("CopperRod",            // its name
                         rodCrossSection/2, rodCrossSection/2, rodLength/2); // its size

          auto rodLV
            = new G4LogicalVolume(
                         rodS,        // its solid
                         copperMaterial, // its material
                         "CopperRod");          // its name

          new G4PVPlacement(
                         nullptr,     // no rotation
                         G4ThreeVector(rodPosX, rodPosY, 0),  // its position
                         rodLV,                                // its logical volume
                         "CopperRod",                                    // its name
                         towerLV,                                   // its mother  volume
                         false,                                     // no boolean operation
                         iRodX * nofRodsPerTowerY + iRodY,         // copy number
                         fCheckOverlaps);                           // checking overlaps

          //
          // Create a vacuum hole LV inside the copper rod and place fibers into that hole
          //
          G4double fiberRadius = fiberDiameter/2;
          G4double holeRadius = centerHoleDiameter/2;

          // hole solid & logical volume (vacuum)
          auto holeS = new G4Tubs("RodHole", 0, holeRadius, rodLength/2 + 1.0*mm, 0, 360*degree);
          auto holeLV = new G4LogicalVolume(holeS, defaultMaterial, "RodHoleLV");
          holeLV->SetVisAttributes(G4VisAttributes::GetInvisible());
          new G4PVPlacement(nullptr, G4ThreeVector(), holeLV, "RodHolePV", rodLV, false, 0, fCheckOverlaps);

          // compute safe radius for fiber placement (outer ring + clearance)
          G4double safety = 0.05*mm; // small clearance from wall
          G4double outerRingR = holeRadius - fiberRadius - safety;
          if (outerRingR < 0) outerRingR = 0.0;

          // center quartz fiber
          auto fiberQuartzCenterS
            = new G4Tubs("QuartzFiberCenter", 0, fiberRadius, rodLength/2, 0, 360*degree);
          auto fiberQuartzCenterLV
            = new G4LogicalVolume(fiberQuartzCenterS, quartzMaterial, "QuartzFiberCenter");
          new G4PVPlacement(nullptr, G4ThreeVector(), fiberQuartzCenterLV, "QuartzFiberCenter",
                            holeLV, false, 0, fCheckOverlaps);

          // outer ring: alternate quartz and scintillator every 60 degrees
          for (G4int iFiber = 0; iFiber < 6; iFiber++) {
            G4double angle = (iFiber * 60.0) * degree; // 60° spacing
            G4double fiberPosX = outerRingR * std::cos(angle);
            G4double fiberPosY = outerRingR * std::sin(angle);

            G4bool isQuartz = (iFiber % 2 == 0); // even: quartz, odd: scintillator
            auto fiberS = new G4Tubs(isQuartz ? "QuartzFiber" : "ScintFiber",
                                     0, fiberRadius, rodLength/2, 0, 360*degree);
            auto fiberLV = new G4LogicalVolume(fiberS,
                                               isQuartz ? quartzMaterial : scintillatorMaterial,
                                               isQuartz ? "QuartzFiber" : "ScintFiber");

            new G4PVPlacement(nullptr, G4ThreeVector(fiberPosX, fiberPosY, 0),
                              fiberLV,
                              isQuartz ? "QuartzFiber" : "ScintFiber",
                              holeLV,
                              false,
                              iFiber + 1, // +1 to distinguish from center fiber
                              fCheckOverlaps);
          }
        }
      }
    }
  }

  //
  // Print parameters
  //
  G4cout
    << G4endl
    << "------------------------------------------------------------" << G4endl
    << "---> DREAM Detector Configuration:" << G4endl
    << "---> Towers: " << nofTowersX << "x" << nofTowersY << G4endl
    << "---> Rods per tower: " << nofRodsPerTowerX << "x" << nofRodsPerTowerY << G4endl
    << "---> Copper rod size: " << rodCrossSection/mm << "mm x " << rodCrossSection/mm
    << " x " << rodLength/m << "m" << G4endl
    << "---> Center hole diameter: " << centerHoleDiameter/mm << "mm" << G4endl
    << "---> Fiber diameter: " << fiberDiameter/mm << "mm" << G4endl
    << "---> Detector rotation: +2° around Y, +0.7° around X" << G4endl
    << "------------------------------------------------------------" << G4endl;

  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto detectorVisAtt = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.3));
  detectorLV->SetVisAttributes(detectorVisAtt);

  auto towerVisAtt = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.5));
  auto copperVisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.2)); // Copper color
  auto scintVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // Green for scintillator
  auto quartzVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // Blue for quartz

  // Apply visualization attributes using string comparison instead of contains()
  G4LogicalVolumeStore* logVolStore = G4LogicalVolumeStore::GetInstance();
  for (auto logVol : *logVolStore) {
    G4String volName = logVol->GetName();
    if (volName.find("Tower") != G4String::npos) {
      logVol->SetVisAttributes(towerVisAtt);
    } else if (volName.find("CopperRod") != G4String::npos) {
      logVol->SetVisAttributes(copperVisAtt);
    } else if (volName.find("ScintFiber") != G4String::npos) {
      logVol->SetVisAttributes(scintVisAtt);
    } else if (volName.find("QuartzFiber") != G4String::npos) {
      logVol->SetVisAttributes(quartzVisAtt);
    }
  }

  //
  // Always return the physical World
  //
  return worldPV;
}

void DetectorConstruction::ConstructSDandField()
{
  // global magnetic field (kept as before)
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  G4AutoDelete::Register(fMagFieldMessenger);

  // --- Sensitive detectors for fibers ---
  auto sdMan = G4SDManager::GetSDMpointer();
  auto scintSD  = new FiberSD("ScintFiberSD", "ScintFiberHitsCollection", 0);  // 0 表示闪烁光纤
  auto quartzSD = new FiberSD("QuartzFiberSD", "QuartzFiberHitsCollection", 1); // 1 表示切伦科夫光纤

  sdMan->AddNewDetector(scintSD);
  sdMan->AddNewDetector(quartzSD);

  // Attach SDs to logical volumes by name
  G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  for (auto lv : *lvStore) {
    G4String name = lv->GetName();
    if (name.find("ScintFiber") != G4String::npos) {
      lv->SetSensitiveDetector(scintSD);
    } else if (name.find("QuartzFiber") != G4String::npos
               || name.find("CkovFiber") != G4String::npos
               || name.find("Ckov") != G4String::npos) {
      lv->SetSensitiveDetector(quartzSD);
    }
  }
}
}
// ConstructSDandField: 挂接全局场（SD 在后续实现）