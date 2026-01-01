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
//
/// \file B4/B4a/src/EventAction.cc
/// \brief Implementation of the B4a::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"  
#include "FiberHit.hh"
#include "G4ParticleGun.hh" 

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include <iomanip>
#include "CLHEP/Units/SystemOfUnits.h"

namespace B4a
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
  // 初始化每事件变量（保留兼容性，实际不再使用）
  fEnergyAbs = 0.;
  fEnergyGap = 0.;
  fTrackLAbs = 0.;
  fTrackLGap = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  // ===== 获取 HitsCollection =====
  auto hce = event->GetHCofThisEvent();
  if(!hce) {
    G4cout << "Warning: No HCofThisEvent in event " << event->GetEventID() << G4endl;
    return;
  }
  
  auto sdManager = G4SDManager::GetSDMpointer();
  G4int scintHCID = sdManager->GetCollectionID("ScintFiberHitsCollection");
  G4int quartzHCID = sdManager->GetCollectionID("QuartzFiberHitsCollection");
  
  // 首次运行时 ID 可能为 -1，获取后会缓存
  if(scintHCID < 0 || quartzHCID < 0) {
    if(event->GetEventID() == 0) {
      G4cout << "Warning: HitsCollection ID not found in first event. "
             << "ScintID=" << scintHCID << ", QuartzID=" << quartzHCID << G4endl;
    }
    return;
  }
  
  auto scintHC = static_cast<FiberHitsCollection*>(hce->GetHC(scintHCID));
  auto quartzHC = static_cast<FiberHitsCollection*>(hce->GetHC(quartzHCID));
  
  // ===== 聚合光子数（所有 tower 的总和）=====
  G4int totalNs = 0;  // 闪烁光子总数
  G4int totalNc = 0;  // 切伦科夫光子总数
  
  if(scintHC) {
    for(size_t i = 0; i < scintHC->entries(); i++) {
      totalNs += (*scintHC)[i]->GetNsPhotons();
    }
  }
  
  if(quartzHC) {
    for(size_t i = 0; i < quartzHC->entries(); i++) {
      totalNc += (*quartzHC)[i]->GetNcPhotons();
    }
  }
  
  // ===== 获取 RunAction 并填充数据到 ROOT =====
  auto runAction = const_cast<RunAction*>(
    static_cast<const RunAction*>(
      G4RunManager::GetRunManager()->GetUserRunAction()));
  
  if(runAction) {
    // 获取 PrimaryGeneratorAction
    auto generatorAction = const_cast<B4::PrimaryGeneratorAction*>(
      static_cast<const B4::PrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction()));
    
    G4double inputEnergy = 0; // 默认值
    if (generatorAction) {
      auto particleGun = generatorAction->GetParticleGun();
      if (particleGun) {
        inputEnergy = particleGun->GetParticleEnergy() / CLHEP::GeV;
      }
    }

    runAction->FillEvent(event->GetEventID(), totalNs, totalNc, inputEnergy);
  }
  
  // ===== 定期打印调试信息 =====
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  
  if((printModulo > 0) && (eventID % printModulo == 0)) {
    G4cout << "Event " << eventID 
           << ": Ns = " << totalNs 
           << ", Nc = " << totalNc << G4endl;
    G4cout << "--> End of event " << eventID << "\n" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}