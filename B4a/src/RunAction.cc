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
/// \file B4/B4a/src/RunAction.cc
/// \brief Implementation of the B4::RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"  
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh" 
// ROOT 头文件
#include "TFile.h"
#include "TTree.h"
#include "CLHEP/Units/SystemOfUnits.h" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction(),
   fRootFile(nullptr),
   fEventTree(nullptr),
   fEventID(0),
   fNsPhotons(0),
   fNcPhotons(0),
   fInputEnergy(0.)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  // 获取粒子枪能量
  auto generatorAction = const_cast<B4::PrimaryGeneratorAction*>(
    static_cast<const B4::PrimaryGeneratorAction*>(
      G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction()
    )
  );
  G4double energy = 0.;
  if (generatorAction) {
    auto particleGun = generatorAction->GetParticleGun();
    if (particleGun) {
      energy = particleGun->GetParticleEnergy() / CLHEP::GeV; 
    }
  }

  // 创建带能量标识的文件名
  G4int runID = run->GetRunID();
  G4String fileName = "dream_pi-_" + std::to_string((int)energy) 
                      + "GeV_run" + std::to_string(runID) + ".root";
  
  fRootFile = new TFile(fileName.c_str(), "RECREATE");
  if (!fRootFile || fRootFile->IsZombie()) {
    G4cerr << "ERROR: Cannot create ROOT file " << fileName << G4endl;
    return;
  }
  
  G4cout << "ROOT file created: " << fileName << G4endl;

  // 创建 TTree
  fEventTree = new TTree("dreamTree", "DREAM Detector Event Data");
  fEventTree->Branch("EventID", &fEventID, "EventID/I");
  fEventTree->Branch("NsPhotons", &fNsPhotons, "NsPhotons/I");
  fEventTree->Branch("NcPhotons", &fNcPhotons, "NcPhotons/I");
  fEventTree->Branch("InputEnergy", &fInputEnergy, "InputEnergy/D");
  
  G4cout << "TTree 'dreamTree' created with branches:" << G4endl;
  G4cout << "  - EventID (Int)" << G4endl;
  G4cout << "  - NsPhotons (Int)" << G4endl;
  G4cout << "  - NcPhotons (Int)" << G4endl;
  G4cout << "  - InputEnergy (Double)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) {
    G4cout << "### Run " << run->GetRunID() << " ended (no events)." << G4endl;
    return;
  }

  // 写入并关闭 ROOT 文件
  if (fRootFile) {
    fRootFile->cd();
    
    if (fEventTree) {
      G4int entries = fEventTree->GetEntries();
      fEventTree->Write();
      G4cout << "TTree 'dreamTree' written with " << entries << " entries." << G4endl;
    }
    
    fRootFile->Close();
    G4cout << "ROOT file closed." << G4endl;
    
    delete fRootFile;
    fRootFile = nullptr;
    fEventTree = nullptr;
  }

  G4cout << "### Run " << run->GetRunID() << " ended." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillEvent(G4int eventID, G4int ns, G4int nc, G4double energy)
{
  if (!fEventTree) {
    G4cerr << "ERROR: TTree not initialized!" << G4endl;
    return;
  }

  // 填充分支数据
  fEventID = eventID;
  fNsPhotons = ns;
  fNcPhotons = nc;
  fInputEnergy = energy;

  // 写入 TTree
  fEventTree->Fill();
  
  // 每10个事件打印一次
  if (eventID % 10 == 0) {
    G4cout << "Event " << eventID << " filled: "
           << "Ns=" << ns << ", Nc=" << nc 
           << ", E=" << energy << " GeV" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......