#include "FiberSD.hh"
#include "FiberHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include <cmath>

// CLHEP Poisson 与单位
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Units/SystemOfUnits.h"

FiberSD::FiberSD(const G4String& name,
                 const G4String& hitsCollectionName,
                 G4int fiberType)
  : G4VSensitiveDetector(name),
    fFiberType(fiberType),  // 0: Scint, 1: Quartz
    fHitsCollection(nullptr)
{
  collectionName.insert(hitsCollectionName);

  // 设置光学参数
  if (fFiberType == 0) {         // 塑料闪烁光纤
    fScintillationYield   = 10000.0 / CLHEP::MeV; // 文档要求
    fCollectionEfficiency = 0.9;                  // 文档要求
  } else {                       // 石英切伦科夫光纤
    fRefractiveIndex      = 1.474;                // 示例值，可按实际介质调整
    fCollectionEfficiency = 0.9;                  // 文档要求
  }
}

FiberSD::~FiberSD() {}

void FiberSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new FiberHitsCollection(SensitiveDetectorName, collectionName[0]);

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

G4bool FiberSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // 共用信息
  G4double edep       = aStep->GetTotalEnergyDeposit();
  G4double stepLength = aStep->GetStepLength();
  G4Track* track      = aStep->GetTrack();
  G4ParticleDefinition* particle = track->GetDefinition();
  G4double charge     = particle->GetPDGCharge();

  G4int nPhotons = 0;

  if (fFiberType == 0) {
    // ===== 闪烁光纤：需有能量沉积 =====
    if (edep <= 0.) return false;
    G4double meanPhotons = edep * fScintillationYield * fCollectionEfficiency;
    nPhotons = (meanPhotons > 0.) ? CLHEP::RandPoisson::shoot(meanPhotons) : 0;
  } else {
    // ===== 石英光纤：切伦科夫光 =====
    // 只有带电粒子
    if (std::abs(charge) < 0.1) return false;

    // 阈值判定
    G4double beta = aStep->GetPreStepPoint()->GetBeta();
    G4double betaThreshold = 1.0 / fRefractiveIndex;
    if (beta <= betaThreshold) return false;

    // 产生率 dN/dx（每厘米）
    if (stepLength <= 0.) return false;
    G4double dNdL = 369.0 * (1.0 - 1.0 / (beta * beta * fRefractiveIndex * fRefractiveIndex));

    // 总期望值并乘收集效率
    G4double meanPhotons = dNdL * (stepLength / CLHEP::cm) * fCollectionEfficiency;
    nPhotons = (meanPhotons > 0.) ? CLHEP::RandPoisson::shoot(meanPhotons) : 0;
  }

  if (nPhotons <= 0) return false;

  // ===== 创建 Hit 对象 =====
  auto* newHit = new FiberHit();

  G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
  newHit->SetPosition(pos);
  newHit->SetEdep(edep);
  
  if (fFiberType == 0) {        // 闪烁
    newHit->SetNsPhotons(nPhotons);
    newHit->SetNcPhotons(0);
  } else {                      // 切伦科夫
    newHit->SetNsPhotons(0);
    newHit->SetNcPhotons(nPhotons);
  }
  
  newHit->SetFiberType(fFiberType);

  // 副本号
  const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int copyNo = touchable->GetCopyNumber();
  newHit->SetCopyID(copyNo);

  fHitsCollection->insert(newHit);
  return true;
}

void FiberSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 0) {
    G4cout << "FiberSD::EndOfEvent - " << fHitsCollection->GetSize()
           << " hits in collection" << G4endl;
  }
}