#ifndef FiberSD_h
#define FiberSD_h 1

#include "G4VSensitiveDetector.hh"
#include "FiberHit.hh"
#include "globals.hh"

class G4Step;
class G4HCofThisEvent;

using FiberHitsCollection = G4THitsCollection<FiberHit>;

class FiberSD : public G4VSensitiveDetector
{
public:
  FiberSD(const G4String& name, const G4String& hitsCollectionName, G4int fiberType);
  virtual ~FiberSD();

  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*);

private:
  FiberHitsCollection* fHitsCollection;
  G4int fFiberType;  // 0: Scintillation, 1: Quartz/Cherenkov
  
  // 光学参数
  G4double fScintillationYield;
  G4double fCollectionEfficiency;
  G4double fRefractiveIndex;
};

#endif
