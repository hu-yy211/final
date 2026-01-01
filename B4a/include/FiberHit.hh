#ifndef FiberHit_h
#define FiberHit_h 1

#include "G4VHit.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "globals.hh"

class FiberHit : public G4VHit
{
public:
  FiberHit();
  ~FiberHit() override;
  FiberHit(const FiberHit&);
  const FiberHit& operator=(const FiberHit&);
  G4int operator==(const FiberHit&) const;

  inline void* operator new(size_t);
  inline void operator delete(void*);

  // setters
  void SetPosition(const G4ThreeVector& pos) { fPosition = pos; }
  void SetEdep(G4double e)                  { fEnergyDeposit = e; }
  void SetNsPhotons(G4int n)                { fNsPhotons = n; }
  void SetNcPhotons(G4int n)                { fNcPhotons = n; }
  void AddNsPhotons(G4int n)                { fNsPhotons += n; }
  void AddNcPhotons(G4int n)                { fNcPhotons += n; }
  void SetFiberType(G4int t)                { fFiberType = t; }
  void SetCopyID(G4int id)                  { fCopyNumber = id; }

  // getters
  const G4ThreeVector& GetPosition() const { return fPosition; }
  G4double GetEdep()   const               { return fEnergyDeposit; }
  G4int GetNsPhotons() const               { return fNsPhotons; }
  G4int GetNcPhotons() const               { return fNcPhotons; }
  G4int GetFiberType() const               { return fFiberType; }
  G4int GetCopyID() const                  { return fCopyNumber; }

private:
  G4ThreeVector fPosition;
  G4double      fEnergyDeposit;
  G4int         fNsPhotons;   // 闪烁光子数
  G4int         fNcPhotons;   // 切伦科夫光子数
  G4int         fFiberType;
  G4int         fCopyNumber;
};

using FiberHitsCollection = G4THitsCollection<FiberHit>;

extern G4ThreadLocal G4Allocator<FiberHit>* FiberHitAllocator;

inline void* FiberHit::operator new(size_t)
{
  if(!FiberHitAllocator) FiberHitAllocator = new G4Allocator<FiberHit>;
  return (void*)FiberHitAllocator->MallocSingle();
}

inline void FiberHit::operator delete(void* hit)
{
  FiberHitAllocator->FreeSingle((FiberHit*)hit);
}

#endif