#include "FiberHit.hh"

G4ThreadLocal G4Allocator<FiberHit>* FiberHitAllocator = nullptr;

FiberHit::FiberHit()
  : G4VHit(),
    fPosition(0.,0.,0.),
    fEnergyDeposit(0.),
    fNsPhotons(0),
    fNcPhotons(0),
    fFiberType(0),
    fCopyNumber(0)
{}

FiberHit::~FiberHit() = default;

FiberHit::FiberHit(const FiberHit& right)
  : G4VHit(),
    fPosition(right.fPosition),
    fEnergyDeposit(right.fEnergyDeposit),
    fNsPhotons(right.fNsPhotons),
    fNcPhotons(right.fNcPhotons),
    fFiberType(right.fFiberType),
    fCopyNumber(right.fCopyNumber)
{}

const FiberHit& FiberHit::operator=(const FiberHit& right)
{
  if (this != &right) {
    fPosition      = right.fPosition;
    fEnergyDeposit = right.fEnergyDeposit;
    fNsPhotons     = right.fNsPhotons;
    fNcPhotons     = right.fNcPhotons;
    fFiberType     = right.fFiberType;
    fCopyNumber    = right.fCopyNumber;
  }
  return *this;
}

G4int FiberHit::operator==(const FiberHit& right) const
{
  return (this == &right) ? 1 : 0;
}