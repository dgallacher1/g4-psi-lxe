#include "PhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"   // most accurate std EM (e+/e-/gamma/ion ionisation, msc)
#include "G4DecayPhysics.hh"                // ordinary decays
#include "G4RadioactiveDecayPhysics.hh"     // needed if the source is a radionuclide (e.g. Am-241)
#include "G4EmParameters.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList()
{
  SetVerboseLevel(1);

  // e+ transport, annihilation, bremsstrahlung, and ion (alpha) ionisation
  // are all handled by the standard EM list. option4 is the high-accuracy
  // variant (Goudsmit-Saunderson msc for e-, tighter step limits).
  RegisterPhysics(new G4EmStandardPhysics_option4());

  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  // Atomic de-excitation. option4 already has fluorescence on, 3-gamma annihilation
  // in flight and Allison positronium at rest. Add Auger electrons (and the cascade)
  // and make them ignore the production cuts, so Xe K/L vacancies are followed
  // explicitly in the LXe and in the foils. The total energy deposited does not
  // change (without Auger the non-radiative energy is deposited locally); what changes
  // is where it goes near boundaries and the number of sub-keV electrons.
  // These calls come after RegisterPhysics(option4) because its constructor
  // calls G4EmParameters::SetDefaults() and sets its own values.
  {
    auto* em = G4EmParameters::Instance();
    em->SetFluo(true);
    em->SetAuger(true);                      // 11.4: also switches on the Auger cascade
    em->SetDeexcitationIgnoreCut(true);
    em->SetDeexActiveRegion("DefaultRegionForTheWorld", true, true, false);   // LXe, steel
    em->SetDeexActiveRegion("WindowRegion",             true, true, false);   // Al, Ti foils
  }

  // No hadronic list: 5 MeV-scale e+/alpha have negligible nuclear
  // interaction. Add G4IonPhysics + a hadronic list only if you push
  // to higher energies or care about (alpha,n) etc.

  // Range cut. Lower it inside the LXe if you need finer track structure
  // near the source (see /run/setCutForRegion, or a G4Region on the LXe LV).
  SetDefaultCutValue(0.1 * mm);
}
