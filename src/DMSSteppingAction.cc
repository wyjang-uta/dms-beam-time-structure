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
/// \file DMSSteppingAction.cc
/// \brief Implementation of the DMSSteppingAction class

#include "DMSSteppingAction.hh"
#include "DMSEventAction.hh"
#include "DMSDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleTypes.hh"
#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DMSSteppingAction::DMSSteppingAction()
: G4UserSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DMSSteppingAction::~DMSSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DMSSteppingAction::UserSteppingAction(const G4Step* step)
{
  auto analysisManager = G4AnalysisManager::Instance();

  G4StepPoint* presp = step->GetPreStepPoint();
  G4StepPoint* postsp = step->GetPostStepPoint();

  G4int evt_type = 0;

  // This routine records energy, momentum, and hit position of neutrons that are coming into the vacuum chamber.
  // In order to identify neutrons coming into the chamber, we compare the step volume of the post-step point and the pre-step point.
  // If the post-step point is "sensor_volume" but the pre-step point is some other volume, we accept those neutrons.
  //
  // To do: is there double counting issue?
  // I think we can prevent this by rejecting the case that the pre-step volume and the post-step volume both are the 'sensor_volume'.
  //
  if( presp == nullptr ) return;
  if( postsp == nullptr ) return;
  if( step->GetTrack()->GetParticleDefinition()->GetPDGEncoding() != 2112 ) return; // only neutrons
  if( postsp->GetTouchableHandle()->GetVolume() == nullptr ) return;
  if( presp->GetTouchableHandle()->GetVolume() == nullptr ) return;
  if( postsp->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() != "sensor_volume" ) return; // save only neutrons passing decay chamber
  if(
      postsp->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "sensor_volume" &&
      presp->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "sensor_volume" )
    return; // we reject this case to prevent double counting issue.

  if( presp->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "Dump")
  {
    evt_type = 1;  // This is for the case that neutron comes from the dump.
  }
  if( presp->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "World" )
  {
    evt_type = 2; // This is for the case that neutrons are coming from the world volume.
  }

  G4double KE=(G4double)postsp->GetKineticEnergy()/CLHEP::MeV;    // Kinetic energy
  G4double TE=(G4double)postsp->GetTotalEnergy()/CLHEP::MeV;      // Total energy
  G4double vx=(G4double)postsp->GetPosition().getX()/CLHEP::m;    // Position where the neutron was produced.
  G4double vy=(G4double)postsp->GetPosition().getY()/CLHEP::m;
  G4double vz=(G4double)postsp->GetPosition().getZ()/CLHEP::m;
  G4double px=(G4double)postsp->GetMomentum().getX()/CLHEP::MeV;  // Momentum of the neutron.
  G4double py=(G4double)postsp->GetMomentum().getY()/CLHEP::MeV;
  G4double pz=(G4double)postsp->GetMomentum().getZ()/CLHEP::MeV;
  G4double gt=(G4double)postsp->GetGlobalTime()/CLHEP::ns;
  G4double lt=(G4double)postsp->GetLocalTime()/CLHEP::ns;
  //G4double theta=acos(pz/sqrt(px*px+py*py+pz*pz));                // Direction of the neutron with respect to the z-axis (beam direction).

  // event type
  analysisManager->FillNtupleDColumn(0, (G4double)evt_type );
  // kineticEnergy
  analysisManager->FillNtupleDColumn(1, KE);
  // total energy
  analysisManager->FillNtupleDColumn(2, TE);
  // vertex coordinates
  analysisManager->FillNtupleDColumn(3, vx);
  analysisManager->FillNtupleDColumn(4, vy);
  analysisManager->FillNtupleDColumn(5, vz);
  // momentum
  analysisManager->FillNtupleDColumn(6, px);
  analysisManager->FillNtupleDColumn(7, py);
  analysisManager->FillNtupleDColumn(8, pz);
  // time info
  analysisManager->FillNtupleDColumn(9, gt);
  analysisManager->FillNtupleDColumn(10, lt);
  analysisManager->AddNtupleRow();






  /*
  G4StepPoint* postStep = step->GetPostStepPoint();
  G4StepPoint* preStep = step->GetPreStepPoint();
  G4Track* track = step->GetTrack();
  if( !postStep ) return;
  if( !preStep  ) return;
  if( !track ) return;
  G4VPhysicalVolume* postVol = postStep->GetPhysicalVolume();
  G4VPhysicalVolume* preVol = postStep->GetPhysicalVolume();
  if( !postVol ) return;
  if( !preVol  ) return;
  if( postStep->GetPhysicalVolume()->GetName() != "physDetector" ) return;
  if( preStep->GetPhysicalVolume()->GetName() == "physDetector" ) return; // deny double counting.

  if( track->GetDefinition()->GetParticleName() == "nu_e" ) return;
  if( track->GetDefinition()->GetParticleName() == "anti_nu_e" ) return;
  if( track->GetDefinition()->GetParticleName() == "nu_mu" ) return;
  if( track->GetDefinition()->GetParticleName() == "anti_nu_mu" ) return;

  analysisManager->FillNtupleSColumn(0, track->GetDefinition()->GetParticleName() );
  analysisManager->FillNtupleSColumn(1, track->GetCreatorProcess()->GetProcessName() );
  analysisManager->FillNtupleDColumn(2, track->GetParentID());
  analysisManager->FillNtupleDColumn(3, track->GetTrackID());
  analysisManager->FillNtupleDColumn(4, track->GetCurrentStepNumber());
  analysisManager->FillNtupleDColumn(5, (G4double)track->GetKineticEnergy()/CLHEP::MeV);
  analysisManager->FillNtupleDColumn(6, (G4double)track->GetTotalEnergy()/CLHEP::MeV);
  analysisManager->FillNtupleDColumn(7, (G4double)track->GetPosition().getX()/CLHEP::cm);
  analysisManager->FillNtupleDColumn(8, (G4double)track->GetPosition().getY()/CLHEP::cm);
  analysisManager->FillNtupleDColumn(9, (G4double)track->GetPosition().getZ()/CLHEP::cm);
  analysisManager->FillNtupleDColumn(10, (G4double)track->GetMomentum().getX()/CLHEP::MeV);
  analysisManager->FillNtupleDColumn(11, (G4double)track->GetMomentum().getY()/CLHEP::MeV);
  analysisManager->FillNtupleDColumn(12, (G4double)track->GetMomentum().getZ()/CLHEP::MeV);
  analysisManager->FillNtupleDColumn(13, (G4double)track->GetGlobalTime()/CLHEP::ns);
  analysisManager->FillNtupleDColumn(14, (G4double)track->GetLocalTime()/CLHEP::ns);
  analysisManager->AddNtupleRow();
  */

  /*
  const std::vector<const G4Track*>* secondary = step->GetSecondaryInCurrentStep();
  for( size_t lp=0; lp < (*secondary).size(); lp++ )
  {
    // secondary particle name
    analysisManager->FillNtupleSColumn(0, (*secondary)[lp]->GetDefinition()->GetParticleName() );
    // production process of secondary particle
    analysisManager->FillNtupleSColumn(1, (*secondary)[lp]->GetCreatorProcess()->GetProcessName());
    // parent id
    analysisManager->FillNtupleDColumn(2, (*secondary)[lp]->GetParentID());
    // track id
    analysisManager->FillNtupleDColumn(3, (*secondary)[lp]->GetTrackID());
    // step number
    analysisManager->FillNtupleDColumn(4, (*secondary)[lp]->GetCurrentStepNumber());
    // kineticEnergy
    analysisManager->FillNtupleDColumn(5, (G4double)(*secondary)[lp]->GetKineticEnergy()/CLHEP::MeV);
    analysisManager->FillNtupleDColumn(6, (G4double)(*secondary)[lp]->GetTotalEnergy()/CLHEP::MeV);
    // vertex coordinates
    analysisManager->FillNtupleDColumn(7, (G4double)(*secondary)[lp]->GetPosition().getX()/CLHEP::cm);
    analysisManager->FillNtupleDColumn(8, (G4double)(*secondary)[lp]->GetPosition().getY()/CLHEP::cm);
    analysisManager->FillNtupleDColumn(9, (G4double)(*secondary)[lp]->GetPosition().getZ()/CLHEP::cm);
    // momentum
    analysisManager->FillNtupleDColumn(10, (G4double)(*secondary)[lp]->GetMomentum().getX()/CLHEP::MeV);
    analysisManager->FillNtupleDColumn(11, (G4double)(*secondary)[lp]->GetMomentum().getY()/CLHEP::MeV);
    analysisManager->FillNtupleDColumn(12, (G4double)(*secondary)[lp]->GetMomentum().getZ()/CLHEP::MeV);
    // time info
    analysisManager->FillNtupleDColumn(13, (G4double)(*secondary)[lp]->GetGlobalTime()/CLHEP::ns);
    analysisManager->FillNtupleDColumn(14, (G4double)(*secondary)[lp]->GetLocalTime()/CLHEP::ns);
    analysisManager->AddNtupleRow();
  }
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

