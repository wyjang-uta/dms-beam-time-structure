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
/// \file DMSEventAction.hh
/// \brief Definition of the DMSEventAction class

#ifndef DMSEventAction_h
#define DMSEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class DMSRunAction;

/// Event action class
///

class DMSEventAction : public G4UserEventAction
{
  public:
    DMSEventAction(DMSRunAction* runAction);
    virtual ~DMSEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

  private:
    DMSRunAction* fRunAction;
    G4int currentBin;
    G4double tlow[620];
    G4double thigh[620];
    G4double time[620];
    G4int bin[620];
    G4int binLimit[620];
    /*
    G4double fBunchDuration;
    G4double fBunchInterval;
    G4double fParticlesPerBunch;
    G4int fNumberOfBunches;
    */
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
