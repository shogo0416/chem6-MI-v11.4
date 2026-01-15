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
// Geant4 class G4ChemDissociationChannels_option1
//
// H. Tran 16.12.2022
//

#include "G4ChemDissociationChannels_option1.hh"

#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4Electron_aq.hh"
#include "G4FakeMolecule.hh"
#include "G4H2.hh"
#include "G4H2O.hh"
#include "G4H2O2.hh"
#include "G4H3O.hh"
#include "G4HO2.hh"
#include "G4Hydrogen.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeTable.hh"
#include "G4O2.hh"
#include "G4O3.hh"
#include "G4OH.hh"
#include "G4Oxygen.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scheduler.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4ChemDissociationChannels_option1::ConstructMolecule()
{
  //-----------------------------------
  //  G4Electron::Definition(); // safety

  //-----------------------------------
  // Create the definition

  G4H2O::Definition();
  G4Hydrogen::Definition();
  G4H3O::Definition();
  G4OH::Definition();
  G4Electron_aq::Definition();
  G4H2O2::Definition();
  G4H2::Definition();

  G4O2::Definition();
  G4HO2::Definition();
  G4Oxygen::Definition();
  G4O3::Definition();

  auto G4OHm = new G4MoleculeDefinition("OH",/*mass*/ 17.00734 * g / Avogadro * c_squared,
                                        2.8e-9 * (m * m / s), -1,
                                          5, 0.958 * angstrom, // radius
                                          2 // number of atoms
                                        );

  auto G4HO2m = new G4MoleculeDefinition("HO_2", 33.0034 * g / Avogadro * c_squared,
                                         2.3e-9 * (m * m / s), -1, 0,
                                        2.1 * angstrom, 3);
  auto G4Om = new G4MoleculeDefinition("O", 15.99773 * g / Avogadro * c_squared,
                                       2.0e-9 * (m * m / s), 0, 0,
                                       2.0 * angstrom, 1);
  //____________________________________________________________________________
  auto molTable = G4MoleculeTable::Instance();
  molTable->CreateConfiguration("H3Op", G4H3O::Definition());
  molTable->GetConfiguration("H3Op")->SetDiffusionCoefficient(9.46e-9
                                                                                 * (m2 / s));
  molTable->GetConfiguration("H3Op")->SetVanDerVaalsRadius(0.25 * nm);

  molTable->CreateConfiguration("°OH", G4OH::Definition());
  molTable->GetConfiguration("°OH")->SetDiffusionCoefficient(2.2e-9 * (m2 / s));
  molTable->GetConfiguration("°OH")->SetVanDerVaalsRadius(0.22 * nm);

  G4MolecularConfiguration* OHm =
    molTable->CreateConfiguration("OHm",  // just a tag to store and retrieve
                                                             // from G4MoleculeTable
                                  G4OHm,
                                                     -1,  // charge
                                                     5.3e-9 * (m2 / s));
  OHm->SetMass(17.0079 * g / Avogadro * c_squared);
  OHm->SetVanDerVaalsRadius(0.33 * nm);

  molTable->CreateConfiguration("e_aq", G4Electron_aq::Definition());
  molTable->GetConfiguration("e_aq")->SetVanDerVaalsRadius(0.50 * nm);

  molTable->CreateConfiguration("H", G4Hydrogen::Definition());
  molTable->GetConfiguration("H")->SetVanDerVaalsRadius(0.19 * nm);

  molTable->CreateConfiguration("H2", G4H2::Definition());
  molTable->GetConfiguration("H2")->SetDiffusionCoefficient(4.8e-9 * (m2 / s));
  molTable->GetConfiguration("H2")->SetVanDerVaalsRadius(0.14 * nm);

  molTable->CreateConfiguration("H2O2", G4H2O2::Definition());
  molTable->GetConfiguration("H2O2")->SetDiffusionCoefficient(2.3e-9 * (m2 / s));
  molTable->GetConfiguration("H2O2")->SetVanDerVaalsRadius(0.21 * nm);

  // molecules extension (RITRACKS)

  molTable->CreateConfiguration("HO2°", G4HO2::Definition());
  molTable->GetConfiguration("HO2°")->SetVanDerVaalsRadius(0.21 * nm);

  G4MolecularConfiguration* HO2m =
    molTable->CreateConfiguration("HO2m",  // just a tag to store and retrieve
                                                              // from G4MoleculeTable
                                                        G4HO2m,
                                                     -1,  // charge
                                                     1.4e-9 * (m2 / s));
  HO2m->SetMass(33.00396 * g / Avogadro * c_squared);
  HO2m->SetVanDerVaalsRadius(0.25 * nm);

  molTable->CreateConfiguration("Oxy", G4Oxygen::Definition());
  molTable->GetConfiguration("Oxy")->SetVanDerVaalsRadius(0.20 * nm);

  G4MolecularConfiguration* Om =
    molTable->CreateConfiguration("Om",  // just a tag to store and retrieve from
                                                            // G4MoleculeTable
                                                      G4Om,
                                                     -1,  // charge
                                                     2.0e-9 * (m2 / s));
  Om->SetMass(15.99829 * g / Avogadro * c_squared);
  Om->SetVanDerVaalsRadius(0.25 * nm);

  molTable->CreateConfiguration("O2", G4O2::Definition());
  molTable->GetConfiguration("O2")->SetVanDerVaalsRadius(0.17 * nm);

  G4MolecularConfiguration* O2m =
    molTable->CreateConfiguration("O2m",  // just a tag to store and retrieve
                                                             // from G4MoleculeTable
                                                     G4O2::Definition(),
                                                     -1,  // charge
                                                     1.75e-9 * (m2 / s));
  O2m->SetMass(31.99602 * g / Avogadro * c_squared);
  O2m->SetVanDerVaalsRadius(0.22 * nm);

  molTable->CreateConfiguration("O3", G4O3::Definition());
  molTable->GetConfiguration("O3")->SetVanDerVaalsRadius(0.20 * nm);

  G4MolecularConfiguration* O3m =
    molTable->CreateConfiguration("O3m",  // just a tag to store and retrieve
                                                             // from G4MoleculeTable
                                                     G4O3::Definition(),
                                                     -1,  // charge
                                                     2.0e-9 * (m2 / s));
  O3m->SetMass(47.99375 * g / Avogadro * c_squared);
  O3m->SetVanDerVaalsRadius(0.20 * nm);

  molTable->CreateConfiguration("H2O(B)",  // just a tag to store and retrieve
                                                              // from G4MoleculeTable
                                                   G4H2O::Definition(),
                                                   0,  // charge
                                                   0 * (m2 / s));

  molTable->CreateConfiguration("H3Op(B)",  // just a tag to store and retrieve
                                                               // from G4MoleculeTable
                                                   G4H3O::Definition(),
                                                   1,  // charge
                                                   0 * (m2 / s));

  molTable->CreateConfiguration("OHm(B)",  // just a tag to store and retrieve
                                                              // from G4MoleculeTable
                                                    G4OHm,
                                                   -1,  // charge
                                                   0 * (m2 / s));

  molTable->CreateConfiguration("NoneM", G4FakeMolecule::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChemDissociationChannels_option1::ConstructDissociationChannels()
{
  //-----------------------------------
  // Get the molecular configuration
  auto molTable = G4MoleculeTable::Instance();
  G4MolecularConfiguration* OH = molTable->GetConfiguration("°OH");
  G4MolecularConfiguration* OHm = molTable->GetConfiguration("OHm");
  G4MolecularConfiguration* e_aq = molTable->GetConfiguration("e_aq");
  G4MolecularConfiguration* H2 = molTable->GetConfiguration("H2");
  G4MolecularConfiguration* H3O = molTable->GetConfiguration("H3Op");
  G4MolecularConfiguration* H = molTable->GetConfiguration("H");
  G4MolecularConfiguration* O = molTable->GetConfiguration("Oxy");

  //-------------------------------------
  // Define the decay channels
  G4MoleculeDefinition* water = G4H2O::Definition();
  G4MolecularDissociationChannel* decCh1;
  G4MolecularDissociationChannel* decCh2;
  G4MolecularDissociationChannel* decCh3;
  G4MolecularDissociationChannel* decCh4;
  G4MolecularDissociationChannel* decCh5;

  G4ElectronOccupancy* occ = new G4ElectronOccupancy(*(water->GetGroundStateElectronOccupancy()));

  //////////////////////////////////////////////////////////
  //            EXCITATIONS                               //
  //////////////////////////////////////////////////////////
  G4DNAWaterExcitationStructure waterExcitation;
  //--------------------------------------------------------
  //---------------Excitation on the fifth layer------------

  decCh1 = new G4MolecularDissociationChannel("A^1B_1_Relax");
  decCh2 = new G4MolecularDissociationChannel("A^1B_1_DissociDecay");
  // Decay 1 : OH + H
  decCh1->SetEnergy(waterExcitation.ExcitationEnergy(0));
  decCh1->SetProbability(0.35);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::NoDisplacement);

  decCh2->AddProduct(OH);
  decCh2->AddProduct(H);
  decCh2->SetProbability(0.65);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

  occ->RemoveElectron(4, 1);  // this is the transition form ground state to
  occ->AddElectron(5, 1);  // the first unoccupied orbital: A^1B_1

  water->NewConfigurationWithElectronOccupancy("A^1B_1", *occ);
  water->AddDecayChannel("A^1B_1", decCh1);
  water->AddDecayChannel("A^1B_1", decCh2);

  //--------------------------------------------------------
  //---------------Excitation on the fourth layer-----------
  decCh1 = new G4MolecularDissociationChannel("B^1A_1_Relax_Channel");
  decCh2 = new G4MolecularDissociationChannel("B^1A_1_DissociDecay");
  decCh3 = new G4MolecularDissociationChannel("B^1A_1_AutoIoni_Channel");
  decCh4 = new G4MolecularDissociationChannel("A^1B_1_DissociDecay");
  decCh5 = new G4MolecularDissociationChannel("B^1A_1_DissociDecay2");

  // Decay 1 : energy
  decCh1->SetEnergy(waterExcitation.ExcitationEnergy(1));
  decCh1->SetProbability(0.175);

  // Decay 2 : 2OH + H_2
  decCh2->AddProduct(H2);
  decCh2->AddProduct(OH);
  decCh2->AddProduct(OH);
  decCh2->SetProbability(0.0325);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay);

  // Decay 3 : OH + H_3Op + e_aq
  decCh3->AddProduct(OH);
  decCh3->AddProduct(H3O);
  decCh3->AddProduct(e_aq);
  decCh3->SetProbability(0.50);
  decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay 4 :  H + OH
  decCh4->AddProduct(H);
  decCh4->AddProduct(OH);
  decCh4->SetProbability(0.2535);
  decCh4->SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

  // Decay 5 : 2H + O
  decCh5->AddProduct(O);
  decCh5->AddProduct(H);
  decCh5->AddProduct(H);
  decCh5->SetProbability(0.039);
  decCh5->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay2);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3);  // this is the transition form ground state to
  occ->AddElectron(5, 1);  // the first unoccupied orbital: B^1A_1

  water->NewConfigurationWithElectronOccupancy("B^1A_1", *occ);
  water->AddDecayChannel("B^1A_1", decCh1);
  water->AddDecayChannel("B^1A_1", decCh2);
  water->AddDecayChannel("B^1A_1", decCh3);
  water->AddDecayChannel("B^1A_1", decCh4);
  water->AddDecayChannel("B^1A_1", decCh5);

  //-------------------------------------------------------
  //-------------------Excitation of 3rd layer-----------------
  decCh1 = new G4MolecularDissociationChannel("Exci3rdLayer_AutoIoni_Channel");
  decCh2 = new G4MolecularDissociationChannel("Exci3rdLayer_Relax_Channel");

  // Decay channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);

  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(2));
  decCh2->SetProbability(0.5);

  // Electronic configuration of this decay
  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->AddElectron(5, 1);

  // Configure the water molecule
  water->NewConfigurationWithElectronOccupancy("Exci3rdLayer", *occ);
  water->AddDecayChannel("Exci3rdLayer", decCh1);
  water->AddDecayChannel("Exci3rdLayer", decCh2);

  //-------------------------------------------------------
  //-------------------Excitation of 2nd layer-----------------
  decCh1 = new G4MolecularDissociationChannel("Exci2ndLayer_AutoIoni_Channel");
  decCh2 = new G4MolecularDissociationChannel("Exci2ndLayer_Relax_Channel");

  // Decay Channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);

  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(3));
  decCh2->SetProbability(0.5);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->AddElectron(5, 1);

  water->NewConfigurationWithElectronOccupancy("Exci2ndLayer", *occ);
  water->AddDecayChannel("Exci2ndLayer", decCh1);
  water->AddDecayChannel("Exci2ndLayer", decCh2);

  //-------------------------------------------------------
  //-------------------Excitation of 1st layer-----------------
  decCh1 = new G4MolecularDissociationChannel("Exci1stLayer_AutoIoni_Channel");
  decCh2 = new G4MolecularDissociationChannel("Exci1stLayer_Relax_Channel");

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->AddElectron(5, 1);

  // Decay Channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);
  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(4));
  decCh2->SetProbability(0.5);

  water->NewConfigurationWithElectronOccupancy("Exci1stLayer", *occ);
  water->AddDecayChannel("Exci1stLayer", decCh1);
  water->AddDecayChannel("Exci1stLayer", decCh2);

  /////////////////////////////////////////////////////////
  //                  IONISATION                         //
  /////////////////////////////////////////////////////////
  //--------------------------------------------------------
  //------------------- Ionisation -------------------------

  decCh1 = new G4MolecularDissociationChannel("Ioni_Channel");

  // Decay Channel 1 : : OH + H_3Op
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(1);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  // this is a ionized h2O with a hole in its last orbital
  water->NewConfigurationWithElectronOccupancy("Ioni5", *occ);
  water->AddDecayChannel("Ioni5", decCh1);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni4", *occ);
  water->AddDecayChannel("Ioni4", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni3", *occ);
  water->AddDecayChannel("Ioni3", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni2", *occ);
  water->AddDecayChannel("Ioni2", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni1", *occ);
  water->AddDecayChannel("Ioni1", new G4MolecularDissociationChannel(*decCh1));

  //////////////////////////////////////////////////////////
  //            Double ionisation                         //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("DoubleIonisation_Channel1");
  decCh2 = new G4MolecularDissociationChannel("DoubleIonisation_Channel2");
  decCh3 = new G4MolecularDissociationChannel("DoubleIonisation_Channel3");

  // Decay Channel #1: H2O^2+ -> 2H+ + O(3P) -> 2H3O+ + O(3P)
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(O);
  decCh1->SetProbability(0.29);
  decCh1->SetDisplacementType(
    G4DNAWaterDissociationDisplacer::DoubleIonisation_DissociationDecay1);

  // Decay Channel #2: H2O^2+ -> H+ + H* + O+ -> 2H3O+ + H* + *OH + O(3P)
  decCh2->AddProduct(H3O);
  decCh2->AddProduct(H3O);
  decCh2->AddProduct(H);
  decCh2->AddProduct(OH);
  decCh2->AddProduct(O);
  decCh2->SetProbability(0.16);
  decCh2->SetDisplacementType(
    G4DNAWaterDissociationDisplacer::DoubleIonisation_DissociationDecay2);

  // Decay Channel #3: H2O^2+ -> H+ + OH+ -> 2H3O+ + O(3P)
  decCh3->AddProduct(H3O);
  decCh3->AddProduct(H3O);
  decCh3->AddProduct(O);
  decCh3->SetProbability(0.55);
  decCh3->SetDisplacementType(
    G4DNAWaterDissociationDisplacer::DoubleIonisation_DissociationDecay3);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation15", *occ);
  water->AddDecayChannel("DoubleIonisation15", decCh1);
  water->AddDecayChannel("DoubleIonisation15", decCh2);
  water->AddDecayChannel("DoubleIonisation15", decCh3);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation14", *occ);
  water->AddDecayChannel("DoubleIonisation14",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation14",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation14",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation13", *occ);
  water->AddDecayChannel("DoubleIonisation13",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation13",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation13",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation12", *occ);
  water->AddDecayChannel("DoubleIonisation12",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation12",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation12",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation11", *occ);
  water->AddDecayChannel("DoubleIonisation11",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation11",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation11",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation10", *occ);
  water->AddDecayChannel("DoubleIonisation10",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation10",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation10",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation9", *occ);
  water->AddDecayChannel("DoubleIonisation9",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation9",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation9",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation8", *occ);
  water->AddDecayChannel("DoubleIonisation8",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation8",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation8",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation7", *occ);
  water->AddDecayChannel("DoubleIonisation7",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation7",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation7",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation6", *occ);
  water->AddDecayChannel("DoubleIonisation6",
                        new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation6",
                        new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation6",
                        new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation5", *occ);
  water->AddDecayChannel("DoubleIonisation5",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation5",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation5",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation4", *occ);
  water->AddDecayChannel("DoubleIonisation4",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation4",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation4",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation3", *occ);
  water->AddDecayChannel("DoubleIonisation3",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation3",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation3",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation2", *occ);
  water->AddDecayChannel("DoubleIonisation2",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation2",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation2",
                         new G4MolecularDissociationChannel(*decCh3));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("DoubleIonisation1", *occ);
  water->AddDecayChannel("DoubleIonisation1",
                         new G4MolecularDissociationChannel(*decCh1));
  water->AddDecayChannel("DoubleIonisation1",
                         new G4MolecularDissociationChannel(*decCh2));
  water->AddDecayChannel("DoubleIonisation1",
                         new G4MolecularDissociationChannel(*decCh3));

  //////////////////////////////////////////////////////////
  //            Triple ionisation                         //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("TripleIonisation_Channel");

  //Decay Channel 1 : H2O^3+ -> 3H_3Op + OH + O(3P)
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(O);
  decCh1->SetProbability(1.0);
  decCh1->SetDisplacementType(
    G4DNAWaterDissociationDisplacer::TripleIonisation_DissociationDecay);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation30", *occ);
  water->AddDecayChannel("TripleIonisation30", decCh1);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation29", *occ);
  water->AddDecayChannel("TripleIonisation29",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation28", *occ);
  water->AddDecayChannel("TripleIonisation28",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation27", *occ);
  water->AddDecayChannel("TripleIonisation27",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation26", *occ);
  water->AddDecayChannel("TripleIonisation26",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation25", *occ);
  water->AddDecayChannel("TripleIonisation25",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation24", *occ);
  water->AddDecayChannel("TripleIonisation24",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation23", *occ);
  water->AddDecayChannel("TripleIonisation23",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation22", *occ);
  water->AddDecayChannel("TripleIonisation22",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation21", *occ);
  water->AddDecayChannel("TripleIonisation21",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 2);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation20", *occ);
  water->AddDecayChannel("TripleIonisation20",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 2);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation19", *occ);
  water->AddDecayChannel("TripleIonisation19",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 2);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation18", *occ);
  water->AddDecayChannel("TripleIonisation18",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 2);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation17", *occ);
  water->AddDecayChannel("TripleIonisation17",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation16", *occ);
  water->AddDecayChannel("TripleIonisation16",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation15", *occ);
  water->AddDecayChannel("TripleIonisation15",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation14", *occ);
  water->AddDecayChannel("TripleIonisation14",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation13", *occ);
  water->AddDecayChannel("TripleIonisation13",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation12", *occ);
  water->AddDecayChannel("TripleIonisation12",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation11", *occ);
  water->AddDecayChannel("TripleIonisation11",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation10", *occ);
  water->AddDecayChannel("TripleIonisation10",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation9", *occ);
  water->AddDecayChannel("TripleIonisation9",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation8", *occ);
  water->AddDecayChannel("TripleIonisation8",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation7", *occ);
  water->AddDecayChannel("TripleIonisation7",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation6", *occ);
  water->AddDecayChannel("TripleIonisation6",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation5", *occ);
  water->AddDecayChannel("TripleIonisation5",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation4", *occ);
  water->AddDecayChannel("TripleIonisation4",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation3", *occ);
  water->AddDecayChannel("TripleIonisation3",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation2", *occ);
  water->AddDecayChannel("TripleIonisation2",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("TripleIonisation1", *occ);
  water->AddDecayChannel("TripleIonisation1",
                         new G4MolecularDissociationChannel(*decCh1));


  //////////////////////////////////////////////////////////
  //            Quadruple ionisation                      //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("QuadrupleIonisation_Channel");

  //Decay Channel 1 : H2O^4+ -> 4H_3Op + 2OH + O(3P)
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(O);
  decCh1->SetProbability(1.0);
  decCh1->SetDisplacementType(
    G4DNAWaterDissociationDisplacer::QuadrupleIonisation_DissociationDecay);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation1", *occ);
  water->AddDecayChannel("QuadrupleIonisation1", decCh1);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation2", *occ);
  water->AddDecayChannel("QuadrupleIonisation2",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation3", *occ);
  water->AddDecayChannel("QuadrupleIonisation3",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation4", *occ);
  water->AddDecayChannel("QuadrupleIonisation4",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation5", *occ);
  water->AddDecayChannel("QuadrupleIonisation5",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation6", *occ);
  water->AddDecayChannel("QuadrupleIonisation6",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation7", *occ);
  water->AddDecayChannel("QuadrupleIonisation7",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation8", *occ);
  water->AddDecayChannel("QuadrupleIonisation8",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation9", *occ);
  water->AddDecayChannel("QuadrupleIonisation9",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation10", *occ);
  water->AddDecayChannel("QuadrupleIonisation10",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation11", *occ);
  water->AddDecayChannel("QuadrupleIonisation11",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation12", *occ);
  water->AddDecayChannel("QuadrupleIonisation12",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation13", *occ);
  water->AddDecayChannel("QuadrupleIonisation13",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation14", *occ);
  water->AddDecayChannel("QuadrupleIonisation14",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation15", *occ);
  water->AddDecayChannel("QuadrupleIonisation15",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation16", *occ);
  water->AddDecayChannel("QuadrupleIonisation16",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation17", *occ);
  water->AddDecayChannel("QuadrupleIonisation17",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation18", *occ);
  water->AddDecayChannel("QuadrupleIonisation18",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation19", *occ);
  water->AddDecayChannel("QuadrupleIonisation19",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation20", *occ);
  water->AddDecayChannel("QuadrupleIonisation20",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation21", *occ);
  water->AddDecayChannel("QuadrupleIonisation21",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation22", *occ);
  water->AddDecayChannel("QuadrupleIonisation22",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation23", *occ);
  water->AddDecayChannel("QuadrupleIonisation23",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation24", *occ);
  water->AddDecayChannel("QuadrupleIonisation24",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation25", *occ);
  water->AddDecayChannel("QuadrupleIonisation25",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation26", *occ);
  water->AddDecayChannel("QuadrupleIonisation26",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation27", *occ);
  water->AddDecayChannel("QuadrupleIonisation27",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation28", *occ);
  water->AddDecayChannel("QuadrupleIonisation28",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation29", *occ);
  water->AddDecayChannel("QuadrupleIonisation29",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation30", *occ);
  water->AddDecayChannel("QuadrupleIonisation30",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation31", *occ);
  water->AddDecayChannel("QuadrupleIonisation31",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation32", *occ);
  water->AddDecayChannel("QuadrupleIonisation32",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(0, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation33", *occ);
  water->AddDecayChannel("QuadrupleIonisation33",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation34", *occ);
  water->AddDecayChannel("QuadrupleIonisation34",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  occ->RemoveElectron(4, 1);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation35", *occ);
  water->AddDecayChannel("QuadrupleIonisation35",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(1, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation36", *occ);
  water->AddDecayChannel("QuadrupleIonisation36",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation37", *occ);
  water->AddDecayChannel("QuadrupleIonisation37",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation38", *occ);
  water->AddDecayChannel("QuadrupleIonisation38",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 2);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation39", *occ);
  water->AddDecayChannel("QuadrupleIonisation39",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(2, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation40", *occ);
  water->AddDecayChannel("QuadrupleIonisation40",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation41", *occ);
  water->AddDecayChannel("QuadrupleIonisation41",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 2);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation42", *occ);
  water->AddDecayChannel("QuadrupleIonisation42",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(3, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation43", *occ);
  water->AddDecayChannel("QuadrupleIonisation43",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 2);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation44", *occ);
  water->AddDecayChannel("QuadrupleIonisation44",
                         new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 2);
  occ->RemoveElectron(4, 2);
  water->NewConfigurationWithElectronOccupancy("QuadrupleIonisation45", *occ);
  water->AddDecayChannel("QuadrupleIonisation45",
                         new G4MolecularDissociationChannel(*decCh1));

  //////////////////////////////////////////////////////////
  //            Dissociative Attachment                   //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("DissociAttachment_ch1");

  // Decay 1 : OHm + H
  decCh1->AddProduct(H2);
  decCh1->AddProduct(OHm);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(1);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::DissociativeAttachment);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->AddElectron(5, 1);  // H_2O^-

  water->NewConfigurationWithElectronOccupancy("DissociAttachment_ch1", *occ);
  water->AddDecayChannel("DissociAttachment_ch1", decCh1);

  //////////////////////////////////////////////////////////
  //            Electron-hole recombination               //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay1");
  decCh2 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay2");
  decCh3 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay3");
  decCh4 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay4");

  // Decay 1 : 2OH + H_2
  decCh1->AddProduct(H2);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(0.1365);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay);

  // Decay 2 : OH + H
  decCh2->AddProduct(OH);
  decCh2->AddProduct(H);
  decCh2->SetProbability(0.3575);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

  // Decay 3 : 2H + O(3p)
  decCh3->AddProduct(O);
  decCh3->AddProduct(H);
  decCh3->AddProduct(H);
  decCh3->SetProbability(0.156);
  decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay2);

  // Decay 4 : relaxation
  decCh4->SetProbability(0.35);

  const auto pH2Ovib = G4H2O::Definition()->NewConfiguration("H2Ovib");
  assert(pH2Ovib != nullptr);

  water->AddDecayChannel(pH2Ovib, decCh1);
  water->AddDecayChannel(pH2Ovib, decCh2);
  water->AddDecayChannel(pH2Ovib, decCh3);
  water->AddDecayChannel(pH2Ovib, decCh4);

  delete occ;
}
