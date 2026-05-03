#include <iostream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include "UMDFP.h"

class RDKitMolecule
{
protected:
    RDKit::ROMol* mol1=nullptr;
public:
    RDKitMolecule(const std::string& smiles)
    {
        mol1 = RDKit::SmilesToMol(smiles);
        if (!mol1) {
            std::cerr << "Error: Invalid SMILES string: " << smiles << std::endl;
        }
    }
    RDKitMolecule(const UMDMolecule& mol) : RDKitMolecule(reinterpret_cast<const char*>(mol.getSMILES().getData())) {}

    ~RDKitMolecule() {delete mol1;}

    std::string toSmiles() const
    {
        if (mol1) {
            return RDKit::MolToSmiles(*mol1);
        } else {
            std::cerr << "Error: Molecule is not initialized." << std::endl;
            return "";
        }
    }
};
