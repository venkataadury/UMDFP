#ifndef MOLFORMAT_CONVERTER_H
#define MOLFORMAT_CONVERTER_H 1
#include "UMDFP.h"
#include <iomanip>

class GenericMoleculeFileFormat
{
public:
    GenericMoleculeFileFormat() {}
    virtual std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out) const = 0;
};

class Mol2Format : public GenericMoleculeFileFormat
{
    std::string resname="UNL"; // Default residue name for all molecules (can be modified in the future to allow different residue names for different molecules if needed) 
public:
    Mol2Format() {}
    std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out) const override
    {
        char atom_line[100];
        out << "# Name: " << molecule.getName() << "\n";
        out << "# SMILES: " << molecule.getSMILES().getData() << "\n";
        out << "# Net charge: " << molecule.computeNetCharge() << "\n";
        out << "@<TRIPOS>MOLECULE\n";
        out << molecule.getName() << "\n";
        out << molecule.getNumAtoms() << " " << molecule.getNumBonds() << " 0 0 0\n"; // Num atoms, num bonds, num substructures, num features, num sets (the last three are set to 0 for simplicity)
        out << "SMALL\n"; // Molecule type 
        out << "NO_CHARGES\n\n"; // Charge type

        out << "@<TRIPOS>ATOM\n";
        for(int i=0;i<molecule.getNumAtoms();i++)
        {
            const UMDAtom& atom = molecule.getAtom(i);
            std::string atom_type=atom.getElement();
            if(atom.isAromatic()) atom_type += ".ar"; // Append .ar to the atom type if the atom is aromatic
            else
            {
                switch(atom.getHybridization())
                {
                    case 1: atom_type += ".1"; break; // sp hybridization
                    case 2: atom_type += ".2"; break; // sp2 hybridization
                    case 3: atom_type += ".3"; break; // sp3 hybridization
                    default: break; // Keep the element symbol as the atom type if hybridization is unknown
                }
            }
            //f.write(f"{a[0]:>6} {a[1]:<4} {a[2]:>10.4f} {a[3]:>10.4f} {a[4]:>10.4f} {a[5]:<5} {a[6]:<4} {a[7]:<5} {a[8]:>10.4f}\n")
            sprintf(atom_line, " %6d %-4s    %10.4f%10.4f%10.4f %-5s   1  %-5s   %10.4f\n", i+1, atom.getElement().c_str(), atom.getX(), atom.getY(), atom.getZ(), atom_type.c_str(), resname.c_str(), atom.getCharge());
            out << atom_line; // Atom ID, atom type, x, y, z, charge
        }

        out << "@<TRIPOS>BOND\n";
        for(int i=0;i<molecule.getNumBonds();i++)
        {
            const UMDBond& bond = molecule.getBond(i);
            int btype= bond.getBondType();
            if(molecule.getAtom(bond.getAtom1ID()).getData().aromatic && molecule.getAtom(bond.getAtom2ID()).getData().aromatic) btype=3; // Set bond type to aromatic if both atoms are aromatic

            std::string btype_str;
            switch(btype)
            {
                case 2: btype_str = "1"; break; // Single bond
                case 4: btype_str = "2"; break; // Double bond
                case 6: btype_str = "3"; break; // Triple bond
                case 3: btype_str = "ar"; break; // Aromatic bond
                default: btype_str = "1"; break; // Default to single bond if unknown
            }
            out << i+1 << " " << bond.getAtom1ID()+1 << " " << bond.getAtom2ID()+1 << " " << btype_str << "\n"; // Bond ID, atom1 ID, atom2 ID, bond type
        }
        return out;
    }

    inline void setResidueName(const std::string& resname) {this->resname = resname;}
    inline std::string getResidueName() const {return resname;}
};

class SDFFormat : public GenericMoleculeFileFormat
{
public:
    SDFFormat() {}
    std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out) const override
    {
        // Comment with SMILES String
        
        out << molecule.getName() << "\n"; // Line 1: Molecule name
        out << "  UMDBabel 0.1 2026\n"; // Line 2: Program information (can be modified in the future to include more specific information if needed)
        out << molecule.getNumAtoms() << " " << molecule.getNumBonds() << " 0 0 0\n"; // Line 3: Number of atoms, number of bonds, and other counts (the last three are set to 0 for simplicity)

        // Atom block
        for(int i=0;i<molecule.getNumAtoms();i++)
        {
            const UMDAtom& atom = molecule.getAtom(i);
            out << std::fixed << std::setprecision(4) << std::setw(10) << atom.getX() << std::setw(10) << atom.getY() << std::setw(10) << atom.getZ() << " " << atom.getElement() << "\n"; // x, y, z, element symbol
        }

        // Bond block
        for(int i=0;i<molecule.getNumBonds();i++)
        {
            const UMDBond& bond = molecule.getBond(i);
            int btype= bond.getBondType();
            if(molecule.getAtom(bond.getAtom1ID()).getData().aromatic && molecule.getAtom(bond.getAtom2ID()).getData().aromatic) btype=4; // Set bond type to aromatic if both atoms are aromatic

            out << std::setw(3) << bond.getAtom1ID()+1 << std::setw(3) << bond.getAtom2ID()+1 << std::setw(3) << btype << "\n"; // Atom1 ID, atom2 ID, bond type
        }
        
        out << "M  END\n"; // End of molecule block for SDF format
        out << ">" << "<SMILES>\n" << molecule.getSMILES().getData() << "\n\n"; // SMILES string as a comment line in the SDF file

        out << "$$$$\n"; // End of molecule block
        return out;
    }
};

#endif