#ifndef MOLFORMAT_CONVERTER_H
#define MOLFORMAT_CONVERTER_H 1
#include "UMDFP.h"
#include "UMFHelpers.h"
#include <iomanip>

class GenericMoleculeFileFormat
{
public:
    GenericMoleculeFileFormat() {}
    virtual std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out, const std::string& charge_method="none") const = 0;
};

class Mol2Format : public GenericMoleculeFileFormat
{
    std::string resname="UNL"; // Default residue name for all molecules (can be modified in the future to allow different residue names for different molecules if needed) 
public:
    Mol2Format() {}
    std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out, const std::string& charge_method="NO_CHARGES") const override
    {
        char atom_line[100];
        out << "# Name: " << molecule.getName() << "\n";
        out << "# SMILES: " << molecule.getSMILES().getData() << "\n";
        out << "# Net charge: " << molecule.computeNetCharge() << "\n";
        out << "@<TRIPOS>MOLECULE\n";
        out << molecule.getName() << "\n";
        out << molecule.getNumAtoms() << " " << molecule.getNumBonds() << " 0 0 0\n"; // Num atoms, num bonds, num substructures, num features, num sets (the last three are set to 0 for simplicity)
        out << "SMALL\n"; // Molecule type 
        out << charge_method << "\n\n"; // Charge type

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
    std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out, const std::string& charge_method="none") const override
    {
        // Comment with SMILES String
        
        //out << "\n" << molecule.getName() << "\n"; // Line 1: Molecule name
        out << molecule.getName() << "\n";
        out << "  UMDBabel0.12026   3D\n"; // Line 2: Program information (can be modified in the future to include more specific information if needed)
        out << "\n";
        out << molecule.getNumAtoms() << " " << molecule.getNumBonds() << " 0 0 0  0  0  0  0  0999 V3000\n"; // Line 3: Number of atoms, number of bonds, and other counts (the last three are set to 0 for simplicity)

        // Atom block
        for(int i=0;i<molecule.getNumAtoms();i++)
        {
            const UMDAtom& atom = molecule.getAtom(i);
            out << std::fixed << std::setprecision(4) << std::setw(10) << atom.getX() << std::setw(10) << atom.getY() << std::setw(10) << atom.getZ() << " " << std::left << std::setw(2) << atom.getElement() << "  0  0  0  0  0  0  0  0  0  0  0  0\n"; // x, y, z, element symbol
        }

        // Bond block
        for(int i=0;i<molecule.getNumBonds();i++)
        {
            const UMDBond& bond = molecule.getBond(i);
            int btype= bond.getBondType();
            if(molecule.getAtom(bond.getAtom1ID()).getData().aromatic && molecule.getAtom(bond.getAtom2ID()).getData().aromatic) btype=3; // Set bond type to aromatic if both atoms are aromatic

            out << std::setw(3) << bond.getAtom1ID()+1 << std::setw(3) << bond.getAtom2ID()+1 << std::setw(3) << ((btype%2==1)?"4":std::to_string(btype/2)) << "\n"; // Atom1 ID, atom2 ID, bond type
        }
        int num_formal_charges=0;
        for(int i=0;i<molecule.getNumAtoms();i++)        {
            int fchg=molecule.getAtom(i).getFormalCharge();
            if(fchg!=0) num_formal_charges++;
        }
        if(num_formal_charges>0)
        {
            out << "M  CHG  " << num_formal_charges << std::right; // Line indicating the number of atoms with formal charges
            for(int i=0;i<molecule.getNumAtoms();i++)
            {
                int formal_charge = molecule.getAtom(i).getFormalCharge();
                if(formal_charge!=0) out << std::setw(4) << (i+1) << std::setw(4) << formal_charge; // Atom ID and formal charge
            }
            out << "\n" << std::left;
        }
        
        out << "M  END\n"; // End of molecule block for SDF format
        out << ">" << "<SMILES>\n" << molecule.getSMILES().getData() << "\n\n"; // SMILES string as a comment line in the SDF file

        out << "$$$$\n"; // End of molecule block
        return out;
    }
};

class PDBQTFormat : public GenericMoleculeFileFormat
{
    bool write_H; // Whether to include hydrogens in the output (can be set to false for docking programs that add hydrogens automatically, or true for programs that require explicit hydrogens)
public:
    PDBQTFormat(bool write_H=true) {this->write_H=write_H;}
    std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& _out, const std::string& charge_method="none") const override
    {
        std::stringstream out;
        std::stringstream remarks;
        remarks << "REMARK Name: " << molecule.getName() << "\n";
        remarks << "REMARK SMILES: " << molecule.getSMILES().getData() << "\n";
        remarks << "REMARK Net charge: " << molecule.computeNetCharge() << "\n";
        remarks << "REMARK Charge method: " << charge_method << "\n";
        // Generate branches
        std::vector<std::pair<int,int>> tree_rep=generateDFSTraversalOrder(molecule, 0, write_H); // Generate a DFS traversal order of the molecule starting from the first atom (index 0), including hydrogens if write_H is true
        // Find all atoms in a ring
        std::pair<std::vector<bool>,std::vector<std::vector<int>>> in_ring_data = computeAtomRingStatus(molecule);
        const std::vector<bool>& in_ring = in_ring_data.first;
        const std::vector<std::vector<int>>& rings = in_ring_data.second;
        
        std::vector<bool> parent_list(molecule.getNumAtoms(), false);
        std::vector<std::string> last_branch={"ROOT"};
        std::string branch_name,atom_type;
        int torsdof=0;
        

        std::vector<int> written_index(molecule.getNumAtoms(), -1);
        int written_count=0;
        
        for(const auto& [atom_index, parent_index] : tree_rep)
        {
            const UMDAtom& atom = molecule.getAtom(atom_index);
            if(AtomIsHydrogen(atom) && !write_H) continue; // Skip hydrogens if write_H is false
            if(parent_index==-1) out << "ROOT\n";
            else if(getNeighborIndices(molecule, atom_index, false).size()==1) {}
            else
            {
                bool is_tors=false;
                if(parent_index!=-1 && getNeighborIndices(molecule, parent_index, false).size()>1 && getNeighborIndices(molecule, atom_index, false).size()>1)
                {
                    bool ring_matched=false;
                    for(const std::vector<int>& ring : rings)
                    {
                        if(contains(ring, parent_index) && contains(ring, atom_index))
                        {
                            ring_matched=true;
                            break;
                        }
                    }
                    // Exceptions
                    // 1. Amide
                    // 2. Non-sigma bonds
                    // 3. In same ring
                    if(ring_matched) {}
                    else if(BondIsAmide(molecule, parent_index, atom_index)) {}
                    else if(getBondBetweenAtoms(molecule, parent_index, atom_index).getBondType()>3) {}
                    else {torsdof++; is_tors=true;}
                }
                
                if(parent_list[parent_index] || (last_branch.back()=="ROOT" && is_tors)) // If the parent already has a child, this is a new branch
                {
                    out << "END" << last_branch.back() << "\n"; // End the last branch
                    last_branch.pop_back();
                    /*char raw_branch_name[32];
                    sprintf(raw_branch_name, "BRANCH   %d  %d", parent_index+1, atom_index+1);
                    raw_branch_name[31]='\0'; // Ensure null termination*/
                    branch_name = "BRANCH   " + std::to_string(written_index[parent_index]+1) + "  " + std::to_string(written_index[atom_index]+1);
                    out << branch_name << "\n";
                    last_branch.push_back(branch_name);
                } 
                else parent_list[parent_index]=true; // Mark the parent as having a child
            }

            char atom_line[256];
            atom_type=(AtomIsCarbon(atom) && atom.isAromatic()) ? "" : atom.getElement();
            if(atom.isAromatic()) atom_type += "A"; // Append A to the atom type if the atom is aromatic (this is a common convention in PDBQT files to indicate aromatic atoms, but can be modified in the future if needed to follow a different convention or include more specific information)
            sprintf(atom_line, "ATOM  %5d  %-3s UNL     0    %8.3f%8.3f%8.3f%6.2f%6.2f    %+3.3f %-3s\n", written_count+1, atom.getElement().c_str(), atom.getX(), atom.getY(), atom.getZ(), 0.0f, 0.0f, atom.getCharge(), atom_type.c_str());
            atom_line[255]='\0'; // Ensure null termination
            out << atom_line; // Atom ID, atom type, x, y, z, charge
            written_index[atom_index]=written_count+1;
            written_count++;
        }
        while(!last_branch.empty())
        {
            out << "END" << last_branch.back() << "\n"; // End any remaining open branches
            last_branch.pop_back();
        }
        out << "TORSDOF " << torsdof << "\n"; // Number of rotatable bonds (this is a simplification, as it assumes that every branch point corresponds to a rotatable bond, which may not always be the case, but can be modified in the future to include a more accurate calculation of rotatable bonds if needed)
        for(const std::vector<int>& ring : rings)
        {
            remarks << "REMARK RING ";
            for(int atom_index : ring)
            {
                remarks << written_index[atom_index] << " (" << molecule.getAtom(atom_index).getElement() << ") ";
            }
            remarks << "\n";
        }
        _out << remarks.str() << out.str(); 
        return _out;
    }
};

class SMILESFormat : public GenericMoleculeFileFormat
{
public:
    SMILESFormat() {}
    std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out, const std::string& charge_method="none") const override {out << molecule.getSMILES().getData() << " " << molecule.getName() << "\n"; return out;}
};

#endif
