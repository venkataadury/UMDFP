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
    virtual std::vector<std::vector<float>> readCoordinates(std::istream& in) const {return std::vector<std::vector<float>>();} // Default implementation returns an empty vector, can be overridden by formats that support reading coordinates
    inline virtual bool hasCoordinates() const {return false;} // Default implementation returns false, can be overridden by formats that support reading coordinates
};

class Mol2Format : public GenericMoleculeFileFormat
{
    std::string resname="UNL"; // Default residue name for all molecules (can be modified in the future to allow different residue names for different molecules if needed) 
    bool write_H=true; // Whether to include hydrogens in the output (can be set to false for docking programs that add hydrogens automatically, or true for programs that require explicit hydrogens)
public:
    Mol2Format(bool write_h=true) : write_H(write_h) {}
    std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out, const std::string& charge_method="NO_CHARGES") const override
    {
        int atom_count=molecule.getNumAtoms(), bond_count=molecule.getNumBonds();
        char atom_line[100];
        std::stringstream header, rest;
        header << "# Name: " << molecule.getName() << "\n";
        header << "# SMILES: " << molecule.getSMILES().getData() << "\n";
        header << "# Net charge: " << molecule.computeNetCharge() << "\n";
        header << "@<TRIPOS>MOLECULE\n";
        header << molecule.getName() << "\n";
        //header << atom_count << " " << bond_count << " 0 0 0\n"; // Num atoms, num bonds, num substructures, num features, num sets (the last three are set to 0 for simplicity)
        rest << "SMALL\n"; // Molecule type 
        rest << charge_method << "\n\n"; // Charge type
        std::vector<int> written_index(molecule.getNumAtoms(), -1);
        int written_count=0;

        rest << "@<TRIPOS>ATOM\n";
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
            if(AtomIsHydrogen(atom) && !write_H) continue; // Skip hydrogens if write_H is false
            //f.write(f"{a[0]:>6} {a[1]:<4} {a[2]:>10.4f} {a[3]:>10.4f} {a[4]:>10.4f} {a[5]:<5} {a[6]:<4} {a[7]:<5} {a[8]:>10.4f}\n")
            sprintf(atom_line, " %6d %-4s    %10.4f%10.4f%10.4f %-5s   1  %-5s   %10.4f\n", written_count+1, atom.getElement().c_str(), atom.getX(), atom.getY(), atom.getZ(), atom_type.c_str(), resname.c_str(), atom.getCharge());
            rest << atom_line; // Atom ID, atom type, x, y, z, charge
            written_index[i]=written_count+1;
            written_count++;
        }

        rest << "@<TRIPOS>BOND\n";
        atom_count=written_count;
        written_count=0;
        for(int i=0;i<molecule.getNumBonds();i++)
        {
            const UMDBond& bond = molecule.getBond(i);
            if(AtomIsHydrogen(molecule.getAtom(bond.getAtom1ID())) || AtomIsHydrogen(molecule.getAtom(bond.getAtom2ID())) && !write_H) continue; // Skip bonds involving hydrogens if write_H is false
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
            rest << written_count+1 << " " << written_index[bond.getAtom1ID()] << " " << written_index[bond.getAtom2ID()] << " " << btype_str << "\n"; // Bond ID, atom1 ID, atom2 ID, bond type
            written_count++;
        }
        bond_count=written_count;
        header << atom_count << " " << bond_count << " 0 0 0\n"; // Update the header with the correct number of atoms and bonds after filtering
        out << header.str() << rest.str();
        return out;
    }

    inline void setResidueName(const std::string& resname) {this->resname = resname;}
    inline std::string getResidueName() const {return resname;}

    std::vector<std::vector<float>> readCoordinates(std::istream& in) const override
    {
        std::vector<std::vector<float>> coordinates;
        std::string line;
        // Read until we find the @<TRIPOS>ATOM section
        while(std::getline(in, line))
        {
            if(line.find("@<TRIPOS>ATOM")!=std::string::npos) break;
        }
        // Read atom lines until we reach the next section
        while(std::getline(in, line))
        {
            if(line.empty()) continue;
            if(line[0]=='@') break; // Next section
            std::istringstream iss(line);
            int atom_id; std::string atom_name; float x, y, z;
            if(!(iss >> atom_id >> atom_name >> x >> y >> z)) continue; // Skip lines that don't match the expected format
            coordinates.push_back({x, y, z});
        }
        return coordinates;
    }
    inline bool hasCoordinates() const override {return true;}
};

class SDFFormat : public GenericMoleculeFileFormat
{
    bool write_H=true; // Whether to include hydrogens in the output (can be set to false for docking programs that add hydrogens automatically, or true for programs that require explicit hydrogens)
public:
    SDFFormat(bool write_h=true) : write_H(write_h) {}
    std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out, const std::string& charge_method="none") const override
    {
        // Comment with SMILES String
        std::stringstream header, rest;
        
        //out << "\n" << molecule.getName() << "\n"; // Line 1: Molecule name
        header << molecule.getName() << "\n";
        header << "  UMDBabel0.12026   3D\n"; // Line 2: Program information (can be modified in the future to include more specific information if needed)
        header << "\n";
        // header << molecule.getNumAtoms() << " " << molecule.getNumBonds() << " 0 0 0  0  0  0  0  0999 V3000\n"; // Line 3: Number of atoms, number of bonds, and other counts (the last three are set to 0 for simplicity)

        // Atom block
        int atom_count=0, bond_count=0;
        for(int i=0;i<molecule.getNumAtoms();i++)
        {
            const UMDAtom& atom = molecule.getAtom(i);
            if(AtomIsHydrogen(atom) && !write_H) continue; // Skip hydrogens if write_H is false
            rest << std::fixed << std::setprecision(4) << std::right << std::setw(10) << atom.getX() << std::setw(10) << atom.getY() << std::setw(10) << atom.getZ() << " " << std::left << std::setw(2) << atom.getElement() << "  0  0  0  0  0  0  0  0  0  0  0  0\n"; // x, y, z, element symbol
            atom_count++;
        }

        // Bond block
        for(int i=0;i<molecule.getNumBonds();i++)
        {
            const UMDBond& bond = molecule.getBond(i);
            if(AtomIsHydrogen(molecule.getAtom(bond.getAtom1ID())) || AtomIsHydrogen(molecule.getAtom(bond.getAtom2ID())) && !write_H) continue; // Skip bonds involving hydrogens if write_H is false
            int btype= bond.getBondType();
            if(molecule.getAtom(bond.getAtom1ID()).getData().aromatic && molecule.getAtom(bond.getAtom2ID()).getData().aromatic) btype=3; // Set bond type to aromatic if both atoms are aromatic

            rest << std::setw(3) << bond.getAtom1ID()+1 << std::setw(3) << bond.getAtom2ID()+1 << std::setw(3) << ((btype%2==1)?"4":std::to_string(btype/2)) << "\n"; // Atom1 ID, atom2 ID, bond type
            bond_count++;
        }
        header << atom_count << " " << bond_count << " 0 0 0  0  0  0  0  0999 V3000\n"; // Line 3: Number of atoms, number of bonds, and other counts (the last three are set to 0 for simplicity)
        out << header.str() << rest.str();

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
    inline bool hasCoordinates() const override {return true;}
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
        std::vector<std::vector<int>> branch_members; branch_members.push_back({}); // Start with the root branch containing the first atom (index 0)
        std::string branch_name,atom_type;
        std::vector<int> written_sequence;
        int torsdof=0;
        

        std::vector<int> written_index(molecule.getNumAtoms(), -1);
        int written_count=0;
        
        for(const auto& [atom_index, parent_index] : tree_rep)
        {
            const UMDAtom& atom = molecule.getAtom(atom_index);
            if(AtomIsHydrogen(atom) && !write_H) continue; // Skip hydrogens if write_H is false
            if(parent_index==-1) {out << "ROOT\n"; branch_members.push_back({atom_index});} // Start the root branch with (usually) the first atom (index 0)
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
                    while((!contains(branch_members.back(),parent_index) || last_branch.back()=="ROOT") && last_branch.size()) // Pop branches until we find the branch that contains the parent atom of the current atom
                    {
                        out << "END" << last_branch.back() << "\n"; // End the last branch
                        last_branch.pop_back();
                        branch_members.pop_back();
                    }
                    /*char raw_branch_name[32];
                    sprintf(raw_branch_name, "BRANCH   %d  %d", parent_index+1, atom_index+1);
                    raw_branch_name[31]='\0'; // Ensure null termination*/
                    branch_name = "BRANCH   " + std::to_string(written_index[parent_index]) + "  " + std::to_string(written_count+1);
                    out << branch_name << "\n";
                    last_branch.push_back(branch_name);
                    branch_members.push_back({});
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
            written_sequence.push_back(atom_index);
            branch_members.back().push_back(atom_index);
            written_count++;
        }
        while(!last_branch.empty())
        {
            out << "END" << last_branch.back() << "\n"; // End any remaining open branches
            last_branch.pop_back();
        }
        out << "TORSDOF " << torsdof << "\n"; // Number of rotatable bonds (this is a simplification, as it assumes that every branch point corresponds to a rotatable bond, which may not always be the case, but can be modified in the future to include a more accurate calculation of rotatable bonds if needed)
        remarks << "REMARK AtomSequence ";
        for(int aind : written_sequence) remarks << aind << " "; remarks << "\n";
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
    
    inline bool hasCoordinates() const override {return true;}
    std::vector<std::vector<float>> readCoordinates(std::istream& in) const override
    {
        // PDBQT format does not have a standard way to indicate the end of the coordinate section, so we will read until we reach a line that does not start with "ATOM" or "ROOT" or "BRANCH" or "END"
        std::vector<std::vector<float>> coordinates;
        std::string line;
        std::vector<int> atom_reorder;
        while(std::getline(in, line))
        {
            if(line.empty()) continue;
            if(line.find("REMARK AtomSequence")==0)
            {
                std::istringstream iss(line.substr(20)); // Extract the part of the line after "REMARK AtomSequence "
                int atom_index;
                while(iss >> atom_index) atom_reorder.push_back(atom_index); // Store the atom indices in the order they are written in the PDBQT file (subtract 1 to convert from 1-based to 0-based indexing)
                continue;
            }
            if(line.find("ATOM")!=0) continue; // Skip non-atom lines
            std::istringstream iss(line);
            char record_name[7], atom_name[5], res_name[5], element_symbol[5];
            int atom_id, res_id;
            float x=1.0, y=1.0, z=1.0, occupancy, temp_factor, charge;
            std::cout << line.c_str() << "\n";
            try {sscanf(line.c_str(), "%6s      %d  %4s   %3s     %d      %f  %f   %f  %f  %f    %f %3s",record_name, &atom_id, atom_name, res_name, &res_id, &x, &y, &z, &occupancy, &temp_factor, &charge, element_symbol);}
            catch(...)
            {
                std::cerr << "Warning: Failed to parse line in PDBQT file while reading coordinates: " << line << "\n";
                std::cerr << "Bad line: " << line << "\n";
                continue;
            }
            coordinates.push_back({x, y, z});
        }
        // Reorder the coordinates according to the order of atoms in the original molecule (if the REMARK AtomSequence line was present)
        if(atom_reorder.size()==coordinates.size())
        {
            std::vector<std::vector<float>> reordered_coordinates(coordinates.size());
            for(size_t i=0;i<atom_reorder.size();i++)
            {
                int original_index = atom_reorder[i];
                if(original_index<0 || original_index>=coordinates.size())
                {
                    std::cerr << "Warning: Atom index in REMARK AtomSequence line is out of bounds: " << original_index << "\n";
                    std::cerr << "Skipping reordering of coordinates\n";
                    return coordinates; // Return the coordinates in the order they were read from the file if there is an error with the REMARK AtomSequence line
                }
                reordered_coordinates[original_index] = coordinates[i];
            }
            coordinates = reordered_coordinates;
        }
        else
        {
            if(!atom_reorder.empty()) std::cerr << "Warning: REMARK AtomSequence line was found but the number of atom indices does not match the number of atoms, returning coordinates in the order they were read from the file\n";
            else std::cerr << "Warning: REMARK AtomSequence line was not found, returning coordinates in the order they were read from the file\n";
            std::cerr << "If this PDBQt file was not written from UMD format, the atom ordering might not match. Please make sure that you manually ensure the same order of atoms in the PDBQt and any template molecule loaded into UMD format.\n";
        }
        return coordinates;
    }
};

class SMILESFormat : public GenericMoleculeFileFormat
{
public:
    SMILESFormat() {}
    std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out, const std::string& charge_method="none") const override {out << molecule.getSMILES().getData() << " " << molecule.getName() << "\n"; return out;}
    std::vector<std::vector<float>> readCoordinates(std::istream& in) const override {return std::vector<std::vector<float>>();} // SMILES format does not contain coordinate information, so we return an empty vector
    inline bool hasCoordinates() const override {return false;}
};

class UMDFormat : public GenericMoleculeFileFormat
{
public:
    UMDFormat() {}
    std::ostream& formatMolecule(const UMDMolecule& molecule, std::ostream& out, const std::string& charge_method="none") const override
    {
        out << "START\n";
        out << molecule.getName() << "\n";
        out << ((char*)molecule.getSMILES().getData()) << "\n";
        out << molecule.getNumAtoms() << "\n";
        out << molecule.getNumBonds() << "\n";
        out << "; Index Element X Y Z Q Hyb Aromatic FormalCharge\n";
        for(int i=0;i<molecule.getNumAtoms();i++)
        {
            const UMDAtom& atom = molecule.getAtom(i);
            out << i << " " << atom.getElement() << " " << atom.getX() << " " << atom.getY() << " " << atom.getZ() << " " << atom.getCharge() << " " << atom.getHybridization() << " " << atom.isAromatic() << " " << atom.getFormalCharge() << "\n";
        }
        for(int i=0;i<molecule.getNumBonds();i++)
        {
            const UMDBond& bond = molecule.getBond(i);
            out << i << " " << bond.getAtom1ID() << " " << bond.getAtom2ID() << " " << bond.getBondType() << "\n";
        }
        //if(molecule.extras.getLength()>0) out.write(molecule.extras.getData(), sizeof(char)*molecule.extras.getLength());
        out << "END\n";
        return out;
    }
    std::vector<std::vector<float>> readCoordinates(std::istream& in) const override
    {
        std::vector<std::vector<float>> coordinates;
        std::string line;
        // Read until we find the START line
        while(std::getline(in, line))
        {
            if(line.find("START")!=std::string::npos) break;
        }
        // Read number of atoms and bonds
        int natoms=-1, nbonds=-1;
        bool read_name=false, read_smiles=false;
        while(std::getline(in, line))
        {
            if(line.empty()) continue;
            if(line[0]==';') continue; // Skip comment lines
            if(!read_name) {read_name=true; continue;} // Skip name line
            else if(!read_smiles) {read_smiles=true; continue;} // Skip SMILES line
            else if(natoms==-1) natoms = std::stoi(line);
            else if(nbonds==-1) {nbonds = std::stoi(line); break;}
        }
        // Read atom lines until we reach the bond block
        while(std::getline(in, line))
        {
            if(line.empty()) continue;
            if(line[0]==';') continue; // Skip comment lines
            if(line.find("END")!=std::string::npos) throw std::runtime_error("Unexpected END line while reading coordinates");
            if(coordinates.size()>=natoms) break; // If we have read all atoms, break out of the loop to avoid reading bond lines as atom lines
            std::istringstream iss(line);
            int index; std::string element; float x, y, z; float charge; unsigned short int hyb; bool aromatic; short int formal_charge;
            if(!(iss >> index >> element >> x >> y >> z >> charge >> hyb >> aromatic >> formal_charge)) continue; // Skip lines that don't match the expected format
            coordinates.push_back({x, y, z});
        }
        return coordinates;
    }
    inline bool hasCoordinates() const override {return true;}
};


static bool loadCoordinatesFromFormat(UMDMolecule& molecule, std::istream& crdfile, GenericMoleculeFileFormat& formatter)
{
    if(!formatter.hasCoordinates()) {std::cerr << "Error: Formatter does not support coordinates" << std::endl; return false;}
    std::vector<std::vector<float>> coords = formatter.readCoordinates(crdfile);
    
    if(coords.size()!=molecule.getNumAtoms())
    {
        // Indicates bad match between coordinate file and loaded molecule. Throw an error in this case since this is likely a user error
        std::cerr << "Error: Number of coordinate sets in file does not match number of atoms in molecule" << std::endl;
        return false;
    }
    for(int i=0;i<molecule.getNumAtoms();i++)
    {
        molecule.getAtom(i).getData().x = coords[i][0];
        molecule.getAtom(i).getData().y = coords[i][1];
        molecule.getAtom(i).getData().z = coords[i][2];
    }
    return true;
}
static bool loadCoordinatesFromFormat(UMDMolecule& molecule, const std::string& filename, GenericMoleculeFileFormat& formatter)
{
    if(!formatter.hasCoordinates()) return false;
    std::ifstream crdfile(filename);
    if(!crdfile)
    {
        std::cerr << "Error opening coordinate file: " << filename << std::endl;
        return false;
    }
    return loadCoordinatesFromFormat(molecule, crdfile, formatter);
}

#endif
