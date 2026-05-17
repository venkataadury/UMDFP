#ifndef INCLUDED_UMDFP_HELPERS_H
#define INCLUDED_UMDFP_HELPERS_H
#include "UMDFP.h"
#include <functional>
#include <sstream>
#include <map>

static inline bool AtomIsHydrogen(const UMDAtom& atom) {return atom.getElement()=="H";}
static inline bool AtomIsCarbon(const UMDAtom& atom) {return atom.getElement()=="C";}
static inline bool AtomIsNitrogen(const UMDAtom& atom) {return atom.getElement()=="N";}
static bool BondIsAmide(const UMDMolecule& molecule, int atom_index1, int atom_index2)
{
    if(AtomIsCarbon(molecule.getAtom(atom_index2))) std::swap(atom_index1, atom_index2); // Ensure that atom 1 is the carbon and atom 2 is the nitrogen (if it's an amide, the carbon will be bonded to the nitrogen)
    const UMDAtom& atom1 = molecule.getAtom(atom_index1);
    const UMDAtom& atom2 = molecule.getAtom(atom_index2);
    
    if(!AtomIsCarbon(atom1)) return false; // If neither atom is carbon, this can't be an amide
    if(!AtomIsNitrogen(atom2)) return false; // If atom 2 is not nitrogen, this can't be an amide
    // Check if atom 1 is a carbonyl carbon (i.e. has a double bond to an oxygen) and atom 2 is a nitrogen bonded to the carbon
    bool carbonyl_carbon=false, nitrogen=false;
    for(int i=0;i<molecule.getNumBonds();i++)
    {
        const UMDBond& bond = molecule.getBond(i);
        if(bond.getAtom1ID()==atom_index1 || bond.getAtom2ID()==atom_index1)
        {
            int neighbor_index = (bond.getAtom1ID()==atom_index1) ? bond.getAtom2ID() : bond.getAtom1ID();
            if(molecule.getAtom(neighbor_index).getElement()=="O" && bond.getBondType()==4) carbonyl_carbon=true; // Check for double bond to oxygen
        }
    }
    bool verdict = (carbonyl_carbon && atom2.getElement()=="N");
    return verdict; // Return true if atom 1 is a carbonyl carbon and atom 2 is a nitrogen
}
template<class T> static bool contains(const std::vector<T>& vec, const T& value) {return std::find(vec.begin(), vec.end(), value) != vec.end();}

static int computeNumNeighbors(const UMDMolecule& molecule, int atom_index, bool include_hydrogens=true)
{
    int count=0;
    for(int i=0;i<molecule.getNumBonds();i++)
    {
        const UMDBond& bond = molecule.getBond(i);
        if(bond.getAtom1ID()==atom_index || bond.getAtom2ID()==atom_index)
        {
            if(!include_hydrogens && (AtomIsHydrogen(molecule.getAtom(bond.getAtom1ID())) || AtomIsHydrogen(molecule.getAtom(bond.getAtom2ID())))) continue;
            count++;
        }
    }
    return count;
}
static std::vector<int> getNeighborIndices(const UMDMolecule& molecule, int atom_index, bool include_hydrogens=true)
{
    std::vector<int> neighbors;
    for(int i=0;i<molecule.getNumBonds();i++)
    {
        const UMDBond& bond = molecule.getBond(i);
        if(bond.getAtom1ID()==atom_index || bond.getAtom2ID()==atom_index)
        {
            int neighbor_index = (bond.getAtom1ID()==atom_index) ? bond.getAtom2ID() : bond.getAtom1ID();
            if(!include_hydrogens && AtomIsHydrogen(molecule.getAtom(neighbor_index))) continue;
            neighbors.push_back(neighbor_index);
        }
    }
    return neighbors;
}

static inline bool AtomIsPolar(const UMDAtom& atom)
{
    std::string element = atom.getElement();
    return (element=="N" || element=="O" || element=="S"); // Consider N, O, and S as polar atoms
}
static inline bool AtomIsPolarHydrogen(const UMDAtom& atom, const UMDMolecule& molecule)
{
    if(!AtomIsHydrogen(atom)) return false;
    int atom_index=-1;
    for(int i=0;i<molecule.getNumAtoms();i++)
    {
        if(&molecule.getAtom(i)==&atom) {atom_index=i; break;}
    }
    if(atom_index==-1) return false; // Atom not found in molecule (shouldn't happen)
    std::vector<int> neighbors = getNeighborIndices(molecule, atom_index, true);
    for(int neighbor : neighbors)
    {
        if(AtomIsPolar(molecule.getAtom(neighbor))) return true; // If the hydrogen is bonded to a polar atom, consider it a polar hydrogen
    }
    return false;
}

static UMDBond getBondBetweenAtoms(const UMDMolecule& molecule, int atom_index1, int atom_index2)
{
    for(int i=0;i<molecule.getNumBonds();i++)
    {
        const UMDBond& bond = molecule.getBond(i);
        if((bond.getAtom1ID()==atom_index1 && bond.getAtom2ID()==atom_index2) || (bond.getAtom1ID()==atom_index2 && bond.getAtom2ID()==atom_index1)) return bond;
    }
    throw std::runtime_error("No bond found between the specified atoms");
}


static std::vector<std::pair<int,int>> generateDFSTraversalOrder(const UMDMolecule& molecule, int start_atom_index, bool include_hydrogens=true, bool include_polar_hydrogens=true) // each pair is (atom index, parent index)
{
    std::vector<std::pair<int,int>> traversal_order;
    std::vector<bool> visited(molecule.getNumAtoms(), false);
    std::function<void(int,int)> dfs = [&](int current, int parent)
    {
        visited[current]=true;
        traversal_order.push_back({current, parent});
        std::vector<int> neighbors=getNeighborIndices(molecule, current, AtomIsPolar(molecule.getAtom(current))?include_polar_hydrogens:include_hydrogens);
        std::vector<int> neighbor_order;
        for(int neighbor : neighbors) neighbor_order.push_back(-computeNumNeighbors(molecule, neighbor, AtomIsPolar(molecule.getAtom(neighbor))?include_polar_hydrogens:include_hydrogens)); // Compute the number of neighbors for each neighbor to determine the order of traversal (e.g. prioritize atoms with fewer neighbors to create branches earlier)
        std::vector<size_t> srt=argsort(neighbor_order);
        
        for(int neighbor_idx : srt)
        {
            int neighbor=neighbors[neighbor_idx];
            if(!visited[neighbor]) dfs(neighbor, current);
        }
    };
    dfs(start_atom_index, -1); // Start DFS with the specified start atom index and no parent
    return traversal_order;
}

static std::pair<std::vector<bool>,std::vector<std::vector<int>>> computeAtomRingStatus(const UMDMolecule& molecule)
{
    std::vector<bool> in_ring(molecule.getNumAtoms(), false);
    std::vector<std::pair<int,int>> tree=generateDFSTraversalOrder(molecule, 0, true, true); // Generate a DFS traversal order of the molecule starting from the first atom (index 0), including hydrogens and polar hydrogens
    std::vector<int> visited;
    std::vector<std::vector<int>> rings; // List of rings, where each ring is represented as a list of atom indices in the ring
    for(const auto& [atom_index, parent_index] : tree)
    {
        if(parent_index==-1 || AtomIsHydrogen(molecule.getAtom(atom_index))) {visited.push_back(atom_index);}
        else
        {
            while(visited.back()!=parent_index) // Pop from the visited list until we reach the parent atom of the current atom (this means we are backtracking along the DFS path)
                visited.pop_back();
            visited.push_back(atom_index);

            std::vector<int> neighbors = getNeighborIndices(molecule, atom_index, false); // Get non-hydrogen neighbors for ring detection
            for(int neighbor : neighbors)
            {
                if(neighbor==parent_index) continue; // Skip the parent atom
                if(contains(visited, neighbor)) // If the neighbor has already been visited and is not the parent, this is a ring closure
                {
                    std::vector<int> current_ring;
                    for(int i=visited.size()-1;i>=0;i--)
                    {
                        current_ring.push_back(visited[i]);
                        in_ring[visited[i]]=true; // Mark all atoms in the current path as being in a ring
                        if(visited[i]==neighbor) break;
                    }
                    rings.push_back(current_ring);
                }
            }
        }
    }
    return {in_ring, rings};
}

static std::vector<bool> computeBondRingStatus(const UMDMolecule& molecule)
{
    std::vector<bool> in_ring_bond(molecule.getNumBonds(), false);
    std::vector<int> bond_traversal;
    std::vector<std::pair<int,int>> tree=generateDFSTraversalOrder(molecule, 0, true, true); // Generate a DFS traversal order of the molecule starting from the first atom (index 0), including hydrogens and polar hydrogens
    std::vector<int> visited;
    for(const auto& [atom_index, parent_index] : tree)
    {
        if(parent_index==-1 || AtomIsHydrogen(molecule.getAtom(atom_index))) {visited.push_back(atom_index);}
        else
        {
            while(visited.back()!=parent_index) // Pop from the visited list until we reach the parent atom of the current atom (this means we are backtracking along the DFS path)
                visited.pop_back();
            visited.push_back(atom_index);

            std::vector<int> neighbors = getNeighborIndices(molecule, atom_index, false); // Get non-hydrogen neighbors for ring detection
            for(int neighbor : neighbors)
            {
                if(neighbor==parent_index) continue; // Skip the parent atom
                if(contains(visited, neighbor)) // If the neighbor has already been visited and is not the parent, this is a ring closure
                {
                    std::vector<int> current_ring;
                    for(int i=visited.size()-1;i>=0;i--)
                    {
                        current_ring.push_back(visited[i]);
                        if(visited[i]==neighbor) break;
                    }
                    for(int i=0;i<current_ring.size();i++)
                    {
                        int atom1=current_ring[i];
                        int atom2=current_ring[(i+1)%current_ring.size()];
                        for(int j=0;j<molecule.getNumBonds();j++)
                        {
                            const UMDBond& bond = molecule.getBond(j);
                            if((bond.getAtom1ID()==atom1 && bond.getAtom2ID()==atom2) || (bond.getAtom1ID()==atom2 && bond.getAtom2ID()==atom1))
                            {
                                in_ring_bond[j]=true; // Mark the bond as being in a ring
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    return in_ring_bond;
}

static std::vector<std::vector<int>> computeRigidSegments(const UMDMolecule& molecule)
{
    std::vector<std::pair<int,int>> tree=generateDFSTraversalOrder(molecule, 0, true, true); // Generate a DFS traversal order of the molecule starting from the first atom (index 0), excluding hydrogens and polar hydrogens
    std::map<int,int> parent_map; // Map of atom index to its parent index in the DFS tree
    std::vector<std::vector<int>> rigid_segments;
    std::map<int,int> h_map;

    for(const auto& [atom_index, parent_index] : tree)
    {
        if(AtomIsHydrogen(molecule.getAtom(atom_index)) && !AtomIsPolarHydrogen(molecule.getAtom(atom_index), molecule)) continue; // Skip hydrogens
        parent_map[atom_index]=parent_index;
        std::cout << "Parent of "<<atom_index<<" is "<<parent_index<<"\n";
        for(int neighbor : getNeighborIndices(molecule, atom_index, true))
        {
            if(AtomIsHydrogen(molecule.getAtom(neighbor)))
            {
                if(AtomIsPolarHydrogen(molecule.getAtom(neighbor), molecule))
                {
                    std::cout << "Neighbor " << neighbor << " is a polar hydrogen with parent " << atom_index << "\n";
                    h_map[neighbor]=atom_index; // Map polar hydrogen to its parent atom
                }
                continue;
            }
            if(neighbor==parent_index) continue; // Skip the parent atom
            if(parent_map.find(neighbor)!=parent_map.end()) // If the neighbor has already been visited and is not the parent, this is a ring closure
            {
                std::vector<int> new_segment={neighbor}; // Start the current segment with the neighbor atom (the one that closes the ring)
                std::vector<int>* current_segment=&new_segment;
                for(std::vector<int>& seg : rigid_segments) // Check if the neighbor atom is already part of an existing segment
                {
                    if(contains(seg, neighbor))
                    {
                        current_segment=&seg; // If the neighbor atom is already part of an existing segment, use that segment as the current segment
                        break;
                    }
                }
                int current_atom=atom_index;
                std::cout << "Ring closure detected between atoms " << atom_index << " and " << neighbor << ". Backtracking to find the ring path:\n";
                while(!contains(*current_segment, current_atom)) // Backtrack along the DFS tree until we reach the neighbor atom (the one that closes the ring)
                {
                    std::cout << "\t"<<current_atom;
                    current_segment->push_back(current_atom);
                    current_atom=parent_map[current_atom];
                    std::cout << " <- " << current_atom << "\n";
                }
                if(new_segment.size()>1) rigid_segments.push_back(new_segment);
            }
        }
    }
    std::cout << "No. of atoms in molecule: " << molecule.getNumAtoms() << "\n";
    for(int i=0;i<molecule.getNumAtoms();i++)
    {
        const UMDAtom& atom = molecule.getAtom(i);
        if(AtomIsHydrogen(atom)) continue; // Skip hydrogens
        bool in_existing_segment=false;
        for(std::vector<int>& seg : rigid_segments)
        {
            if(contains(seg, i)) {in_existing_segment=true; break;}
        }
        if(!in_existing_segment) rigid_segments.push_back({i}); // If the atom is not part of any existing segment, create a new segment with just that atom
        else std::cout << "Atom " << i << " is already part of an existing segment\n";
    }
    for(const auto& [h_index, parent_index] : h_map)
    {
        for(std::vector<int>& seg : rigid_segments)
        {
            if(contains(seg, parent_index)) {seg.push_back(h_index); break;}
        }
    }
    return rigid_segments;
}

static std::pair<std::string,file_pointer> getQueryByIndex(const std::string& pointer_filename, unsigned long long index)
{
    FILE* pointer_file = fopen(pointer_filename.c_str(), "rb");
    if(!pointer_file)
    {
        std::cerr << "Error opening pointer file for reading: " << pointer_filename << std::endl;
        throw NoPointerFileFoundException();
    }
    fseek(pointer_file, index * (POINTER_NAME_SIZE+sizeof(file_pointer)), SEEK_SET);
    char name_buffer[POINTER_NAME_SIZE];
    fread(name_buffer, sizeof(char), POINTER_NAME_SIZE, pointer_file);
    file_pointer position;
    fread(&position, sizeof(file_pointer), 1, pointer_file);
    fclose(pointer_file);
    return {std::string(name_buffer), position};
}

static std::map<std::string, UMDMolecule> drainUMDReader(UMDReader& reader)
{
    std::map<std::string, UMDMolecule> molecules;
    while(true)
    {
        try
        {
            UMDMolecule mol = reader.getNextMolecule();
            molecules[mol.getName()]=mol;
        }
        catch(const MoleculeDataEndedException& e) {break;}
    }
    return molecules;
}

#endif
