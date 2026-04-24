#ifndef INCLUDED_UMDFP_HELPERS_H
#define INCLUDED_UMDFP_HELPERS_H
#include "UMDFP.h"
#include <functional>
#include <sstream>

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

static UMDBond getBondBetweenAtoms(const UMDMolecule& molecule, int atom_index1, int atom_index2)
{
    for(int i=0;i<molecule.getNumBonds();i++)
    {
        const UMDBond& bond = molecule.getBond(i);
        if((bond.getAtom1ID()==atom_index1 && bond.getAtom2ID()==atom_index2) || (bond.getAtom1ID()==atom_index2 && bond.getAtom2ID()==atom_index1)) return bond;
    }
    throw std::runtime_error("No bond found between the specified atoms");
}


static std::vector<std::pair<int,int>> generateDFSTraversalOrder(const UMDMolecule& molecule, int start_atom_index, bool include_hydrogens=true) // each pair is (atom index, parent index)
{
    std::vector<std::pair<int,int>> traversal_order;
    std::vector<bool> visited(molecule.getNumAtoms(), false);
    std::function<void(int,int)> dfs = [&](int current, int parent)
    {
        visited[current]=true;
        traversal_order.push_back({current, parent});
        std::vector<int> neighbors=getNeighborIndices(molecule, current, include_hydrogens);
        std::vector<int> neighbor_order;
        for(int neighbor : neighbors) neighbor_order.push_back(computeNumNeighbors(molecule, neighbor,include_hydrogens));
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
    std::vector<std::pair<int,int>> tree=generateDFSTraversalOrder(molecule, 0, true); // Generate a DFS traversal order of the molecule starting from the first atom (index 0)
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


#endif
