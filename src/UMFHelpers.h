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

static std::vector<std::vector<float>> invertMatrix(const std::vector<std::vector<float>>& matrix)
{
    // This function inverts a 4x4 matrix using Gaussian elimination
    int n = matrix.size();
    std::vector<std::vector<float>> augmented(n, std::vector<float>(2*n, 0.0f));
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++) augmented[i][j]=matrix[i][j];
        augmented[i][n+i]=1.0f; // Append identity matrix
    }
    for(int i=0;i<n;i++)
    {
        float pivot = augmented[i][i];
        float abs_pivot = (pivot>0)? pivot : -pivot;
        if(abs_pivot<1e-6) throw std::runtime_error("Matrix is singular and cannot be inverted with pivot at row "+std::to_string(i)+" having value "+std::to_string(abs_pivot));
        for(int j=0;j<2*n;j++) augmented[i][j]/=pivot; // Normalize pivot row
        for(int k=0;k<n;k++)
        {
            if(k==i) continue;
            float factor = augmented[k][i];
            for(int j=0;j<2*n;j++) augmented[k][j]-=factor*augmented[i][j]; // Eliminate column
        }
    }
    std::vector<std::vector<float>> inverse(n, std::vector<float>(n, 0.0f));
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++) inverse[i][j]=augmented[i][n+j]; // Extract inverse matrix
    }
    return inverse;
}
static std::vector<std::vector<float>> multiplyMatrices(const std::vector<std::vector<float>>& A, const std::vector<std::vector<float>>& B)
{
    int n = A.size();
    int m = B[0].size();
    int p = B.size();
    if(A[0].size()!=p) throw std::runtime_error("Matrix dimensions do not match for multiplication");
    std::vector<std::vector<float>> C(n, std::vector<float>(m, 0.0f));
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            for(int k=0;k<p;k++) C[i][j]+=A[i][k]*B[k][j];
        }
    }
    return C;
}
static std::vector<std::vector<std::vector<float>>> localAlignAtoms(std::vector<std::vector<float>> old_coords, std::vector<std::vector<float>> new_coords, std::vector<bool> use_idxs, const UMDMolecule& molecule)
{
    std::vector<std::vector<std::vector<float>>> ret;
    for(int i=0;i<use_idxs.size();i++)
    {
        if(use_idxs[i]) {ret.push_back({}); continue;} // Don't transform the coordinates that are used for alignment; transform the rest
        std::vector<int> neighbors = getNeighborIndices(molecule, i, true);
        std::vector<int> align_atoms;
        int qidx=0;
        constexpr int N_NEIGH=4;
        while(align_atoms.size()<N_NEIGH)
        {
            for(int neighbor : neighbors)
            {
                if(use_idxs[neighbor] && !contains(align_atoms, neighbor)) align_atoms.push_back(neighbor);
                if(align_atoms.size()>=N_NEIGH) break;
            }
            if(align_atoms.size()>=N_NEIGH) break;
            if(align_atoms.size()==0) break; // If no neighbors are available for alignment, break
            neighbors=getNeighborIndices(molecule, align_atoms[qidx++], true);
        }
        if(align_atoms.size()<N_NEIGH) {ret.push_back({}); continue;} // If less than 4 atoms are available for alignment, skip this atom
        
        std::vector<std::vector<float>> old_align_coords, new_align_coords;
        for(int idx : align_atoms)
        {
            old_align_coords.push_back(old_coords[idx]);
            new_align_coords.push_back(new_coords[idx]);
        }
        std::vector<std::vector<float>> transform; // 4x4 transformation matrix
        for(int i=0;i<4;i++) transform.push_back(std::vector<float>(4, 0.0f));

        // Compute the transformation matrix by matrix inversion using old and new coordinates
        std::vector<std::vector<float>> old_matrix(4, std::vector<float>(4, 0.0f));
        std::vector<std::vector<float>> new_matrix(4, std::vector<float>(4, 0.0f));
        for(int i=0;i<N_NEIGH;i++)
        {
            old_matrix[i][0]=old_align_coords[i][0];
            old_matrix[i][1]=old_align_coords[i][1];
            old_matrix[i][2]=old_align_coords[i][2];
            old_matrix[i][3]=1.0f;
            new_matrix[i][0]=new_align_coords[i][0];
            new_matrix[i][1]=new_align_coords[i][1];
            new_matrix[i][2]=new_align_coords[i][2];
            new_matrix[i][3]=1.0f;
        }
        old_matrix[3][3]=1.0f; // Set the last row of the old matrix to [0, 0, 0, 1]
        new_matrix[3][3]=1.0f; // Set the last row of the new matrix to [0, 0, 0, 1]

        // Compute the transformation matrix as inverse(old_matrix)*new_matrix
        try{
            std::vector<std::vector<float>> old_matrix_inv = invertMatrix(old_matrix);
            transform = multiplyMatrices(old_matrix_inv, new_matrix);
            ret.push_back(transform);
        }
        catch(const std::runtime_error& e)
        {
            std::cerr << "Error inverting matrix for atom " << i << ": " << e.what() << std::endl;
            ret.push_back({}); // If the matrix is singular, skip this atom
        }
    }
    return ret;
}

static void minimizeTerminalAngles(UMDMolecule& molecule, float tolerance=5.0f, float step_size=0.1f, std::vector<bool> edit_frozen=std::vector<bool>())
{
    std::vector<bool> allow_movement(molecule.getNumAtoms(), true);
    for(int i=0;i<molecule.getNumAtoms();i++)
    {
        if(computeNumNeighbors(molecule, i, true)>1) allow_movement[i]=false; // Only move terminal atoms (with only one non-hydrogen neighbor)
        if(!edit_frozen.empty() && edit_frozen[i]) allow_movement[i]=false; // If edit_limited is provided, only allow movement for atoms marked as true
    }
    for(int i=0;i<molecule.getNumAtoms();i++)
    {
        float TARGET_ANGLE=109.5f; // Target angle for SP3 hybridized atoms
        if(molecule.getAtom(i).getHybridization()==3) // SP3
            TARGET_ANGLE=109.5f;
        else if(molecule.getAtom(i).getHybridization()==2) // SP2
            TARGET_ANGLE=120.0f;
        else if(molecule.getAtom(i).getHybridization()==6) // SP3D2
            TARGET_ANGLE=90.0f;
        else continue;

        std::vector<int> neighbors = getNeighborIndices(molecule, i, true);
        if(neighbors.size()<3) continue; // Need at least 3 neighbors to define an angle
        std::vector<std::vector<float>> forces;
        bool forces_set=true;
        int iter_num=0;
        while(iter_num<1000 && forces_set)
        {
            forces_set=false;
            forces=std::vector<std::vector<float>>(); // Reset forces
            for(int j=0;j<neighbors.size();j++) forces.push_back({0.0f, 0.0f, 0.0f}); // Initialize forces to zero 

            for(int j=0;j<neighbors.size();j++)
            {
                std::vector<float> force={0.0f, 0.0f, 0.0f};
                for(int k=0;k<neighbors.size();k++)
                {
                    if(j==k) continue; // Skip the same neighbor
                    int neighbor1=neighbors[j];
                    int neighbor2=neighbors[k];
                    //if(!allow_movement[neighbor1] && !allow_movement[neighbor2]) continue; // Only adjust angles for terminal atoms

                    std::vector<float> vec1 = {molecule.getAtom(neighbor1).getX()-molecule.getAtom(i).getX(), molecule.getAtom(neighbor1).getY()-molecule.getAtom(i).getY(), molecule.getAtom(neighbor1).getZ()-molecule.getAtom(i).getZ()};
                    std::vector<float> vec2 = {molecule.getAtom(neighbor2).getX()-molecule.getAtom(i).getX(), molecule.getAtom(neighbor2).getY()-molecule.getAtom(i).getY(), molecule.getAtom(neighbor2).getZ()-molecule.getAtom(i).getZ()};
                    float dot_product = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
                    float mag1 = std::sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]);
                    float mag2 = std::sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2]);
                    float angle = std::acos(dot_product/(mag1*mag2))*180.0f/M_PI; // Angle in degrees
                    if(angle<TARGET_ANGLE-tolerance || angle>TARGET_ANGLE+tolerance)
                    {
                        // Compute the force vector to adjust the angle
                        std::vector<float> middle_vector = {vec1[0]+vec2[0], vec1[1]+vec2[1], vec1[2]+vec2[2]};
                        float middle_mag = std::sqrt(middle_vector[0]*middle_vector[0]+middle_vector[1]*middle_vector[1]+middle_vector[2]*middle_vector[2]);
                        if(middle_mag<1e-6) continue; // Skip if the cross product is too small (vectors are parallel)
                        for(int d=0;d<3;d++) middle_vector[d]/=middle_mag; // Normalize
                        float angle_diff = (angle<TARGET_ANGLE) ? (TARGET_ANGLE-angle) : (TARGET_ANGLE-angle);
                        std::vector<float> unproj_vector;
                        for(int d=0;d<3;d++) unproj_vector.push_back(middle_vector[d]-((middle_vector[0]*vec1[0]+middle_vector[1]*vec1[1]+middle_vector[2]*vec1[2])/(mag1*mag1))*vec1[d]);
                        for(int d=0;d<3;d++)
                            force[d]-=unproj_vector[d]*angle_diff*step_size; // Scale the force by the angle difference and step size

                        /*std::vector<float> off_plane={vec1[1]*vec2[2]-vec1[2]*vec2[1], vec1[2]*vec2[0]-vec1[0]*vec2[2], vec1[0]*vec2[1]-vec1[1]*vec2[0]};
                        float off_plane_mag = std::sqrt(off_plane[0]*off_plane[0]+off_plane[1]*off_plane[1]+off_plane[2]*off_plane[2]);
                        if(off_plane_mag>1e-6)
                        {
                            for(int d=0;d<3;d++) off_plane[d]/=off_plane_mag; // Normalize
                            for(int d=0;d<3;d++)
                                force[d]+=off_plane[d]*step_size; // Scale the force by the angle difference and step size
                        }*/
                        /*for(int d=0;d<3;d++)
                            forces[k][d]-=unproj_vector[d]*angle_diff*step_size; // Scale the force by the angle difference and step size*/
                    }
                }
                // Apply the computed force to the neighbor atom if it is allowed to move
                if(allow_movement[neighbors[j]])
                {
                    for(int d=0;d<3;d++)
                        forces[j][d]+=force[d]; // Accumulate the force for this neighbor
                    if(force[0]*force[0]+force[1]*force[1]+force[2]*force[2]>1e-6) forces_set=true; // If the force is significant, set forces_set to true
                }
                else forces.push_back({0.0f, 0.0f, 0.0f}); // If the atom is not allowed to move, set the force to zero
            }
            // Apply the forces to the neighbor atoms
            for(int j=0;j<neighbors.size();j++)
            {
                if(allow_movement[neighbors[j]])
                {
                    molecule.getAtom(neighbors[j]).getData().x += forces[j][0];
                    molecule.getAtom(neighbors[j]).getData().y += forces[j][1];
                    molecule.getAtom(neighbors[j]).getData().z += forces[j][2];
                }
            }
            iter_num++;
        }
        /*for(int j=0;j<neighbors.size();j++)
        {
            int neighbor1=neighbors[j];
            int neighbor2=neighbors[(j+1)%neighbors.size()];
            if(!allow_movement[neighbor2]) continue; // Only adjust angles for terminal atoms
            std::vector<float> vec1 = {molecule.getAtom(neighbor1).getX()-molecule.getAtom(i).getX(), molecule.getAtom(neighbor1).getY()-molecule.getAtom(i).getY(), molecule.getAtom(neighbor1).getZ()-molecule.getAtom(i).getZ()};
            std::vector<float> vec2 = {molecule.getAtom(neighbor2).getX()-molecule.getAtom(i).getX(), molecule.getAtom(neighbor2).getY()-molecule.getAtom(i).getY(), molecule.getAtom(neighbor2).getZ()-molecule.getAtom(i).getZ()};
            float dot_product = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
            float mag1 = std::sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]);
            float mag2 = std::sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2]);
            float angle = std::acos(dot_product/(mag1*mag2))*180.0f/M_PI; // Angle in degrees
            while(angle<TARGET_ANGLE-tolerance || angle>TARGET_ANGLE+tolerance)
            {
                std::cout << "Adjusting angle between atoms " << neighbor1 << ", " << i << ", " << neighbor2 << " from " << angle << " degrees to " << TARGET_ANGLE << " degrees\n";
                // Adjust the position of neighbor2 to achieve the target angle (this is a simple adjustment and may not be physically accurate)
                float scale = mag2*std::tan(TARGET_ANGLE*(M_PI/180.0f-0.001))/mag1;
                molecule.getAtom(neighbor2).getData().x = molecule.getAtom(i).getX() - vec1[0]*scale*step_size;
                molecule.getAtom(neighbor2).getData().y = molecule.getAtom(i).getY() - vec1[1]*scale*step_size;
                molecule.getAtom(neighbor2).getData().z = molecule.getAtom(i).getZ() - vec1[2]*scale*step_size;

                // Recompute the angle after adjustment
                vec2 = {molecule.getAtom(neighbor2).getX()-molecule.getAtom(i).getX(), molecule.getAtom(neighbor2).getY()-molecule.getAtom(i).getY(), molecule.getAtom(neighbor2).getZ()-molecule.getAtom(i).getZ()};
                dot_product = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
                mag1 = std::sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]);
                mag2 = std::sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2]);
                angle = std::acos(dot_product/(mag1*mag2))*180.0f/M_PI; // Angle in degrees
            }
        }*/
    }
}
#endif
