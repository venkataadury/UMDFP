#ifndef UMDFP_H
#define UMDFP_H 1
#define _FILE_OFFSET_BITS 64
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <string>
#include <gzstream.h>
#include <cstdint>
#include <algorithm>

typedef fpos_t file_pointer;


// Add a byte with ascii value 244 to end of header as the 'magic number' to avoid corruption
static const char* UMFHeader = "\x1eUMF1\n"; // UMF version 1.0

// A compresser to take Universal Molecular Description (UMD) files to compress them into easy-to-access Universal Molecular Format (UMF) files.
// Each UMF file can contain over 1 billion molecules and is designed to be read in a streaming fashion, so that it can be used on machines with limited memory.
class ByteString
{
    unsigned char* data;
    size_t length;
public:
    ByteString() : data(nullptr), length(0) {}
    ByteString(size_t length) : length(length)
    {
        data = new unsigned char[length];
    }
    ByteString(const char* str, size_t length) : length(length)
    {
        data = new unsigned char[length];
        memcpy(data, str, length);
    }
    template<typename T> ByteString(const T& obj)
    {
        length = sizeof(T);
        data = new unsigned char[length];
        memcpy(data, &obj, length);
    }
    ByteString(const ByteString& other) : length(other.length)
    {
        data = new unsigned char[length];
        memcpy(data, other.data, length);
    }
    ~ByteString()
    {
        delete[] data;
    }

    inline unsigned char* getData() { return data; }
    inline const unsigned char* getData() const { return data; }
    inline size_t getLength() const { return length; }
    inline void writeToFile(FILE* file) const {fwrite(data, sizeof(unsigned char), length, file);}

    ByteString& operator=(const ByteString& other)
    {
        if (this != &other)
        {
            delete[] data;
            length = other.length;
            data = new unsigned char[length];
            memcpy(data, other.data, length);
        }
        return *this;
    }
};

template<typename T> static std::vector<size_t> argsort(const std::vector<T>& vec)
{
    std::vector<size_t> indices(vec.size());
    for (size_t i = 0; i < vec.size(); i++) indices[i] = i;
    std::sort(indices.begin(), indices.end(), [&vec](size_t a, size_t b) { return vec[a] < vec[b]; });
    return indices;
}

struct UMDAtomData
{
    float x, y, z; // 3D coordinates
    char element[3]; // Element symbol (e.g., C, H, O)
    float charge; // Charge of the atom
    int hybridization; // Hybridization state (e.g., sp, sp2, sp3)
    bool aromatic; // Whether the atom is aromatic
};

class UMDAtom
{
protected:
    UMDAtomData data;

public:
    UMDAtom(std::string line)
    {
        // Parse the line to extract atom information
        int index;
        int arom;
        double posx, posy, posz;
        sscanf(line.c_str(), "%d %s %lf %lf %lf %f %d %d",&index, data.element, &posx, &posy, &posz, &data.charge, &data.hybridization, &arom);
        data.x = posx;
        data.y = posy;
        data.z = posz;
        data.aromatic = (arom != 0);
    }
    UMDAtom(UMDAtomData data) : data(data) {}
    UMDAtom(const ByteString& bs) {this->fromByteString(bs);}

    inline const UMDAtomData& getData() const {return data;}
    inline UMDAtomData& getData()  {return data;}
    inline ByteString toByteString() const {return ByteString(data);}
    inline void fromByteString(ByteString bs) {memcpy(&data, bs.getData(), bs.getLength());}

    inline std::string getElement() const {return std::string(data.element);} // Get element symbol as a string for easier use
    inline float getX() const {return data.x;}
    inline float getY() const {return data.y;}
    inline float getZ() const {return data.z;}
    inline float getCharge() const {return data.charge;}
    inline int getHybridization() const {return data.hybridization;}
    inline bool isAromatic() const {return data.aromatic;}
};

struct UMDBondData
{
    int atom1, atom2; // Indices of the two atoms in the bond
    int bond_type; // Type of bond (e.g., single, double, triple, aromatic represented as 2, 4, 6, 3 respectively)  
};

class UMDBond
{
protected:
    UMDBondData data;
public:
    UMDBond(std::string line)
    {
        // Parse the line to extract bond information
        int index;
        sscanf(line.c_str(), "%d %d %d %d",&index, &data.atom1, &data.atom2, &data.bond_type);
    }
    UMDBond(UMDBondData data) : data(data) {}
    UMDBond(const ByteString& bs) {this->fromByteString(bs);}

    inline const UMDBondData& getData() const {return data;}
    inline UMDBondData& getData()  {return data;}
    inline ByteString toByteString() const {return ByteString(data);}
    inline void fromByteString(ByteString bs) {memcpy(&data, bs.getData(), bs.getLength());}
    inline int getAtom1ID() const {return data.atom1;}
    inline int getAtom2ID() const {return data.atom2;}
    inline int getBondType() const {return data.bond_type;}
};

class UMDMolecule
{
protected:
    std::vector<UMDAtom> atoms;
    std::vector<UMDBond> bonds;
    std::string name;
    ByteString smiles;
    ByteString extras; // For any additional data that may be needed in the future, such as stereochemistry, isotopes, etc.
    int lookahead_256_size=0; // Size of the next 256 molecules in bytes (for look-ahead)
    unsigned long lookahead_65536_size=0; // Size of the next 65536 molecules in bytes (for look-ahead)

private:
    std::string _nextReleventLine(std::vector<std::string>& lines, unsigned long& index, bool allow_end=false)
    {
        index++;
        for(;index<lines.size();index++)
        {
            // Remove trailing whitespaces
            lines[index].erase(lines[index].find_last_not_of(" \n\r\t")+1);
            if(lines[index].empty()) continue;
            if(lines[index][0]==';') continue; // Skip comment lines
            if(lines[index]=="END") {
                if(allow_end) return lines[index];
                throw std::runtime_error("Unexpected END line while parsing molecule data");
            }
            return lines[index];
        }
        std::cerr << "Error: Unexpected end of molecule data" << std::endl;
        throw std::runtime_error("Invalid molecule data");
    }

public:
    UMDMolecule() {}
    UMDMolecule(std::vector<std::string> lines)
    {
        unsigned long i=0;
        for(;i<lines.size();i++)
        {
            // Remove trailing whitespaces
            lines[i].erase(lines[i].find_last_not_of(" \n\r\t")+1);
            if(lines[i].empty()) continue;
            if(lines[i][0]==';') continue; // Skip comment lines
            if(lines[i]=="START") break; // Wait for the START line to begin parsing
        }
        if(i==lines.size())
        {
            std::cerr << "Error: No START line found in molecule data" << std::endl;
            throw std::runtime_error("Invalid molecule data");
        }
        name=this->_nextReleventLine(lines, i);
        std::string smiles_line=this->_nextReleventLine(lines, i);
        smiles = ByteString(smiles_line.c_str(), smiles_line.size()+1); // Include null terminator
        
        int natoms=std::stoi(this->_nextReleventLine(lines, i));
        int nbonds=std::stoi(this->_nextReleventLine(lines, i));
        for(int a=0;a<natoms;a++) atoms.push_back(UMDAtom(this->_nextReleventLine(lines, i)));
        for(int b=0;b<nbonds;b++) bonds.push_back(UMDBond(this->_nextReleventLine(lines, i)));

        // Parse any additional data
        std::string remaining_data;
        std::string line;
        while((line=this->_nextReleventLine(lines, i, true))!="END") remaining_data+=line+"\n";
        extras=ByteString(remaining_data.c_str(), remaining_data.size());
    }

    int computeWrittenSize() const
    {
        // Size of header
        // Size of this molecule block
        // Size of next 256 molecules (look-ahead)
        // Size of next 65536 molecules (look-ahead)
        // Size of int * 2 (# atoms, # bonds)
        // Size of name (Include '\0')
        // Size of SMILES (Include '\0')
        // Size of atoms
        // Size of bonds
        // Size of extras
        
        return sizeof(char)*strlen(UMFHeader) + sizeof(int)*2+ sizeof(unsigned long) + 2*sizeof(int) + name.size()*sizeof(char) + smiles.getLength() + atoms.size()*sizeof(UMDAtomData) + bonds.size()*sizeof(UMDBondData) + extras.getLength();
    }

    int toFile(FILE* file) const
    {
        // Write header
        fwrite(UMFHeader, sizeof(char), strlen(UMFHeader), file);
        int mySize = this->computeWrittenSize();
        fwrite(&mySize, sizeof(int), 1, file);
        fwrite(&lookahead_256_size, sizeof(int), 1, file);
        fwrite(&lookahead_65536_size, sizeof(unsigned long), 1, file);

        int natoms=atoms.size();
        int nbonds=bonds.size();
        fwrite(&natoms, sizeof(int), 1, file);
        fwrite(&nbonds, sizeof(int), 1, file);

        // Get data sizes
        size_t name_size = name.size();
        size_t smiles_size = smiles.getLength();
        size_t atoms_size = atoms.size() * sizeof(UMDAtomData);
        size_t bonds_size = bonds.size() * sizeof(UMDBondData);
        size_t extras_size = extras.getLength();

        // Write data
        fwrite(name.c_str(), sizeof(char), name.size()+1, file);
        smiles.writeToFile(file);

        for(int i=0;i<atoms.size();i++) atoms[i].toByteString().writeToFile(file);
        for(int i=0;i<bonds.size();i++) bonds[i].toByteString().writeToFile(file);

        extras.writeToFile(file);
        return mySize;
    }

    void summary() const
    {
        std::cout << "Molecule: " << name << ", SMILES: " << std::string((char*)smiles.getData()) << ", Atoms: " << atoms.size() << ", Bonds: " << bonds.size() << std::endl;
        for(int i=0;i<atoms.size();i++)
        {
            const UMDAtomData& atom_data = atoms[i].getData();
            std::cout << "  Atom " << i << ": " << atom_data.element << " (" << atom_data.x << ", " << atom_data.y << ", " << atom_data.z << "), Charge: " << atom_data.charge << ", Hybridization: " << atom_data.hybridization << ", Aromatic: " << (atom_data.aromatic ? "Yes" : "No") << std::endl;
        }
        for(int i=0;i<bonds.size();i++)        {
            const UMDBondData& bond_data = bonds[i].getData();
            std::cout << "  Bond " << i << ": Atoms " << bond_data.atom1 << " - " << bond_data.atom2 << ", Type: " << bond_data.bond_type << std::endl; 
        }
    }

    inline const std::string& getName() const {return name;}
    inline const ByteString& getSMILES() const {return smiles;}
    inline int getNumAtoms() const {return atoms.size();}
    inline int getNumBonds() const {return bonds.size();}

    inline const UMDAtom& getAtom(int index) const {return atoms[index];}
    inline UMDAtom& getAtom(int index) {return atoms[index];}
    inline const UMDBond& getBond(int index) const {return bonds[index];}
    inline UMDBond& getBond(int index) {return bonds[index];}
};

#define ATOM_DATA_SIZE sizeof(UMDAtomData)
#define BOND_DATA_SIZE sizeof(UMDBondData)

class MoleculeDataEndedException : public std::runtime_error
{
public:    MoleculeDataEndedException() : std::runtime_error("No more molecule data to read") {}
};

class UMDReader
{
    std::ifstream file;
public:
    UMDReader(std::string filename) : file(filename, std::ios::binary)
    {
        if (!file)
        {
            std::cerr << "Error opening file for reading: " << filename << std::endl;
            throw std::runtime_error("Could not open file");
        }
    }

    UMDMolecule getNextMolecule()
    {
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(file, line))
        {
            if(line.empty()) continue;
            if(line[0]==';') continue; // Skip comment lines
            // Remove trailing whitespaces
            line.erase(line.find_last_not_of(" \n\r\t")+1);
            lines.push_back(line);
            if(line=="END") break; // End of molecule block
        }
        if(!lines.size()) throw MoleculeDataEndedException();
        if(lines[0]!="START")
        {
            std::cerr << "Error: Expected START line at the beginning of molecule data" << std::endl;
            throw std::runtime_error("Invalid molecule data");
        }
        return UMDMolecule(lines);
    }
};

template<class T, int N> class CyclicArray
{
public:
    T data[N];
    int index=0;

    int size() const {return N;}
    void push(const T& value)
    {
        data[index]=value;
        this->step();
    }
    inline void step() {index=(index+1)%N;}
    const T& get(int i) const {return data[(index+i)%N];}
    T& get(int i) {return data[(index+i)%N];}

    const T& getCurrent() const {return data[(index-1+N)%N];}
    T& getCurrent() {return data[(index-1+N)%N];}

    const T& operator[](int i) const {return this->get(i);}
    T& operator[](int i) {return this->get(i);}
};

class UMFWriter
{
    FILE* file;
    CyclicArray<int, 256> lookahead_256_sizes; // Store the sizes of the next 256 molecules for look-ahead
    CyclicArray<unsigned long, 65536> lookahead_65536_sizes; // Store the sizes of the next 65536 molecules for look-ahead
    CyclicArray<file_pointer, 65536> molecule_positions; // Store the file positions of the last 65536 molecules for potential future use (e.g., seeking back to them if needed)
    int molnum=0; // Number of molecules written so far, used to determine when to update look-ahead sizes
public:
    UMFWriter(std::string filename)
    {
        file = fopen(filename.c_str(), "wb");
        if (!file)
        {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            throw std::runtime_error("Could not open file");
        }
    }
    ~UMFWriter()
    {
        if (file)
            fclose(file);
    }

    void writeMolecule(const UMDMolecule& molecule, bool build_lookahead=true)
    {
        fgetpos(file,&(molecule_positions[0])); // Store the file position of the current molecule for potential future use (e.g., seeking back to it if needed)
        molecule_positions.step();
        int molsize=molecule.toFile(file);
        molnum++;
        // Update look-ahead sizes
        if(build_lookahead)
        {
            lookahead_256_sizes.push(molsize);
            lookahead_65536_sizes.push(molsize);
            if(molnum>=256)
            {
                int total_256_size=0;
                file_pointer current_pos;
                fgetpos(file, &current_pos);
                for(int i=0;i<256;i++) total_256_size+=lookahead_256_sizes[i];
                // Update the look-ahead size for the next 256 molecules in the header of the molecule that was written 256 molecules ago
                file_pointer pos = molecule_positions.get(-256);
                fsetpos(file, &pos);
                fseek(file, sizeof(char)*strlen(UMFHeader) + sizeof(int), SEEK_CUR); // Move to the position of the look-ahead size in the header
                fwrite(&total_256_size, sizeof(int), 1, file);
                fsetpos(file, &current_pos); // Move back to the current position to continue
            }
            if(molnum>=65536)
            {
                unsigned long total_65536_size=0;
                file_pointer current_pos;
                fgetpos(file, &current_pos);
                for(int i=0;i<65536;i++) total_65536_size+=lookahead_65536_sizes[i];
                // Update the look-ahead size for the next 65536 molecules in the header of the molecule that was written 65536 molecules ago
                file_pointer pos = molecule_positions.get(0);
                fsetpos(file, &pos);
                fseek(file, sizeof(char)*strlen(UMFHeader) + sizeof(int)*2, SEEK_CUR); // Move to the position of the look-ahead size in the header
                fwrite(&total_65536_size, sizeof(unsigned long), 1, file);
                fsetpos(file, &current_pos); // Move back to the current position to continue
            }
        }
    }
    void flush() const {fflush(file);}
};

static void push_string(char* str, const size_t& n)
{
    for(int i=0;i<n-1;i++)
        str[i]=str[i+1];
    str[n-1]='\0';
}

class NoHeaderRemainingException : public std::runtime_error
{
public:
    NoHeaderRemainingException() : std::runtime_error("No UMF header remaining in file") {}
};
class NoPointerFileFoundException : public std::runtime_error
{public:
    NoPointerFileFoundException() : std::runtime_error("No pointer file found for UMF file") {}
};
class NoSuchMoleculeException : public std::runtime_error
{
public:
    NoSuchMoleculeException() : std::runtime_error("No molecule matching the query was found in the UMF file") {}
};

struct CurrentHeader
{
    int size;
    int lookahead_256_size;
    unsigned long lookahead_65536_size;
    int natoms;
    int nbonds;
    std::string name;
    std::string smiles;
};
class UMFReader
{
    FILE* file;
    file_pointer ref;
    CurrentHeader current_header;
    bool ended=false;
    std::string filename;
    bool has_pointer=false;
    FILE* pointer_file;
public:
    UMFReader(std::string filename, const std::string& pointer_filename="")
    {
        this->filename=filename;
        file = fopen(filename.c_str(), "rb");
        if (!file)
        {
            std::cerr << "Error opening file for reading: " << filename << std::endl;
            throw std::runtime_error("Could not open file");
        }
        this->_findPointerFile(pointer_filename);
        this->_verifyAndReadHeader();
        ended=false;
    }
    ~UMFReader()
    {
        if (file)
            fclose(file);
    }
private:
    void _findPointerFile(std::string pointer_filename)
    {
        if(pointer_filename.empty()) pointer_filename=filename+"p";
        pointer_file = fopen(pointer_filename.c_str(), "rb");
        if(!pointer_file)
        {
            std::cout << "NOTE: No pointer file found. Quick access will be restricted\n";
            has_pointer=false;
            return; // If pointer file doesn't exist, just start reading from the beginning of the UMF file and hope for the best (i.e., that the desired molecule is near the beginning of the file). This allows the UMFReader to still function even without the pointer file, albeit with potentially much slower access for molecules that are far into the file.
        }
        has_pointer=true;
    }
    void _verifyAndReadHeader()
    {
        char header[7];
        fread(header, sizeof(char), 6, file);
        header[6] = '\0';
        bool warned=false;
        while (strcmp(header, UMFHeader) != 0)
        {
            if (feof(file)) throw NoHeaderRemainingException();
            
            if(!warned)
            {
                std::cerr << "Warning: UMF header not found at expected position, attempting to resync..." << std::endl;
                std::cerr << "Read header: '" << header << "'" << std::endl;
                warned=true;
            }
            // Shift the header by one byte and read the next byte
            push_string(header, 7);
            fread(header + 5, sizeof(char), 1, file);
        }
        fgetpos(file, &ref); // Store the position of the start of the molecule block for potential future use (e.g., seeking back to it if needed)
        
        // Read the sizes of the molecule block and look-ahead blocks
        fread(&current_header.size, sizeof(int), 1, file);
        fread(&current_header.lookahead_256_size, sizeof(int), 1, file);
        fread(&current_header.lookahead_65536_size, sizeof(unsigned long), 1, file);
        
        // Read the number of atoms and bonds
        fread(&current_header.natoms, sizeof(int), 1, file);
        fread(&current_header.nbonds, sizeof(int), 1, file);

        //Read name character-by-character until '\0' is found (since name can be of variable length)
        char c;
        current_header.name="";
        while(fread(&c, sizeof(char), 1, file) == 1 && c != '\0') {current_header.name+=c;}

        // Similarly, read the SMILES string character-by-character until '\0' is found (since it can also be of variable length)
        current_header.smiles="";
        while(fread(&c, sizeof(char), 1, file) == 1 && c != '\0') {current_header.smiles+=c;}
    }

public:
    inline int getNextMoleculeSize() const {return current_header.size;}
    inline int getNext256MoleculesSize() const {return current_header.lookahead_256_size;}
    inline unsigned long getNext65536MoleculesSize() const {return current_header.lookahead_65536_size;}
    inline std::string getNextMoleculeName() const {return current_header.name;}
    inline std::string getNextMoleculeSMILES() const {return current_header.smiles;}

    void dumpHeader()
    {
        std::cout << "Next molecule size: " << current_header.size << " bytes" << std::endl;
        std::cout << "Next 256 molecules size: " << current_header.lookahead_256_size << " bytes" << std::endl;
        std::cout << "Next 65536 molecules size: " << current_header.lookahead_65536_size << " bytes" << std::endl;
        std::cout << "Next molecule name: " << current_header.name << std::endl;
        std::cout << "Next molecule SMILES: " << current_header.smiles << std::endl;
    }

    void skipMolecule(unsigned long n=1)
    {
        for(unsigned long i=0;i<n;)
        {
            // Move the file pointer back to 'ref'
            fsetpos(file, &ref); // Move to the position stored in 'ref'
            int offset=sizeof(char)*strlen(UMFHeader);
            if((n-i)>=65536 && current_header.lookahead_65536_size>0)
            {
                fseek(file, current_header.lookahead_65536_size-offset+65536, SEEK_CUR); // Skip the next 65536 molecules
                i+=65536;
            }
            else if((n-i)>=256 && current_header.lookahead_256_size>0)
            {
                fseek(file, current_header.lookahead_256_size-offset+256, SEEK_CUR); // Skip the next 256 molecules
                i+=256;
            }
            else
            {
                fseek(file, current_header.size-offset+1, SEEK_CUR); // Skip the current molecule
                i++;
            }

            try{this->_verifyAndReadHeader();} // Read the header of the next molecule
            catch(const NoHeaderRemainingException& e)
             {
                ended=true;
                break;
            }
        }
    }

    UMDMolecule readMolecule()
    {
        // Read the molecule data based on the information in the header
        std::vector<std::string> lines;
        lines.push_back("START");
        lines.push_back(current_header.name);
        lines.push_back(current_header.smiles);
        lines.push_back(std::to_string(current_header.natoms));
        lines.push_back(std::to_string(current_header.nbonds));

        // Read atoms
        for(int i=0;i<current_header.natoms;i++)
        {
            UMDAtomData atom_data;
            fread(&atom_data, sizeof(UMDAtomData), 1, file);
            UMDAtom atom(atom_data);
            lines.push_back(std::to_string(i)+" "+std::string(atom.getData().element)+" "+std::to_string(atom.getData().x)+" "+std::to_string(atom.getData().y)+" "+std::to_string(atom.getData().z)+" "+std::to_string(atom.getData().charge)+" "+std::to_string(atom.getData().hybridization)+" "+std::to_string(atom.getData().aromatic));
        }

        // Read bonds
        for(int i=0;i<current_header.nbonds;i++)
        {
            UMDBondData bond_data;
            fread(&bond_data, sizeof(UMDBondData), 1, file);
            UMDBond bond(bond_data);
            lines.push_back(std::to_string(i)+" "+std::to_string(bond.getData().atom1)+" "+std::to_string(bond.getData().atom2)+" "+std::to_string(bond.getData().bond_type));
        }

        // Read extras
        int remaining_bytes = current_header.size - sizeof(char)*strlen(UMFHeader) - sizeof(int)*2 - sizeof(unsigned long) - 2*sizeof(int) - (current_header.name.size()+1)*sizeof(char) - current_header.smiles.size() - current_header.natoms*sizeof(UMDAtomData) - current_header.nbonds*sizeof(UMDBondData);
        
        char* extras_data = new char[remaining_bytes];
        fread(extras_data, sizeof(char), current_header.size - sizeof(char)*strlen(UMFHeader) - sizeof(int)*2 - sizeof(unsigned long) - 2*sizeof(int) - (current_header.name.size()+1)*sizeof(char) - current_header.smiles.size() - current_header.natoms*sizeof(UMDAtomData) - current_header.nbonds*sizeof(UMDBondData), file);
        std::string extras_str(extras_data, current_header.size - sizeof(char)*strlen(UMFHeader) - sizeof(int)*2 - sizeof(unsigned long) - 2*sizeof(int) - (current_header.name.size()+1)*sizeof(char) - current_header.smiles.size() - current_header.natoms*sizeof(UMDAtomData) - current_header.nbonds*sizeof(UMDBondData));
        lines.push_back(extras_str);
        delete[] extras_data;
        lines.push_back("END");
        
        try{this->_verifyAndReadHeader();} // Read the header of the next molecule for future reads
        catch(const NoHeaderRemainingException& e)
        {
            ended=true;
        }
        return UMDMolecule(lines);
    }

    inline bool hasEnded() const {return ended;}
    inline const CurrentHeader& getHeader() const {return current_header;}
    inline file_pointer getCurrentMoleculePosition() const {return ref;}

    file_pointer query_by_name(const std::string& molname)
    {
        if(!has_pointer) throw NoPointerFileFoundException();
        // Perform a binary search on the pointer file to find the file position of the molecule with the given name
        fseek(pointer_file, 0, SEEK_END);
        long num_entries = ftell(pointer_file) / (64*sizeof(char) + sizeof(file_pointer));
        long left=0, right=num_entries-1;
        char name_buffer[64];
        file_pointer position;
        while(left<=right)
        {
            long mid = left + (right-left)/2;
            fseek(pointer_file, mid*(64*sizeof(char) + sizeof(file_pointer)), SEEK_SET);
            fread(name_buffer, sizeof(char), 64, pointer_file);
            fread(&position, sizeof(file_pointer), 1, pointer_file);
            std::string current_name(name_buffer);
            if(current_name==molname) return position;
            else if(current_name<molname) left=mid+1;
            else right=mid-1;
        }
        throw NoSuchMoleculeException();
    }

    inline void jumpToMoleculeAtPosition(file_pointer position)
    {
        fsetpos(file, &position); // Move to the position of the header of the molecule block
        fseek(file, -sizeof(char)*strlen(UMFHeader), SEEK_CUR); // Move back to the position of the UMF header for that molecule
        try{this->_verifyAndReadHeader();} // Read the header of the molecule at the new position for future reads
        catch(const NoHeaderRemainingException& e)
        {
            ended=true;
        }
    }
};


static std::vector<std::string> readTextFileToStringList(const std::string& filename)
{
    std::vector<std::string> lines;
    std::ifstream infile(filename);
    if (!infile)
    {
        std::cerr << "Error opening input file: " << filename << std::endl;
        throw std::runtime_error("Could not open file");
    }
    std::string line;
    while (std::getline(infile, line)) lines.push_back(line);
    if(lines.empty())
    {
        std::cerr << "Error: Input file is empty: " << filename << std::endl;
        throw std::runtime_error("Empty file");
    }
    return lines;
}

static void rebuildLookAhead(const std::string& filename, std::string output_filename="", const std::string& temp_filename="temp.umf")
{
    // This function can be used to rebuild the look-ahead sizes in an existing UMF file by reading through the file and recalculating the sizes, then writing them back to the appropriate positions in the file. This can be useful if there was an error during writing that caused the look-ahead sizes to be incorrect, or if you want to add look-ahead sizes to an existing UMF file that was written without them.
    UMFReader reader(filename);
    UMFWriter writer(temp_filename);
    if(output_filename.empty()) output_filename=filename; // If no output filename is provided, overwrite the original file
    while(!reader.hasEnded())
    {
        try
        {
            UMDMolecule mol = reader.readMolecule();
            writer.writeMolecule(mol);
        }
        catch(const NoHeaderRemainingException& e)
        {
            break;
        }
    }
    writer.flush();

    // Replace the original file with the new file (if overwriting only) that has the rebuilt look-ahead sizes
    if(output_filename==filename)
    {
        if (remove(filename.c_str()) != 0)
        {
            std::cerr << "Error deleting original file: " << filename << std::endl;
            throw std::runtime_error("Could not delete original file");
        }
    }
    if (rename(temp_filename.c_str(), output_filename.c_str()) != 0)
    {
        std::cerr << "Error renaming temp file: " << temp_filename <<" to " << output_filename << std::endl;
        throw std::runtime_error("Could not rename temp file");
    }
}

static void buildUMFPointerFile(const std::string& umf_filename, std::string pointer_filename="")
{
    if(pointer_filename.empty()) pointer_filename=umf_filename+"p"; // If no pointer filename is provided, use the UMF filename with .ptr extension
    UMFReader reader(umf_filename);
    std::vector<std::string> molecule_names;
    std::vector<file_pointer> molecule_positions;
    while(!reader.hasEnded())
    {
        try
        {
            molecule_names.push_back(reader.getNextMoleculeName());
            molecule_positions.push_back(reader.getCurrentMoleculePosition());
            reader.skipMolecule();
        }
        catch(const NoHeaderRemainingException& e)
        {
            break;
        }
    }

    // Sort the positions and names in ascending order of names for efficient binary search
    std::vector<size_t> sorted_indices = argsort(molecule_names);
    std::ofstream pointer_file(pointer_filename, std::ios::binary);
    if (!pointer_file)
    {
        std::cerr << "Error opening pointer file for writing: " << pointer_filename << std::endl;
        throw std::runtime_error("Could not open pointer file");
    }
    char name_buffer[64];
    for(size_t i=0;i<sorted_indices.size();i++)
    {
        strcpy(name_buffer, molecule_names[sorted_indices[i]].c_str());
        name_buffer[63]='\0'; // Ensure null termination and limit to 63 characters
        pointer_file.write(name_buffer, sizeof(char)*64); // Write fixed-size name (64 bytes)
        pointer_file.write((char*)&molecule_positions[sorted_indices[i]], sizeof(file_pointer)); // Write corresponding file position (8 bytes on 64-bit systems)
    }
    pointer_file.flush();
    pointer_file.close();
}
#endif