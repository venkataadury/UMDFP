#ifndef RDKITINTERFACE_HPP
#define RDKITINTERFACE_HPP 1
#include <iostream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitVectUtils.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include "UMDFP.h"
#include <cstdint>

#define RDKIT_FINGERPRINT_BITS 2048

template<size_t N=RDKIT_FINGERPRINT_BITS> class CPPBitVector
{
public:
    unsigned char data[(N+7)/8]; // Allocate enough bytes to hold N bits
public:
    CPPBitVector() {memset(data, 0, sizeof(data));} // Initialize all bits to 0
    CPPBitVector(const unsigned char* bytes)
    {
        size_t length =(N+7)/8;
        memcpy(data, bytes, sizeof(data));
    }
    CPPBitVector(ExplicitBitVect* rdkit_fp)
    {
        memset(data, 0, sizeof(data)); // Initialize all bits to 0
        if(rdkit_fp->getNumBits()>N)
        {
            std::cerr << "Error: RDKit fingerprint has more bits than the CPPBitVector can hold" << std::endl;
            throw std::runtime_error("Fingerprint size exceeds CPPBitVector capacity");
        }
        for(size_t i=0;i<rdkit_fp->getNumBits();i++)
        {
            if(rdkit_fp->getBit(i)) data[i/8] |= (1 << (i%8)); // Set the corresponding bit in the data array
        }
    }

    int getNumBits() const {return N;}
    bool getBit(size_t i) const {
        if(i>=N)
        {
            std::cerr << "Error: Bit index out of range" << std::endl;
            throw std::out_of_range("Bit index exceeds CPPBitVector capacity");
        }
        return (data[i/8] >> (i%8)) & 1; // Return the value of the corresponding bit
    }

    int operator&(const CPPBitVector& other) const
    {
        // Number of common on bits in both vectors
        int count=0;
        for(size_t i=0;i<(N+7)/8;i++)
        {
            count+=__builtin_popcount(data[i] & other.data[i]);
        }
        return count;
    }

    int operator|(const CPPBitVector& other) const
    {
        // Number of on bits in either vector
        int count=0;
        for(size_t i=0;i<(N+7)/8;i++)
        {
            count+=__builtin_popcount(data[i] | other.data[i]);
        }
        return count;
    }

    int operator^(const CPPBitVector& other) const
    {
        // Number of on bits in one vector but not the other
        int count=0;
        for(size_t i=0;i<(N+7)/8;i++)
        {
            count+=__builtin_popcount(data[i] ^ other.data[i]);
        }
        return count;
    }

    float jaccardSimilarity(const CPPBitVector& other) const
    {
        int intersection = (*this) & other;
        int union_count = (*this) | other;
        if (union_count == 0) return 1.0; // Both vectors are empty, define similarity as 1
        return static_cast<float>(intersection) / union_count;
    }

    inline ByteString toByteString() const {return ByteString(reinterpret_cast<const char*>(data), sizeof(data));}
};

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

    ExplicitBitVect* getRDKFingerprint() const
    {
        if (mol1) {
            return RDKit::RDKFingerprintMol(*mol1, /*minPath=*/1, /*maxPath=*/7, /*fpSize=*/RDKIT_FINGERPRINT_BITS, /*nBitsPerHash=*/2, /*useHs=*/true);
        } else {
            std::cerr << "Error: Molecule is not initialized." << std::endl;
            return nullptr;
        }
    }
    ExplicitBitVect* getMorganFingerprint(int radius=2) const
    {
        std::unique_ptr<RDKit::FingerprintGenerator<std::uint32_t>> generator(RDKit::MorganFingerprint::getMorganGenerator<std::uint32_t>(radius, RDKIT_FINGERPRINT_BITS));
        if (mol1) {
            return generator->getFingerprint(*mol1);
        } else {
            std::cerr << "Error: Molecule is not initialized." << std::endl;
            return nullptr;
        }
    }
};

static void buildUMFFingerprintFile(const std::string& umf_filename, std::string fingerprint_filename="", const std::string& fingerprint_type="morgan")
{
    if(fingerprint_filename.empty()) fingerprint_filename=umf_filename+"b"; // If no fingerprint filename is provided, use the UMF filename with .bit extension
    std::ofstream fingerprint_file(fingerprint_filename, std::ios::binary);
    if (!fingerprint_file)
    {
        std::cerr << "Error opening fingerprint file for writing: " << fingerprint_filename << std::endl;
        throw std::runtime_error("Could not open fingerprint file");
    }

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
    for(size_t i=0;i<sorted_indices.size();i++)
    {
        size_t idx = sorted_indices[i];
        reader.jumpToMoleculeAtPosition(molecule_positions[idx]);
        RDKitMolecule rdkit_mol(reader.getNextMoleculeSMILES());
        
        ExplicitBitVect* fp = (fingerprint_type == "morgan") ? rdkit_mol.getMorganFingerprint() : rdkit_mol.getRDKFingerprint();
        if(fp)
        {
            CPPBitVector cpp_fp(fp);
            ByteString byte_string = cpp_fp.toByteString();
            byte_string.dumpToStream(fingerprint_file);
        }
        else
        {
            CPPBitVector cpp_fp; // Create an empty fingerprint (all bits 0) if there was an error computing the fingerprint for this molecule
            ByteString byte_string = cpp_fp.toByteString();
            byte_string.dumpToStream(fingerprint_file);
        }
        delete fp; // Clean up the RDKit fingerprint to free memory
    }
}

template<int N> static void _blockSimilarity(const CPPBitVector<N>* batch, const CPPBitVector<N>& query_fp, float threshold, std::vector<std::pair<unsigned long long, float>>& ret, int batch_filling, const std::string& similarity_metric="tanimoto", unsigned long long total_index=0)
{
    //#pragma omp parallel for
    for(size_t i=0;i<batch_filling;i++)
    {
        float similarity;
        if(similarity_metric=="tanimoto")
            similarity = query_fp.jaccardSimilarity(batch[i]);
        else
        {
            std::cerr << "Error: Unsupported similarity metric: " << similarity_metric << std::endl;
            throw std::runtime_error("Unsupported similarity metric");
        }
        if(similarity>=threshold) ret.push_back({total_index + i, similarity});
    }
}

template<int N> static std::vector<std::pair<unsigned long long, float>> similaritySearch(const std::string& fingerprint_filename, const CPPBitVector<N>& query_fp, float threshold, size_t top_k=-1, int batch_size=4096, const std::string& similarity_metric="tanimoto")
{
    FILE* fp_file = fopen(fingerprint_filename.c_str(), "rb");
    if (!fp_file)
    {
        std::cerr << "Error opening fingerprint file for reading: " << fingerprint_filename << std::endl;
        throw std::runtime_error("Could not open fingerprint file");
    }
    std::vector<std::pair<unsigned long long, float>> results;
    CPPBitVector<N>* batch = new CPPBitVector<N>[batch_size];
    constexpr size_t fp_size_bytes = (N+7)/8;
    unsigned long long total_index=0;
    while (true)
    {
        for(int i=0; i<batch_size; i++)
        {
            size_t read_bytes = fread(batch[i].data, sizeof(unsigned char), fp_size_bytes, fp_file);
            if (read_bytes != fp_size_bytes)
            {
                if (feof(fp_file)) // End of file reached, process the remaining batch and break
                {
                    _blockSimilarity<N>(batch, query_fp, threshold, results, i, similarity_metric, total_index);
                    break;
                }
                else
                {
                    std::cerr << "Error reading fingerprint file: " << fingerprint_filename << std::endl;
                    delete[] batch;
                    fclose(fp_file);
                    throw std::runtime_error("Error reading fingerprint file");
                }
            }
        }
        _blockSimilarity<N>(batch, query_fp, threshold, results, batch_size, similarity_metric, total_index);
        if (feof(fp_file)) break; // End of file reached, break the loop
        total_index += batch_size;
    }
    delete[] batch;
    fclose(fp_file);
    // Sort results by similarity in descending order
    std::sort(results.begin(), results.end(),
        [](const std::pair<unsigned long long, float>& a, const std::pair<unsigned long long, float>& b) {return a.second > b.second;} // Sort in descending order of similarity
    );
    
    if (top_k != -1 && results.size() > top_k) results.resize(top_k); // Keep only the top_k results

    return results;
}
#endif
