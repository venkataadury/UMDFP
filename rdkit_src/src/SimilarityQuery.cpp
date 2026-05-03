#include "RDKitInterface.hpp"
#include "UMFHelpers.h"
#include <iostream>

int main(int argc, char** argv)
{
    if(argc<3)
    {
        std::cerr << "Usage: " << argv[0] << " <UMF file> <Query Smiles> [threshold (default: 0.7)] [top_k (default: all)]" << std::endl;
        return 1;
    }
    std::string umf_filename = argv[1];
    std::string query_smiles = argv[2];
    float threshold = (argc > 3) ? std::stof(argv[3]) : 0.7;
    size_t top_k = (argc > 4) ? std::stoull(argv[4]) : -1;
    RDKitMolecule query_mol(query_smiles);
    ExplicitBitVect* query_fp = query_mol.getMorganFingerprint();
    if(!query_fp)
    {
        std::cerr << "Error computing fingerprint for query molecule" << std::endl;
        return 1;
    }
    std::cout << "Query molecule loaded\n\n";

    try {std::ifstream test_file(umf_filename+"b", std::ios::binary); if(!test_file) throw std::runtime_error("Fingerprint file not found"); else test_file.close(); std::cout << "Fingerprint file found\n\n";}
    catch(const std::exception& e)
    {
        std::cout << "Building the fingerprint file ... "; std::cout.flush();
        buildUMFFingerprintFile(umf_filename);
        std::cout << "Done\n\n";
    }

    CPPBitVector query_cpp_fp(query_fp);
    std::ifstream fingerprint_file(umf_filename+"b", std::ios::binary);
    if (!fingerprint_file)
    {
        std::cerr << "Error opening fingerprint file for reading: " << umf_filename+"b" << std::endl;
        return 1;
    }
    fingerprint_file.close();

    std::string fpfile=umf_filename+"b";
    std::string pfile=umf_filename+"p";
    std::vector<std::pair<unsigned long long, float>> results = similaritySearch<RDKIT_FINGERPRINT_BITS>(fpfile, query_cpp_fp, threshold, top_k);
    for(const auto& [index, similarity] : results)
    {
        std::pair<std::string, file_pointer> query_result = getQueryByIndex(pfile, index);
        std::cout << "Index: " << index << ", Name: " << query_result.first << ", Similarity: " << similarity << std::endl;
    }
    
    return 0;
}