//#define ENABLE_CPU_PARALLELISM 1
#include "UMFHelpers.h"
#include "RDKitInterface.hpp"
#include <iostream>

int main(int argc, char** argv)
{
    if(argc<3)
    {
        std::cerr << "Usage: " << argv[0] << " <UMF file> <Query Smiles File> [threshold (default: 0.7)] [top_k (default: all)]" << std::endl;
        return 1;
    }
    std::string umf_filename = argv[1];
    std::string query_smiles_file = argv[2];
    float threshold = (argc > 3) ? std::stof(argv[3]) : 0.7;
    size_t top_k = (argc > 4) ? std::stoull(argv[4]) : -1;
    
    std::ifstream qsmi(query_smiles_file);
    if(!qsmi.good())
    {
        std::cerr << "Error opening query smiles file "<<query_smiles_file<<"\n";
        return 1;
    }
    
    while(!qsmi.eof())
    {
        std::string line;
        if(!std::getline(qsmi,line)) break;
        std::stringstream lss(line);
        std::string query_smiles,name;
        lss >> query_smiles >> name;
        
        RDKitMolecule query_mol(query_smiles);
        ExplicitBitVect* query_fp = query_mol.getMorganFingerprint();
        if(!query_fp)
        {
            std::cerr << "Error computing fingerprint for query molecule" << std::endl;
            return 1;
        }
        //std::cout << "Query molecule loaded\n\n";

        try {std::ifstream test_file(umf_filename+"b", std::ios::binary); if(!test_file) throw std::runtime_error("Fingerprint file not found"); else test_file.close();}
        catch(const std::exception& e)
        {
            std::cerr << "Building the fingerprint file ... "; std::cerr.flush();
            buildUMFFingerprintFile(umf_filename);
            std::cerr << "Done\n\n";
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
        std::cout << name << " " << results.size() << "\t";
        for(const auto& [index, similarity] : results)
        {
            std::pair<std::string, file_pointer> query_result = getQueryByIndex(pfile, index);
            //std::cout << "Index: " << index << ", Name: " << query_result.first << ", Similarity: " << similarity << std::endl;
            std::cout << query_result.first << " " << similarity << " ";
        }
        std::cout << "\n";
    }
    
    return 0;
}
