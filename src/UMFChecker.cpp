#include "UMDFP.h"
#include "UMDBabel.h"

int main(int argc, char** argv)
{
    if(argc<2)
    {
        std::cerr << "Usage: " << argv[0] << " <UMF file>\n";
        return 1;
    }
    std::string umf_file = argv[1];
    // Check if the file exists
    std::ifstream infile(umf_file);
    if (!infile.good())
    {
        std::cerr << "Error: UMF file not found: " << umf_file << "\n";
        return 1;
    }
    else infile.close();

    std::cout << "Testing the UMF file\n";
    UMFReader reader(umf_file);
    if(!reader.hasPointerFile())
    {
        try {buildUMFPointerFile(umf_file); if(!reader.refreshPointerFile()) throw std::runtime_error("Pointer file still not found after building");}
        catch(const std::exception& e) {std::cerr << "WARN: Couldn't building pointer file (do you have write permissions to " << umf_file << "p?). This file can only be used for sequential access." << std::endl;}
    }
    else std::cout << "Pointer file found as '"<<umf_file<<"p'\n";

    reader.dumpHeader();
    auto first_mol = reader.readMolecule();
    std::cout << "First molecule name: " << first_mol.getName() << std::endl;
    first_mol.summary();
    std::cout << "\n\n";

    std::cout << "Testing all molecules ... "; std::cout.flush();
    long long total_mols=1; // Start at 1 since we already read the first molecule
    while(!reader.hasEnded())
    {
        auto next_mol = reader.readMolecule();
        total_mols++;
        for(int i=0;i<next_mol.getNumBonds();i++)
        {
            UMDBond bond_data = next_mol.getBond(i);
            if(bond_data.getBondType()>=8) throw std::runtime_error("Bad molecule: "+next_mol.getName());
        }
    }
    std::cout << "done\n";
    std::cout << "Total molecules read: " << total_mols << std::endl;
    std::cout << "Status: OK\n";
    return 0;
}