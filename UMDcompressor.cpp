#include "UMDFP.h"
#include "UMDBabel.h"

int main(int argc, char** argv)
{
    std::string output_file = "test.umf";
    if(argc>=3) output_file = argv[2];
    std::cout << "Output file: " << output_file << std::endl;
    /*std::cout << "Size of UMDAtomData: " << sizeof(UMDAtomData) << " bytes, Size of UMDBondData: " << sizeof(UMDBondData) << " bytes" << std::endl; // Debug line (will be removed in the future)
    exit(0);*/
    UMFWriter writer(output_file);

    UMDReader base_reader("GastUMD.umd");
    while(true)
    {
        try
        {
            UMDMolecule mol = base_reader.getNextMolecule();
            writer.writeMolecule(mol);
        }
        catch(const MoleculeDataEndedException& e)
        {
            std::cout << "No more molecules to read from input file, finished writing UMF file." << std::endl;
            break;
        }
    }
    writer.flush();
    std::cout << "Finished writing UMF file, now building pointer file for quick access by name ..." << std::endl;
    buildUMFPointerFile(output_file);
    std::cout << "done\n";
    std::cout << "Write is complete, now testing the reader...\n\n\n" << std::endl;


    // Test the reader
    std::cout << "Testing the UMF file\n";
    UMFReader reader(output_file);
    //reader._forceVerifyNextHeader(); // Force reading the first header to initialize the reader for testing
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