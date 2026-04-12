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
    buildUMFPointerFile(output_file);
    std::cout << "Write is complete, now testing the reader...\n\n\n" << std::endl;

    Mol2Format formatter;

    // Test the reader
    std::string query="Z1643219758";
    UMFReader reader(output_file);
    while(!reader.hasEnded())
    {
        try
        {
            std::string mol_name = reader.getHeader().name;
            if(mol_name==query)
            {
                //for(int i=0;i<258;i++) reader.skipMolecule(); // Skip 258 molecules to test both the 256 and 65536 look-aheads
                UMDMolecule mol = reader.readMolecule();
                mol.summary();
                std::cout << "\n";
                formatter.formatMolecule(mol, std::cout);
                break; // Remove this break statement in the future to read through the entire file
            }
            else reader.skipMolecule();
        }
        catch(const NoHeaderRemainingException& e)
        {
            std::cout << "No more molecules to read, exiting." << std::endl;
            break;
        }
    }
    return 0;
}