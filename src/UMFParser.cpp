#include "UMDFP.h"
#include "UMDBabel.h"

int main(int argc, char** argv)
{
    // Pick your formatter
    Mol2Format formatter;

    std::string UMF_file = argv[1];
    std::cout << "Reading UMF file: " << UMF_file << std::endl;
    
    UMFReader reader(UMF_file);
    while(true)
    {
        std::cout << "Enter query (Ligand name or index, or 'exit' to quit): ";
        std::string query;
        std::getline(std::cin, query);
        if(query=="exit") break;
        int ligid;
        bool ligindex=true;
        try{ligid = std::stoi(query);}
        catch(const std::exception& e) {ligindex=false;}
        if(ligindex) std::cout << "Querying for molecule with index: " << ligid << std::endl;
        else std::cout << "Querying for molecule with name: " << query << std::endl;

        if(ligindex)
        {
            try
            {
                if(ligid!=0) reader.skipMolecule(ligid);
                UMDMolecule mol = reader.readMolecule();
                formatter.formatMolecule(mol, std::cout);
            }
            catch(const NoHeaderRemainingException& e)
            {
                std::cout << "Didn't find molecule with index: " << ligid << std::endl;
                break;
            }
        }
        else
        {
            file_pointer position;
            try {position = reader.query_by_name(query);}
            catch(const NoPointerFileFoundException& e)
            {
                std::cout << "No pointer file found, cannot perform name query. Please provide a pointer file for quick access by name.\n";
                break;
            }
            catch(const NoSuchMoleculeException& e)
            {
                std::cout << "Didn't find molecule with name: " << query << std::endl;
                continue;
            }

            reader.jumpToMoleculeAtPosition(position);
            try
            {
                UMDMolecule mol = reader.readMolecule();
                formatter.formatMolecule(mol, std::cout, "GASTEIGER");
            }
            catch(const NoHeaderRemainingException& e)
            {
                std::cout << "Didn't find molecule with name: " << query << std::endl;
                break;
            }
        }
    }
    return 0;
}