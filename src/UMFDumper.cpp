#include "UMDFP.h"
#include "UMDBabel.h"

static std::ofstream outfile;
int main(int argc, char** argv)
{
    if(argc<2)
    {
        std::cerr << "Usage: " << argv[0] << " <UMF file> [<Output file>] [-f Format (mol2/sdf/pdbqt/smi)] [-s (separate files?)]\n";
        return 1;
    }
    std::string umf_file = argv[1];
    std::string dump_file="output.mol2"; 
    if(argc>2) dump_file=argv[2];
    std::string output_fmt="mol2";
    
    bool split=false;
    for(int i=3;i<argc;i++)
    {
        std::string arg = argv[i];
        if(arg=="-f" && i+1<argc) output_fmt = argv[++i];
        else if(arg=="-s") split=true;
        else
        {   
            std::cerr << "Unknown argument: " << arg << "\n";
            std::cerr << "Output file is needed before more arguments\n";
            return 1;
        }
    }

    GenericMoleculeFileFormat* formatter;
    if(output_fmt=="mol2") formatter = new Mol2Format();
    else if(output_fmt=="sdf") formatter = new SDFFormat();
    else if(output_fmt=="pdbqt") formatter = new PDBQTFormat();
    else if(output_fmt=="smi") formatter=new SMILESFormat();
    else
    {
        std::cerr << "Unsupported output format: " << output_fmt << ". Supported formats are: mol2, sdf, pdbqt\n";
        return 1;
    }

    // Check if the file exists
    std::ifstream infile(umf_file);
    if (!infile.good())
    {
        std::cerr << "Error: UMF file not found: " << umf_file << "\n";
        return 1;
    }
    else infile.close();
    
    std::cout << "Chosen format: "<<output_fmt<<"\n";
    if(split) std::cout << "Will write each molecule in a new file\n";

    std::string prefix = dump_file.substr(0, dump_file.find_last_of('.'));
    if(!split) outfile.open(dump_file);
    else outfile.open(prefix+"_0."+output_fmt);
    if(!outfile.good())
    {
        std::cerr << "Failed to open output file: " << dump_file<<"\n";
        return 1;
    }

    std::cout << "Testing the UMF file\n";
    UMFReader reader(umf_file);
    if(!reader.hasPointerFile())
    {
        try {buildUMFPointerFile(umf_file); if(!reader.refreshPointerFile()) throw std::runtime_error("Pointer file still not found after building");}
        catch(const std::exception& e) {std::cerr << "WARN: Couldn't building pointer file (do you have write permissions to " << umf_file << "p?). This file can only be used for sequential access." << std::endl;}
    }
    else std::cout << "Pointer file found as '"<<umf_file<<"p'\n";

    reader.dumpHeader();
    std::cout << "Writing all molecules ... "; std::cout.flush();
    long long total_mols=0; // Start at 1 since we already read the first molecule
    long long total_skip=0;
    while(!reader.hasEnded())
    {
        if(split) outfile.close();
        auto next_mol = reader.readMolecule();
        bool flag=true;
        for(int i=0;i<next_mol.getNumBonds();i++)
        {
            UMDBond bond_data = next_mol.getBond(i);
            if(bond_data.getBondType()>=8) {std::cout << "Skipping bad molecule: "<< next_mol.getName() << "\n"; flag=false; break; total_skip++;} //throw std::runtime_error("Bad molecule: "+next_mol.getName());
        }
        if(!flag) continue;
        if(split) outfile.open(prefix+"_"+std::to_string(total_mols)+"."+output_fmt);
        formatter->formatMolecule(next_mol, outfile, "GASTEIGER");
        total_mols++;
    }
    outfile.close();
    std::cout << "done\n";
    std::cout << "Total molecules read: " << total_mols << std::endl;
    std::cout << "Total molecules skipped: "<<total_skip << std::endl;
    return 0;
}
