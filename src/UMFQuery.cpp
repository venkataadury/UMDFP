#include "UMDFP.h"
#include "UMDBabel.h"

static std::ostream* out_stream = &std::cout; // Global output stream pointer (used for printing status messages to stderr while allowing formatted molecule output to be printed to stdout or a file)
static long file_counter=0; // Global file counter for generating unique filenames when printing to separate files
static std::ofstream current_out_file; // Global ofstream object for the current output file when printing to separate files
static void updateWritePointer(bool separate_files, const std::string& output_file, std::string extension)
{
    if(output_file.empty()) out_stream = &std::cout;
    else
    {
        if(separate_files)
        {
            std::string prefix = output_file.substr(0, output_file.find_last_of('.')); // Get the prefix of the output file (without extension) to use for generating unique filenames for each molecule
            std::string filename = prefix + "_" + std::to_string(file_counter++) + "." + extension; // Generate a unique filename for each molecule
            if(current_out_file.is_open()) current_out_file.close();
            current_out_file.open(filename, std::ios::binary);
            if(!current_out_file)
            {
                std::cerr << "Error opening output file: " << filename << "\n";
                throw std::runtime_error("Could not open output file");
            }
            out_stream = &current_out_file;
        }
        else out_stream = &current_out_file;
    }
}

int main(int argc, char** argv)
{
    // Pick your formatter
    //Mol2Format formatter;
    GenericMoleculeFileFormat* formatter;
    if(argc<3)
    {
        std::cerr << "Usage: " << argv[0] << " <UMF file> <molecule name or index> [-f output_format: mol2 (default) or sdf] [-o output_file] [-s]\n"; //-s flag determines if each molecule is printed to a separate file (depending on the output format) instead of all being printed to the same file. Needs '-o'
        return 1;
    }

    std::string UMF_file = argv[1];
    std::string query = argv[2];
    std::string output_fmt = "mol2";
    std::string output_file="";
    bool separate_files=false;
    for(int i=3;i<argc;i++)
    {
        std::string arg = argv[i];
        if(arg=="-f" && i+1<argc) output_fmt = argv[++i];
        else if(arg=="-o" && i+1<argc) output_file = argv[++i];
        else if(arg=="-s") separate_files=true;
        else
        {   
            std::cerr << "Unknown argument: " << arg << "\n";
            return 1;
        }
    }

    if(output_fmt=="mol2") formatter = new Mol2Format();
    else if(output_fmt=="sdf") formatter = new SDFFormat();
    else if(output_fmt=="pdbqt") formatter = new PDBQTFormat();
    else
    {
        std::cerr << "Unsupported output format: " << output_fmt << ". Supported formats are: mol2, sdf, pdbqt\n";
        return 1;
    }

    if(separate_files && output_file.empty())    {
        std::cerr << "Error: Output file prefix must be provided with '-o' when using '-s' flag for separate files\n";
        return 1;
    }
    else if(!output_file.empty() && !separate_files)
    {
        current_out_file.open(output_file, std::ios::binary);
        if(!current_out_file)
        {
            std::cerr << "Error opening output file: " << output_file << "\n";
            return 1;
        }
    }


    int ligid;
    bool ligindex=true;
    std::vector<std::string> queries;
    // If query is a valid file, then load queries from the file (one per line)
    std::ifstream query_file(query);
    if(query_file.good())
    {
        std::string line;
        while(std::getline(query_file, line))
        {
            line.erase(line.find_last_not_of(" \n\r\t")+1); // Remove trailing whitespaces
            if(line.empty()) continue;
            queries.push_back(line);
        }
        query_file.close();
        ligindex=false; // If we're loading queries from a file, then they must be molecule names, not indices
    }
    else queries.push_back(query); // Otherwise, treat the query as a single molecule name or index
    
    try{ligid = std::stoi(query);}
    catch(const std::exception& e) {ligindex=false;}
    std::cerr << "Reading UMF file: " << UMF_file << std::endl; // Keep the status messages separate from the formatted molecule output (which may be piped to another program) by writing them to stderr instead of stdout
    
    UMFReader reader(UMF_file);
    if(!reader.hasPointerFile())
    {
        try {buildUMFPointerFile(UMF_file); if(!reader.refreshPointerFile()) throw std::runtime_error("Pointer file still not found after building");}
        catch(const std::exception& e)
        {
            std::cerr << "Error building pointer file. Resorting to sequential search." << std::endl;
            int index=0;
            while(!reader.hasEnded())
            {
                if(ligindex)
                {
                    if(index==ligid)
                    {
                        updateWritePointer(false, output_file, output_fmt); // For single molecule output, we can just use the provided output file directly without generating unique filenames for each molecule
                        UMDMolecule mol = reader.readMolecule();
                        formatter->formatMolecule(mol, *out_stream, "GASTEIGER");
                        return 0;
                    }
                }
                else
                {
                    std::string name = reader.getNextMoleculeName();
                    for(const auto& q : queries)
                    {
                        if(name==q)
                        {
                            updateWritePointer(separate_files, output_file, output_fmt);
                            UMDMolecule mol = reader.readMolecule();
                            formatter->formatMolecule(mol, *out_stream, "GASTEIGER");
                        }
                    }
                }
                reader.getNextMoleculeName();
                index++;
            }
        }
    }

    file_pointer position;
    if(ligindex)
    {
        try
        {
            if(ligid!=0) reader.skipMolecule(ligid);
            UMDMolecule mol = reader.readMolecule();
            updateWritePointer(false, output_file, output_fmt); // For single molecule output, we can just use the provided output file directly without generating unique filenames for each molecule
            formatter->formatMolecule(mol, *out_stream, "GASTEIGER");
        }
        catch(const NoHeaderRemainingException& e)
        {
            std::cout << "Didn't find molecule with index: " << ligid << std::endl;
            return 0;
        }
    }
    else
    {
        file_pointer position;
        for(const auto& q : queries)
        {
            try {position = reader.query_by_name(q);}
            catch(const NoPointerFileFoundException& e)
            {
                std::cout << "No pointer file found, cannot perform name query. Please provide a pointer file for quick access by name.\n";
                return 1;
            }
            catch(const NoSuchMoleculeException& e)
            {
                std::cerr << "Didn't find molecule with name: " << q << std::endl;
                continue;
            }

            reader.jumpToMoleculeAtPosition(position);
            try
            {
                UMDMolecule mol = reader.readMolecule();
                updateWritePointer(separate_files, output_file, output_fmt);
                formatter->formatMolecule(mol, *out_stream, "GASTEIGER");
            }
            catch(const NoHeaderRemainingException& e)
            {
                std::cerr << "Didn't find molecule with name: " << q << std::endl;
                continue;
            }
        }
    }
    return 0;
}