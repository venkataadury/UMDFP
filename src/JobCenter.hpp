#ifndef INCLUDED_JOBCENTER_HPP
#define INCLUDED_JOBCENTER_HPP 1
#include "UMDFP.h"
#include "UMFHelpers.h"
#include "UMDBabel.h"
#include <filesystem>

#define STANDARD_SUCCESS 0 // Exit code for successful execution
#define STANDARD_FAILURE 1 // Exit code for premature stopping, but not technically a failure (e.g., no more molecules to read from the file)
#define USER_INPUT_ERROR 2 // Exit code for invalid user input (e.g., unsupported format)
#define USER_EXIT 4 // Exit code for usage (input syntax) error
#define FILE_NOT_FOUND 8 // Exit code for file not found error
#define IO_ERROR 16 // Exit code for general I/O error
#define PROCESSING_ERROR 32 // Exit code for errors during processing (e.g., error inverting matrix for alignment)
#define LOGIC_ERROR 64
#define USAGE_ERROR 128 // Exit code for usage (input syntax) error

class JobOutput;
class JSONEntry
{
    std::vector<std::string> value;
    std::vector<JSONEntry*> nested_objects;
public:
    JSONEntry(std::string value) : value({value}) {}
    JSONEntry(std::vector<std::string> value) : value(value) {}
    JSONEntry(std::vector<JSONEntry*> nested_objects) : nested_objects(nested_objects) {}
    
    std::string toString() const
    {
        if(nested_objects.size())
        {
            std::string result = "{\n";
            if(nested_objects.size()>1) result = "[";
            for(size_t i=0;i<nested_objects.size();i++)
            {
                result += nested_objects[i]->toString();
                if(i<nested_objects.size()-1) result += ",";
                result += "\n";
            }
            if(nested_objects.size()>1) result += "]";
            else result += "}";
            return result;
        }

        if(value.size()==1) return "\"" + value[0] + "\"";
        else
        {
            std::string result = "[";
            for(size_t i=0;i<value.size();i++)
            {
                result += "\"" + value[i] + "\"";
                if(i<value.size()-1) result += ",";
            }
            result += "]";
            return result;
        }
    }

    inline void addNestedObject(JSONEntry* obj) {nested_objects.push_back(obj);}
    inline void addValue(std::string val) {value.push_back(val);}

    friend class JobOutput;
};
class JobOutput
{
    protected:
        int exit_code=0;
        std::map<std::string, JSONEntry*> key_value_pairs; // For storing any additional data that may be needed in the future
    public:
        JobOutput() {}
        virtual ~JobOutput() {for(auto& kv : key_value_pairs) delete kv.second;}

        virtual void toJSON(std::ostream& out) const
        {
            out << "{\n";
            out << "  \"exit_code\": " << exit_code << ",\n";
            int count=0;
            for(const auto& kv : key_value_pairs)
            {
                out << "  \"" << kv.first << "\": " << kv.second->toString() << "";
                if(count<key_value_pairs.size()-1) out << ",";
                out << "\n";
                count++;
            }
            out << "}\n";
        }

        virtual inline int getExitCode() const =0;
        inline void store(const std::string& key, const std::string& value)
        {
            if(this->contains(key)) key_value_pairs[key]->addValue(value);
            else key_value_pairs[key]=new JSONEntry(value);
        }
        inline void store(const std::string& key, const std::vector<std::string>& value)
        {
            if(this->contains(key))
            {
                if(key_value_pairs[key]->nested_objects.size()) key_value_pairs[key]->addNestedObject(new JSONEntry(value));
                else
                {
                    std::vector<JSONEntry*> nested_objs={key_value_pairs[key],new JSONEntry(value)};
                    key_value_pairs[key]=new JSONEntry(nested_objs);
                }
            }
            else key_value_pairs[key]=new JSONEntry(value);
        }
        inline void store(const std::string& key, const std::vector<JSONEntry*>& value)
        {
            if(this->contains(key))
            {
                throw std::logic_error("Cannot store a nested object in a key that already contains a non-nested value");
            }
            key_value_pairs[key]=new JSONEntry(value);
        }
        inline void store(const std::string& key, int value) {this->store(key, std::to_string(value));}
        inline JSONEntry* get(const std::string& key) const {return key_value_pairs.at(key);}
        inline std::string getString(const std::string& key) const {return key_value_pairs.at(key)->toString();}
        inline bool contains(const std::string& key) const {return key_value_pairs.find(key)!=key_value_pairs.end();}
};

static std::ostream& operator<<(std::ostream& out, const JobOutput& jo)
{
    jo.toJSON(out);
    return out;
}

template<class T=JobOutput> class Job
{
    protected:
        const std::string name;
        std::ostream* output_stream=&std::cout;
        std::ostream* error_stream=&std::cerr;
    public:
        Job(std::string name) : name(name)
        {
            // Ensure that T is a subclass of JobOutput
            static_assert(std::is_base_of<JobOutput, T>::value, "Template parameter T must be a subclass of JobOutput");
        }
        virtual ~Job() {}
    protected:
        virtual T* execute(int argc, char** argv) = 0;
        int processGlobalArguments(int& argc, char**& argv)
        {
            for(int i=1;i<argc;i++)
            {
                std::string arg = argv[i];
                bool match=false;
                if(arg=="-silent") {this->makeSilent(); match=true;}
                else if(arg=="-hide_errors") {this->hideErrors(); match=true;}
                else if(arg=="-help" || arg=="--help" || arg=="-h") {match=true; std::cerr << this->getUsageString() << "\n"; return USER_EXIT;}
                if(match)
                {
                    // Remove the argument from argv
                    for(int j=i;j<argc-1;j++) argv[j]=argv[j+1];
                    argc--;
                    i--; // Stay at the same index to check for multiple global arguments in a row
                }
            }
            return 0;
        }
        virtual inline const std::string& getUsageString() const =0;
    public:
        inline T* run(int argc, char** argv)
        {
            int follow =this->processGlobalArguments(argc, argv);
            if(follow) return new T(follow);
            
            T* ret=this->execute(argc, argv);
            if(ret->getExitCode()==USAGE_ERROR) *error_stream << this->getUsageString() << "\n";
            return ret;
        }
        inline const std::string& getName() const {return name;}

        inline void setOutputStream(std::ostream& out) {output_stream=&out;}
        inline void setErrorStream(std::ostream& err) {error_stream=&err;}
        inline void makeSilent()
        {
            std::ofstream* devNull=new std::ofstream("/dev/null");
            output_stream=devNull;
        }
        inline void hideErrors()
        {
            std::ofstream* devNull=new std::ofstream("/dev/null");
            error_stream=devNull;
        }
};

class UMFJobOutput : public JobOutput
{
    public:
        UMFJobOutput() : UMFJobOutput(0) {}
        UMFJobOutput(int exit_status) {exit_code=exit_status;}

        inline int getExitCode() const {return exit_code;}
};

class UMFDumpJob : public Job<UMFJobOutput>
{
    std::ofstream outfile;
    public:
        UMFDumpJob() : Job("Dump UMF") {}
        UMFDumpJob(std::string name) : Job(name) {}
    protected:
        inline const std::string& getUsageString() const override
        {
            static const std::string usage = "Usage: " + this->getName() + " <UMF file> [<Output file>] [-f Format (mol2/sdf/pdbqt/smi/umd)] [-s (separate files?)] [-b (blocks?)] [-k (skip mols)] [-l (write limit)] [-noh (no hydrogens?)]\n(Note that output file must appear before any optional arguments)\n";
            return usage;
        }
        UMFJobOutput* execute(int argc, char** argv) override
        {
            if(argc<2) return new UMFJobOutput(USAGE_ERROR);
           
            std::string umf_file = argv[1];
            std::string dump_file="output.mol2"; 
            if(argc>2) dump_file=argv[2];
            std::string output_fmt="mol2";
            
            bool split=false;
            long block_size=-1, num_entries=-1, write_limit=-1;
            long skip_mols=0;
            bool write_H=true;
            for(int i=3;i<argc;i++)
            {
                std::string arg = argv[i];
                if(arg=="-f" && i+1<argc) output_fmt = argv[++i];
                else if(arg=="-s") split=true;
                else if(arg=="-b" && i+1<argc) block_size=std::stol(argv[++i]);
                else if(arg=="-k" && i+1<argc) skip_mols=std::stol(argv[++i]);
                else if(arg=="-l" && i+1<argc) write_limit=std::stol(argv[++i]);
                else if(arg=="-noh") write_H=false;
                else
                {   
                    *error_stream << "Unknown argument: " << arg << "\n";
                    //*error_stream << "Output file is needed before more arguments\n";
                    return new UMFJobOutput(USAGE_ERROR);
                }
            }

            GenericMoleculeFileFormat* formatter;
            if(output_fmt=="mol2") formatter = new Mol2Format(write_H);
            else if(output_fmt=="sdf") formatter = new SDFFormat(write_H);
            else if(output_fmt=="pdbqt") formatter = new PDBQTFormat(write_H);
            else if(output_fmt=="smi") formatter=new SMILESFormat(); // SMILES format can never have explicit hydrogens, so we ignore the write_H flag if the user specified smi as the output format
            else if(output_fmt=="umd")
            {
                formatter=new UMDFormat(); // UMD Format must always have hydrogens
                if(!write_H) *error_stream << "Warning: UMD format must always include hydrogens. The -noh flag is ignored.\n";
            }
            else
            {
                *error_stream << "Unsupported output format: " << output_fmt << ". Supported formats are: mol2, sdf, pdbqt, smi, umd\n";
                return new UMFJobOutput(USER_INPUT_ERROR);
            }

            // Check if the file exists
            std::ifstream infile(umf_file);
            if (!infile.good())
            {
                *error_stream << "Error: UMF file not found: " << umf_file << "\n";
                return new UMFJobOutput(FILE_NOT_FOUND);
            }
            else infile.close();
            
            *output_stream << "Chosen format: "<<output_fmt<<"\n";
            bool subdir=false;
            if(split)
            {
                if(block_size>0) subdir=true; 
                else *output_stream << "Will write each molecule in a new file\n";
            }
            if(block_size>0) num_entries=countMoleculesInUMF(umf_file);
            long mols_per_block=(block_size>0 && num_entries>0)?(num_entries/block_size + (num_entries%block_size==0?0:1)):(-1);
            if(mols_per_block<0 && block_size>0)
            {
                *error_stream << "Invalid block size or unable to count number of entries in UMF file. Please check the -b flag and ensure that the UMF file is valid and not corrupted." << "\n";
                return new UMFJobOutput(LOGIC_ERROR);
            }
            std::string prefix = dump_file.substr(0, dump_file.find_last_of('.'));
            std::string old_prefix=prefix;
            std::string open_file;
            if(!split && block_size<=0)
            {
                outfile.open(dump_file);
                open_file=dump_file;
            }
            else if(split && block_size>0)
            {
                std::string subdir_name = prefix + "_blocks_0";
                if(!std::filesystem::exists(subdir_name))
                {
                    if(!std::filesystem::create_directory(subdir_name))
                    {
                        *error_stream << "Failed to create subdirectory for blocks: " << subdir_name << "\n";
                        return new UMFJobOutput(IO_ERROR);
                    }
                }
                old_prefix= prefix;
                prefix=subdir_name+"/"+(prefix.substr(prefix.find_last_of("/\\")+1));
                *output_stream << "\tPrefix: "<< prefix << "\n";
                outfile.open(prefix+"_0."+output_fmt);
                open_file=prefix+"_0."+output_fmt;
            }
            else
            {
                outfile.open(prefix+"_0."+output_fmt);
                open_file=prefix+"_0."+output_fmt;
            }
            if(!outfile.good())
            {
                *error_stream << "Failed to open output file: " << open_file <<"\n";
                return new UMFJobOutput(IO_ERROR);
            }

            *output_stream << "Testing the UMF file\n";
            UMFReader reader(umf_file);
            if(!reader.hasPointerFile())
            {
                try {buildUMFPointerFile(umf_file); if(!reader.refreshPointerFile()) *error_stream << "Pointer file still not found after building\n"; return new UMFJobOutput(IO_ERROR);}
                catch(const std::exception& e) { *error_stream << "WARN: Couldn't building pointer file (do you have write permissions to " << umf_file << "p?). This file can only be used for sequential access." << std::endl; return new UMFJobOutput(IO_ERROR); }
            }
            else *output_stream << "Pointer file found as '"<<umf_file<<"p'\n";

            reader.dumpHeader();
            *output_stream << "Writing all molecules ... "; output_stream->flush();
            long long total_mols=0; // Start at 1 since we already read the first molecule
            long long total_skip=0;
            long long next_target=mols_per_block;
            int block_id=0;
            while(!reader.hasEnded())
            {
                if(skip_mols>0)
                {
                    skip_mols--;
                    try {reader.skipMolecule();}
                    catch(const MoleculeDataEndedException& e) {break;}
                    continue;
                }
                if(write_limit>0 && total_mols>=write_limit) break;
                if(split) outfile.close();
                if(block_size>0 && total_mols>=next_target)
                {
                    block_id++;
                    next_target+=mols_per_block;
                    outfile.close();
                    if(split)
                    {
                        std::string subdir_name = old_prefix + "_blocks_" + std::to_string(block_id);
                        if(!std::filesystem::exists(subdir_name))
                        {
                            if(!std::filesystem::create_directory(subdir_name))
                            {
                                *error_stream << "Failed to create subdirectory for blocks: " << subdir_name << "\n";
                                return new UMFJobOutput(1);
                            }
                        }
                        prefix=subdir_name+"/"+(old_prefix.substr(old_prefix.find_last_of("/\\")+1));
                    }
                    else
                    {
                        outfile.open(prefix+"_"+std::to_string(block_id)+"."+output_fmt);
                        if(!outfile.good())
                        {
                            *error_stream << "Failed to open output file for block " << block_id << ": " << prefix+"_"+std::to_string(block_id)+"."+output_fmt<<"\n";
                            return new UMFJobOutput(1);
                        }
                    }
                }
                auto next_mol = reader.readMolecule();
                bool flag=true;
                for(int i=0;i<next_mol.getNumBonds();i++)
                {
                    UMDBond bond_data = next_mol.getBond(i);
                    if(bond_data.getBondType()>=8) {*output_stream << "Skipping bad molecule: "<< next_mol.getName() << "\n"; flag=false; break; total_skip++;} //throw std::runtime_error("Bad molecule: "+next_mol.getName());
                }
                if(!flag) continue;
                if(split) outfile.open(prefix+"_"+std::to_string(total_mols)+"."+output_fmt);
                formatter->formatMolecule(next_mol, outfile, "GASTEIGER");
                total_mols++;
            }
            outfile.close();
            *output_stream << "done\n";
            *output_stream << "Total molecules read: " << total_mols << std::endl;
            *output_stream << "Total molecules skipped: "<<total_skip << std::endl;
            return new UMFJobOutput(0);
        }
};

class UMFCheckOutput : public JobOutput
{
    public:
        UMFCheckOutput() : UMFCheckOutput(0) {}
        UMFCheckOutput(int exit_status)  {exit_code=exit_status;}
        void toJSON(std::ostream& out) const override
        {
            out << "{\n";
            out << "  \"exit_code\": " << exit_code << "\n";
            out << "}\n";
        }

        inline int getExitCode() const {return exit_code;}
};

class UMFCheckJob : public Job<UMFCheckOutput>
{
public:
    UMFCheckJob() : Job("Check UMF") {}
    UMFCheckJob(std::string name) : Job(name) {}
protected:
    inline const std::string& getUsageString() const override
    {
        static const std::string usage = "Usage: " + this->getName() + " <UMF file>\n";
        return usage;
    }

    UMFCheckOutput* execute(int argc, char** argv) override
    {
        if(argc<2) return new UMFCheckOutput(USAGE_ERROR);
        
        std::string umf_file = argv[1];
        // Check if the file exists
        std::ifstream infile(umf_file);
        if (!infile.good())
        {
            *error_stream << "Error: UMF file not found: " << umf_file << "\n";
            return new UMFCheckOutput(FILE_NOT_FOUND);
        }
        else infile.close();

        *output_stream << "Testing the UMF file\n";
        UMFReader reader(umf_file);
        if(!reader.hasPointerFile())
        {
            try {buildUMFPointerFile(umf_file); if(!reader.refreshPointerFile()) *error_stream << "Pointer file still not found after building\n"; return new UMFCheckOutput(IO_ERROR);}
            catch(const std::exception& e) { *error_stream << "WARN: Couldn't building pointer file (do you have write permissions to " << umf_file << "p?). This file can only be used for sequential access." << std::endl; }
        }
        else *output_stream << "Pointer file found as '"<<umf_file<<"p'\n";

        reader.dumpHeader(*output_stream);
        auto first_mol = reader.readMolecule();
        *output_stream << "First molecule name: " << first_mol.getName() << std::endl;
        first_mol.summary(*output_stream);
        *output_stream << "\n\n";
        *output_stream << "NOTE: UMFChecker only checks for file corruption/incomplete writes.\nIt does not check for chemical validity of the molecules (e.g., valence errors, unrealistic bond lengths, etc.), so even if the file passes this check, it may still contain chemically invalid molecules. For a more thorough check of chemical validity, consider using a cheminformatics toolkit like RDKit or Open Babel to read through the UMF file and validate each molecule.\n";

        *output_stream << "Testing all molecules ... "; output_stream->flush();
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
        *output_stream << "done\n";
        *output_stream << "Total molecules read: " << total_mols << std::endl;
        *output_stream << "Status: OK\n";

        return new UMFCheckOutput(0);
    }
};

class UMFQueryJobOutput : public JobOutput
{
    public:
        UMFQueryJobOutput() : UMFQueryJobOutput(0) {}
        UMFQueryJobOutput(int exit_status)  {exit_code=exit_status;}

        inline int getExitCode() const {return exit_code;}
};

class UMFQueryJob : public Job<UMFQueryJobOutput>
{
protected:
    std::ostream* out_stream = nullptr;
    long file_counter=0; // Global file counter for generating unique filenames when printing to separate files
    std::ofstream current_out_file; // Global ofstream object for the current output file when printing to separate files

    int updateWritePointer(bool separate_files, const std::string& output_file, std::string extension)
    {
        if(output_file.empty()) out_stream = output_stream; // If no output file is specified, write to the Job's output stream (which may have been set to /dev/null if -silent was used)
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
                    *error_stream << "Error opening output file: " << filename << "\n";
                    return IO_ERROR;
                }
                out_stream = &current_out_file;
            }
            else out_stream = &current_out_file;
        }
        return 0;
    }

public:
    UMFQueryJob() : Job("Query UMF") {}
    UMFQueryJob(std::string name) : Job(name) {}
protected:
    inline const std::string& getUsageString() const override
    {
        static const std::string usage = "Usage: " + this->getName() + " <UMF file> <molecule name or index> [-f output_format: mol2 (default), sdf, pdbqt, smi, or umd] [-o output_file] [-s] [-noh (no hydrogens)] \n(Note: molecule index should be a non-negative integer, and molecule name should match exactly with the name in the UMF file.\nIf both a valid index and a valid name are provided, the index will take precedence.)\n";
        return usage;
    }

    UMFQueryJobOutput* execute(int argc, char** argv) override
    {
        out_stream=output_stream; // Default to the Job's output stream (which may have been set to /dev/null if -silent was used)
        GenericMoleculeFileFormat* formatter;
        if(argc<3) return new UMFQueryJobOutput(USAGE_ERROR);

        std::string UMF_file = argv[1];
        std::string query = argv[2];
        std::string output_fmt = "mol2";
        std::string output_file="";
        bool separate_files=false;
        bool write_H=true;
        for(int i=3;i<argc;i++)
        {
            std::string arg = argv[i];
            if(arg=="-f" && i+1<argc) output_fmt = argv[++i];
            else if(arg=="-o" && i+1<argc) output_file = argv[++i];
            else if(arg=="-s") separate_files=true;
            else if(arg=="-noh") write_H=false;
            else
            {   
                *error_stream << "Unknown argument: " << arg << "\n";
                return new UMFQueryJobOutput(USAGE_ERROR);
            }
        }

        if(output_fmt=="mol2") formatter = new Mol2Format(write_H);
        else if(output_fmt=="sdf") formatter = new SDFFormat(write_H);
        else if(output_fmt=="pdbqt") formatter = new PDBQTFormat(write_H);
        else if(output_fmt=="smi") formatter=new SMILESFormat(); // SMILES format can never have explicit hydrogens, so we ignore the write_H flag if the user specified smi as the output format
        else if(output_fmt=="umd")
        {
            formatter=new UMDFormat(); // UMD Format must always have hydrogens
            if(!write_H) *error_stream << "Warning: UMD format must always include hydrogens. The -noh flag is ignored.\n";
        }
        else
        {
            *error_stream << "Unsupported output format: " << output_fmt << ". Supported formats are: mol2, sdf, pdbqt, smi, umd\n";
            return new UMFQueryJobOutput(USER_INPUT_ERROR);
        }

        if(separate_files && output_file.empty())    {
            *error_stream << "Error: Output file prefix must be provided with '-o' when using '-s' flag for separate files\n";
            return new UMFQueryJobOutput(USER_INPUT_ERROR);
        }
        else if(!output_file.empty() && !separate_files)
        {
            current_out_file.open(output_file, std::ios::binary);
            if(!current_out_file)
            {
                *error_stream << "Error opening output file: " << output_file << "\n";
                return new UMFQueryJobOutput(IO_ERROR);
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
        *error_stream << "Reading UMF file: " << UMF_file << std::endl; // Keep the status messages separate from the formatted molecule output (which may be piped to another program) by writing them to stderr instead of stdout

        UMFReader reader(UMF_file);
        if(!reader.hasPointerFile())
        {
            try {buildUMFPointerFile(UMF_file); if(!reader.refreshPointerFile()) throw std::runtime_error("Pointer file still not found after building");}
            catch(const std::exception& e)
            {
                *error_stream << "Error building pointer file. Resorting to sequential search." << std::endl;
                int index=0;
                while(!reader.hasEnded())
                {
                    if(ligindex)
                    {
                        if(index==ligid)
                        {
                            int wp_exitstat = updateWritePointer(false, output_file, output_fmt); // For single molecule output, we can just use the provided output file directly without generating unique filenames for each molecule
                            if(wp_exitstat) return new UMFQueryJobOutput(wp_exitstat);
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
                                int wp_exitstat = updateWritePointer(separate_files, output_file, output_fmt);
                                if(wp_exitstat) return new UMFQueryJobOutput(wp_exitstat);
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
                *error_stream << "Didn't find molecule with index: " << ligid << std::endl;
                return new UMFQueryJobOutput(STANDARD_FAILURE);
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
                    *error_stream << "No pointer file found, cannot perform name query. Please provide a pointer file for quick access by name.\n";
                    return new UMFQueryJobOutput(FILE_NOT_FOUND);
                }
                catch(const NoSuchMoleculeException& e)
                {
                    *output_stream << "Didn't find molecule with name: " << q << std::endl;
                    continue;
                }

                reader.jumpToMoleculeAtPosition(position);
                try
                {
                    UMDMolecule mol = reader.readMolecule();
                    int wp_exitstat = updateWritePointer(separate_files, output_file, output_fmt);
                    if(wp_exitstat) return new UMFQueryJobOutput(wp_exitstat);
                    formatter->formatMolecule(mol, *out_stream, "GASTEIGER");
                }
                catch(const NoHeaderRemainingException& e)
                {
                    *output_stream << "Didn't find molecule with name: " << q << std::endl;
                    continue;
                }
            }
        }
        return new UMFQueryJobOutput(0);
    }
};

class UMDBasedConvJobOutput : public JobOutput
{
    public:
        UMDBasedConvJobOutput() : UMDBasedConvJobOutput(0) {}
        UMDBasedConvJobOutput(int exit_status)  {exit_code=exit_status;}

        inline int getExitCode() const {return exit_code;}
};

class UMDBasedConvJob : public Job<UMDBasedConvJobOutput>
{
protected:
    inline int setup_formatter(GenericMoleculeFileFormat*& formatter, const std::string& fmt)
    {
        if(fmt=="mol2") formatter = new Mol2Format();
        else if(fmt=="sdf") formatter = new SDFFormat();
        else if(fmt=="pdbqt") formatter = new PDBQTFormat();
        else if(fmt=="smi") formatter=new SMILESFormat();
        else if(fmt=="umd") formatter=new UMDFormat();
        else
        {
            *error_stream << "Unsupported format: " << fmt << ". Supported formats are: mol2, sdf, pdbqt, smi, and umd\n";
            return USER_INPUT_ERROR;
        }
        return 0;
    }
public:
    UMDBasedConvJob() : Job("Convert to/from UMF") {}
    UMDBasedConvJob(std::string name) : Job(name) {}
protected:
    inline const std::string& getUsageString() const override
    {
        static const std::string usage = "Usage: " + this->getName() + " <input file> [<template UMD file> (Default: input file)] [-if <input format> (Default: umd)] [-of <output format> (Default: mol2)] [-o <output file>] [-e (Read names from stdin. Default: false)] [-noh (Do not expect hydrogens)]\n(Note: If the -e flag is used to read molecule names from stdin, the program will read each line from stdin as a molecule name and attempt\nto find a matching molecule in the template UMD file (or input file if no template is provided) to use as a template for the conversion. The converted molecule will then\nbe written to the output file or stdout. If no matching molecule is found for a given name, that name will be skipped with a warning message printed to stderr.)\nSupported formats are: umd, mol2, sdf, pdbqt, smi\n";
        return usage;
    }

    UMDBasedConvJobOutput* execute(int argc, char** argv) override
    {
        GenericMoleculeFileFormat* input_formatter;
        GenericMoleculeFileFormat* output_formatter;
        if(argc<2) return new UMDBasedConvJobOutput(USAGE_ERROR);
        
        std::string input_file = argv[1];
        std::string template_file = input_file;
        std::string input_fmt = "umd";
        std::string output_fmt = "mol2";
        std::string output_file = "";
        bool read_names_from_stdin=false;
        bool template_chosen=false;
        bool has_H=true;
        for(int i=2;i<argc;i++)
        {
            std::string arg = argv[i];
            if(arg=="-if" && i+1<argc) input_fmt = argv[++i];
            else if(arg=="-of" && i+1<argc) output_fmt = argv[++i];
            else if(arg=="-o" && i+1<argc) output_file = argv[++i];
            else if(arg=="-noh") has_H=false;
            else if(arg=="-e") read_names_from_stdin=true;
            else if(arg[0]!='-' && !template_chosen) {
                template_file = arg;
                template_chosen=true;
            }
            else
            {
                *error_stream << "Unknown argument: " << arg << "\n";
                return new UMDBasedConvJobOutput(USER_INPUT_ERROR);
            }
        }

        int fmt_setup_status=0;
        if((fmt_setup_status=setup_formatter(input_formatter, input_fmt))) return new UMDBasedConvJobOutput(fmt_setup_status);
        if((fmt_setup_status=setup_formatter(output_formatter, output_fmt))) return new UMDBasedConvJobOutput(fmt_setup_status);

        if(!template_chosen && input_fmt!="umd")
        {
            *error_stream << "FATAL: Cannot use a non-UMD template file. Please use either a UMD input file or specify a UMD template file.\n";
            return new UMDBasedConvJobOutput(USER_INPUT_ERROR);
        }

        UMDReader reader(template_file);
        std::map<std::string, UMDMolecule> template_molecules = drainUMDReader(reader);
        std::ofstream current_out_file;
        //std::ostream* out_stream;
        if(!output_file.empty())
        {
            current_out_file.open(output_file, std::ios::binary);
            if(!current_out_file)        {
                *error_stream << "Error opening output file: " << output_file << "\n";
                return new UMDBasedConvJobOutput(IO_ERROR);
            }
        }

        std::ifstream infile(input_file);

        if(read_names_from_stdin)
        {
            std::string name;
            while(std::getline(std::cin, name))
            {
                name.erase(name.find_last_not_of(" \n\r\t")+1); // Remove trailing whitespaces
                if(name.empty()) continue;
                if(template_molecules.find(name)==template_molecules.end())
                {
                    *error_stream << "Warning: Molecule with name '" << name << "' not found in template file. Skipping.\n";
                    continue;
                }
                UMDMolecule mol = template_molecules[name];
                std::string temp_input_file=name+"."+input_fmt;
                std::ifstream temp_in(temp_input_file);
                std::istream* in_stream;
                if(!temp_in)
                {
                    *error_stream << "Warning: Could not open input file for molecule '" << name << "' with expected filename '" << temp_input_file << "'. Falling back to 'input_file' given: " << input_file << "\n";
                    in_stream = &infile;
                    if(!infile.good())
                    {
                        *error_stream << "Error opening input file: " << temp_input_file << "\n";
                        return new UMDBasedConvJobOutput(IO_ERROR);
                    }
                }
                else in_stream = &temp_in;
                
                loadCoordinatesFromFormat(mol, *in_stream, *input_formatter);
                output_formatter->formatMolecule(mol, current_out_file, "GASTEIGER");
            }
        }
        else
        {
            *error_stream << "No names provided in stdin (use '-e' to do this). Processing all molecules in template file: " << template_file << " in order\n";
            for(auto& pair : template_molecules)
            {
                const std::string& name = pair.first;
                UMDMolecule& mol = pair.second;
                std::string temp_input_file=name+"."+input_fmt;
                std::ifstream temp_in(temp_input_file);
                std::istream* in_stream;
                if(!temp_in)
                {
                    in_stream = &infile;
                    if(!infile.good())
                    {
                        *error_stream << "Error opening input file: " << temp_input_file << "\n";
                        return new UMDBasedConvJobOutput(IO_ERROR);
                    }
                }
                else in_stream = &temp_in;
                
                loadCoordinatesFromFormat(mol, *in_stream, *input_formatter);
                output_formatter->formatMolecule(mol, current_out_file, "GASTEIGER");
            }
        }
        return new UMDBasedConvJobOutput(0);
    }
};
#endif