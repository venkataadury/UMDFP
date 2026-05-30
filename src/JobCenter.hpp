#ifndef INCLUDED_JOBCENTER_HPP
#define INCLUDED_JOBCENTER_HPP 1
#include "UMDFP.h"
#include "UMFHelpers.h"
#include "UMDBabel.h"
#include <filesystem>

#define USER_INPUT_ERROR 1 // Exit code for invalid user input (e.g., unsupported format)
#define USER_EXIT 2 // Exit code for usage (input syntax) error
#define FILE_NOT_FOUND 3 // Exit code for file not found error
#define IO_ERROR 4 // Exit code for general I/O error
#define LOGIC_ERROR 15
#define USAGE_ERROR 222 // Exit code for usage (input syntax) error


class JobOutput
{
    protected:
        int exit_code=0;
    public:
        JobOutput() {}
        virtual ~JobOutput() {}
        virtual void toJSON(std::ostream& out) const = 0;
        std::ostream& operator<<(std::ostream& out) const
        {
            this->toJSON(out);
            return out;
        }

        virtual inline int getExitCode() const =0;
};

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
            if(((JobOutput*)ret)->getExitCode()==USAGE_ERROR) *error_stream << this->getUsageString() << "\n";
            return ret;
        }
        inline const std::string& getName() const {return name;}

        inline void setOutputStream(std::ostream& out) {output_stream=&out;}
        inline void setErrorStream(std::ostream& err) {error_stream=&err;}
        inline void makeSilent()
        {
            std::ofstream devNull("/dev/null");
            output_stream=&devNull;
        }
        inline void hideErrors()
        {
            std::ofstream devNull("/dev/null");
            error_stream=&devNull;
        }
};

class UMFJobOutput : public JobOutput
{
    protected:
        std::map<std::string, std::string> key_value_pairs;
    public:
        UMFJobOutput() : UMFJobOutput(0) {}
        UMFJobOutput(int exit_status) {exit_code=exit_status;}
        void toJSON(std::ostream& out) const override
        {
            out << "{\n";
            out << "  \"exit_code\": " << exit_code << ",\n";
            out << "  \"data\": {\n";
            int count=0;
            for(const auto& kv : key_value_pairs)
            {
                out << "    \"" << kv.first << "\": \"" << kv.second << "\"";
                if(count<key_value_pairs.size()-1) out << ",";
                out << "\n";
                count++;
            }
            out << "  }\n";
            out << "}\n";
        }

        inline int getExitCode() const {return exit_code;}
};

class UMFDumpJob : public Job<UMFJobOutput>
{
    std::ofstream outfile;
    public:
        UMFDumpJob() : Job("Dump UMF") {}
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
#endif