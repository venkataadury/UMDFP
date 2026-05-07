#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "UMDBabel.h"
#include "UMFHelpers.h"

inline void setup_formatter(GenericMoleculeFileFormat*& formatter, const std::string& fmt)
{
    if(fmt=="mol2") formatter = new Mol2Format();
    else if(fmt=="sdf") formatter = new SDFFormat();
    else if(fmt=="pdbqt") formatter = new PDBQTFormat();
    else if(fmt=="smi") formatter=new SMILESFormat();
    else if(fmt=="umd") formatter=new UMDFormat();
    else
    {
        std::cerr << "Unsupported format: " << fmt << ". Supported formats are: mol2, sdf, pdbqt, smi, and umd\n";
        throw std::runtime_error("Unsupported format");
    }
}
int main(int argc, char** argv)
{
    GenericMoleculeFileFormat* input_formatter;
    GenericMoleculeFileFormat* output_formatter;
    if(argc<2)
    {
        std::cerr << "Usage: " << argv[0] << " <input file> [<template UMD file> (Default: input file)] [-if <input format> (Default: umd)] [-of <output format> (Default: mol2)] [-o <output file>] [-e (Read names from stdin. Default: false)]\n";
        return 1;
    }
    std::string input_file = argv[1];
    std::string template_file = input_file;
    std::string input_fmt = "umd";
    std::string output_fmt = "mol2";
    std::string output_file = "";
    bool read_names_from_stdin=false;
    bool template_chosen=false;
    for(int i=2;i<argc;i++)
    {
        std::string arg = argv[i];
        if(arg=="-if" && i+1<argc) input_fmt = argv[++i];
        else if(arg=="-of" && i+1<argc) output_fmt = argv[++i];
        else if(arg=="-o" && i+1<argc) output_file = argv[++i];
        else if(arg=="-e") read_names_from_stdin=true;
        else if(arg[0]!='-' && !template_chosen) {
            template_file = arg;
            template_chosen=true;
        }
        else
        {
            std::cerr << "Unknown argument: " << arg << "\n";
            return 1;
        }
    }

    setup_formatter(input_formatter, input_fmt);
    setup_formatter(output_formatter, output_fmt);

    if(!template_chosen && input_fmt!="umd")
    {
        std::cerr << "FATAL: Cannot use a non-UMD template file. Please use either a UMD input file or specify a UMD template file.\n";
        return 1;
    }

    UMDReader reader(template_file);
    std::map<std::string, UMDMolecule> template_molecules = drainUMDReader(reader);
    std::ofstream current_out_file;
    //std::ostream* out_stream;
    if(!output_file.empty())
    {
        current_out_file.open(output_file, std::ios::binary);
        if(!current_out_file)        {
            std::cerr << "Error opening output file: " << output_file << "\n";
            return 1;
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
                std::cerr << "Warning: Molecule with name '" << name << "' not found in template file. Skipping.\n";
                continue;
            }
            UMDMolecule mol = template_molecules[name];
            std::string temp_input_file=name+"."+input_fmt;
            std::ifstream temp_in(temp_input_file);
            std::istream* in_stream;
            if(!temp_in)
            {
                std::cerr << "Warning: Could not open input file for molecule '" << name << "' with expected filename '" << temp_input_file << "'. Falling back to 'input_file' given: " << input_file << "\n";
                in_stream = &infile;
                if(!infile.good())
                {
                    std::cerr << "Error opening input file: " << temp_input_file << "\n";
                    return 1;
                }
            }
            else in_stream = &temp_in;
            
            loadCoordinatesFromFormat(mol, *in_stream, *input_formatter);
            output_formatter->formatMolecule(mol, current_out_file, "GASTEIGER");
        }
    }
    else
    {
        std::cerr << "No names provided in stdin. Processing all molecules in template file: " << template_file << " in order\n";
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
                    std::cerr << "Error opening input file: " << temp_input_file << "\n";
                    return 1;
                }
            }
            else in_stream = &temp_in;
            
            loadCoordinatesFromFormat(mol, *in_stream, *input_formatter);
            output_formatter->formatMolecule(mol, current_out_file, "GASTEIGER");
        }
    }
}