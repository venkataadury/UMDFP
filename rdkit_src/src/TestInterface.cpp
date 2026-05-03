#include "RDKitInterface.hpp"
#include <iostream>

int main()
{
    RDKitMolecule mol("C1CCCCC1C(=O)N");
    std::cout << "SMILES: " << mol.toSmiles() << std::endl;
    return 0;
}