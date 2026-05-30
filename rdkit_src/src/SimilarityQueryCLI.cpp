//#define ENABLE_CPU_PARALLELISM 1
#include "UMFHelpers.h"
#include "RDKitInterface.hpp"
#include <iostream>

int main(int argc, char** argv)
{
    SimilarityQueryJobOutput* output = SimilarityQueryJob(std::string(argv[0])).run(argc, argv);
    std::cout << "\n\n";
    std::cout << (*output) << "\n";
    return output->getExitCode();
}
