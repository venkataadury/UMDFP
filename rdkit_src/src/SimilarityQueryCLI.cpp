//#define ENABLE_CPU_PARALLELISM 1
#include "UMFHelpers.h"
#include "RDKitInterface.hpp"
#include <iostream>

int main(int argc, char** argv)
{
    SimilarityQueryJobOutput* output = SimilarityQueryJob(std::string(argv[0])).run(argc, argv);
    return output->getExitCode();
}
