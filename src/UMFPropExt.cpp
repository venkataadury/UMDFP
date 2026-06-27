#include "JobCenter.hpp"

int main(int argc, char** argv)
{
    UMFPropertyExtractorJobOutput* output = UMFPropertyExtractorJob(std::string(argv[0])).run(argc, argv);
    return output->getExitCode();
}
