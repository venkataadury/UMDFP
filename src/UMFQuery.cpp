#include "JobCenter.hpp"

int main(int argc, char** argv)
{
    UMFQueryJobOutput* output = UMFQueryJob(std::string(argv[0])).run(argc, argv);
    return output->getExitCode();
}
