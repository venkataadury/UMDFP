#include "JobCenter.hpp"

static std::ofstream outfile;
int main(int argc, char** argv)
{
    UMFJobOutput* output = UMFDumpJob().run(argc, argv);
    return output->getExitCode();
}
