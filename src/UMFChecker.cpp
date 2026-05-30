#include "JobCenter.hpp"

int main(int argc, char** argv)
{
    UMFCheckOutput* output = UMFCheckJob().run(argc, argv);
    return output->getExitCode();
}