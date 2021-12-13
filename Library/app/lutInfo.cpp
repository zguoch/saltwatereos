#include "H2ONaCl.H"

int main(int argc, char** argv)
{
    if(argc!=2)ERROR("Usage: lutinfo myLUT.bin");
    H2ONaCl::cH2ONaCl sw;
    sw.loadLUT(argv[1]);
    return 0;
}