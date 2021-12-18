#include "H2ONaCl.H"

int main(int argc, char** argv)
{
    if(argc!=2)ERROR("Usage: lutinfo myLUT.bin");
    // H2ONaCl::cH2ONaCl sw;
    // sw.loadLUT(argv[1]);
    int m_dim_lut = LOOKUPTABLE_FOREST::get_dim_from_binary(argv[1]);
        
    switch (m_dim_lut)
    {
    case 2:
        {
            H2ONaCl::LookUpTableForest_2D lut_2d(argv[1], NULL);
        }
        break;
    case 3:
        {
            H2ONaCl::LookUpTableForest_3D lut_3d(argv[1], NULL);
        }
        break;
    default:
        ERROR("The dim in the binary file is neither 2 nor 3, it is not a valid LUT file: "+string(argv[1]));
        break;
    }
    return 0;
}