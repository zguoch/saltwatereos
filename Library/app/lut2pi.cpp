#include "H2ONaCl.H"
#include "LookUpTableForestI.H"

int main(int argc, char** argv)
{
    if(argc!=2)ERROR("Construct and write point index to binary file for a existed LUT file\nUsage: lutinfo myLUT.bin");
    H2ONaCl::cH2ONaCl sw;
    sw.loadLUT(argv[1]);
    if(sw.m_pLUT)
    {
        if(sw.m_dim_lut==2)((H2ONaCl::LookUpTableForest_2D*)sw.m_pLUT)->write_point_index(argv[1]);
        else ((H2ONaCl::LookUpTableForest_3D*)sw.m_pLUT)->write_point_index(argv[1]);
    }
    return 0;
}