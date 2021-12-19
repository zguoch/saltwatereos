// generate random points in PTX/HPX range of a LUT

#include "H2ONaCl.H"

int main(int argc, char** argv)
{
    if(argc!=3)ERROR("Usage: "+string(argv[0])+" myLUT.bin 1000");
    
    
    int m_dim_lut = LOOKUPTABLE_FOREST::get_dim_from_binary(argv[1]);
    int n_randSample = atoi(argv[2]);
    double x, y, z;
    switch (m_dim_lut)
    {
    case 2:
        {
            H2ONaCl::LookUpTableForest_2D lut_2d(argv[1], NULL);
            FILE* fp = NULL;
            fp = fopen("Random.xyz", "w");
            if(!fp)ERROR("Open file failed: random.xyz");
            for (size_t i = 0; i < n_randSample; i++)
            {
                x = (rand()/(double)RAND_MAX)*(lut_2d.m_xyz_max[0] - lut_2d.m_xyz_min[0]) + lut_2d.m_xyz_min[0];
                y = (rand()/(double)RAND_MAX)*(lut_2d.m_xyz_max[1] - lut_2d.m_xyz_min[1]) + lut_2d.m_xyz_min[1];
                fprintf(fp, "%f %f %f\n", x, y, lut_2d.m_constZ);
            }
            fclose(fp);
        }
        break;
    case 3:
        {
            H2ONaCl::LookUpTableForest_3D lut_3d(argv[1], NULL);
            FILE* fp = NULL;
            fp = fopen("Random.xyz", "w");
            if(fp == NULL)ERROR("Open file failed: Random.xyz");
            for (size_t i = 0; i < n_randSample; i++)
            {
                x = (rand()/(double)RAND_MAX)*(lut_3d.m_xyz_max[0] - lut_3d.m_xyz_min[0]) + lut_3d.m_xyz_min[0];
                y = (rand()/(double)RAND_MAX)*(lut_3d.m_xyz_max[1] - lut_3d.m_xyz_min[1]) + lut_3d.m_xyz_min[1];
                z = (rand()/(double)RAND_MAX)*(lut_3d.m_xyz_max[2] - lut_3d.m_xyz_min[2]) + lut_3d.m_xyz_min[2];
                fprintf(fp, "%f %f %f\n", x, y, z);
            }
            fclose(fp);
        }
        break;
    default:
        ERROR("The dim in the binary file is neither 2 nor 3, it is not a valid LUT file: "+string(argv[1]));
        break;
    }
    return 0;
}