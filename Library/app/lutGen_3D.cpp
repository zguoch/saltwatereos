#include "H2ONaCl.H"

int main(int argc, char** argv)
{
    if(argc!=11)
    {
        STATUS("Usage: "+string(argv[0])+" [T|H] [xmin] [xmax] [ymin] [ymax] [zmin] [zmax] [min_level] [max_level] [nThread]");
        STATUS("Unit: T[deg.C], P[bar], X[wt. NaCl, 0-1], H[MJ/kg]");
        STATUS_color("Example 1: T 1    1000  5 600 0.001 1  4   7   8 (TPX)", COLOR_BLUE);
        STATUS_color("Example 2: H 0.1  3.9   5 600 0.001 1  4   7   8 (HPX)", COLOR_PURPLE);
        return 0;
    }
    string TorH = string(argv[1]);
    double xmin = atof(argv[2]);
    double xmax = atof(argv[3]);
    double ymin = atof(argv[4]);
    double ymax = atof(argv[5]);
    double zmin = atof(argv[6]);
    double zmax = atof(argv[7]);
    int min_level = atoi(argv[8]);
    int max_level = atoi(argv[9]);
    int n_threads = atoi(argv[10]);
    H2ONaCl::cH2ONaCl eos;
    #if USE_OMP == 1
        eos.set_num_threads(n_threads);
    #endif
    if(min_level > max_level)
    {
        WARNING("min_level > max_level: " + to_string(min_level)+" "+to_string(max_level));
        int tmp = min_level;
        min_level = max_level;
        max_level = tmp;
    }
    if(min_level<0)
    {
        WARNING("min_level<0: "+ to_string(min_level));
        min_level = 0;
    }
    if(max_level>(MAX_FOREST_LEVEL-3))
    {
        WARNING("min_level is too big: " + to_string(min_level));
        min_level = MAX_FOREST_LEVEL - 3;
    }
    if (TorH == "T")
    {
        eos.createLUT_3D(xmin + 273.15, xmax + 273.15, ymin * 1E5, ymax * 1E5, xmin, xmax, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level, Update_prop_rho | Update_prop_h);
        eos.save_lut_to_vtk("lut_TPX_"+std::to_string(max_level)+".vtu");
        eos.save_lut_to_binary("lut_TPX_"+std::to_string(max_level)+".bin");
    }else if (TorH == "H")
    {
        eos.createLUT_3D(xmin*1E6, xmax*1E6, ymin*1E5, ymax*1E5, zmin, zmax, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level, Update_prop_rho | Update_prop_h);
        eos.save_lut_to_vtk("lut_HPX_"+std::to_string(max_level)+".vtu");
        eos.save_lut_to_binary("lut_HPX_"+std::to_string(max_level)+".bin");
    }else
    {
        ERROR("The first arg must be T or H.");
    }

    return 0;
}