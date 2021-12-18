#include "H2ONaCl.H"

int main(int argc, char** argv)
{
    if(argc!=11)
    {
        STATUS("Usage: "+string(argv[0])+" [T|H] [T|H|P|X] [z0] [xmin] [xmax] [ymin] [ymax] [min_level] [max_level] [nThread]");
        STATUS("Unit: T[deg.C], P[bar], X[wt. NaCl, 0-1], H[MJ/kg]");
        STATUS_color("Example 1: T X 0.2 1 1000 5      600  4   7   8 (constX, var TP)", COLOR_BLUE);
        STATUS_color("Example 2: T P 200 0.001 1 1     1000 4   7   8 (constP, var XT)", COLOR_BLUE);
        STATUS_color("Example 3: T T 100 0.001 1 5     600  4   7   8 (constT, var XP)", COLOR_BLUE);
        STATUS_color("Example 5: H X 0.2 0.1 3.9 5     600  4   7   8 (constX, var HP)", COLOR_PURPLE);
        STATUS_color("Example 6: H P 200 0.001 1 0.1   3.9  4   7   8 (constP, var XH)", COLOR_PURPLE);
        STATUS_color("Example 7: H H 100 0.001 1 5     600  4   7   8 (constH, var XP)", COLOR_PURPLE);
        return 0;
    }
    string TorH = string(argv[1]);
    string const_var = string(argv[2]);
    double constZ = atof(argv[3]);
    double xmin = atof(argv[4]);
    double xmax = atof(argv[5]);
    double ymin = atof(argv[6]);
    double ymax = atof(argv[7]);
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
        if(const_var == "X")
        {
            eos.createLUT_2D(xmin + 273.15, xmax + 273.15, ymin * 1E5, ymax * 1E5, constZ, LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level, Update_prop_rho | Update_prop_drhodh | Update_prop_h);
            eos.save_lut_to_vtk("lut_constX_TP_"+std::to_string(max_level)+".vtu");
            eos.save_lut_to_binary("lut_constX_TP_"+std::to_string(max_level)+".bin");
        }else if(const_var == "P")
        {
            eos.createLUT_2D(xmin, xmax, ymin + 273.15, ymax+273.15, constZ * 1E5, LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level, Update_prop_rho | Update_prop_drhodh | Update_prop_h);
            eos.save_lut_to_vtk("lut_constP_XT_"+std::to_string(max_level)+".vtu");
            eos.save_lut_to_binary("lut_constP_XT_"+std::to_string(max_level)+".bin");
        }else if(const_var == "T")
        {
            eos.createLUT_2D(xmin, xmax, ymin*1E5, ymax*1E5, constZ + 273.15, LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level, Update_prop_rho | Update_prop_drhodh | Update_prop_h);
            eos.save_lut_to_vtk("lut_constT_XP_"+std::to_string(max_level)+".vtu");
            eos.save_lut_to_binary("lut_constT_XP_"+std::to_string(max_level)+".bin");
        }
        else
        {
            ERROR("If the first arg is T, the second arg must be one of T, P, X");
        }
    }else if (TorH == "H")
    {
        if(const_var == "X")
        {
            eos.createLUT_2D(xmin*1E6, xmax*1E6, ymin*1E5, ymax*1E5, constZ, LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level, Update_prop_rho | Update_prop_drhodh | Update_prop_h);
            eos.save_lut_to_vtk("lut_constX_HP_"+std::to_string(max_level)+".vtu");
            eos.save_lut_to_binary("lut_constX_HP_"+std::to_string(max_level)+".bin");
        }else if(const_var == "P")
        {
            eos.createLUT_2D(xmin, xmax, ymin*1E6, ymax*1E6, constZ*1E5, LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level, Update_prop_rho | Update_prop_drhodh | Update_prop_h);
            eos.save_lut_to_vtk("lut_constP_XH_"+std::to_string(max_level)+".vtu");
            eos.save_lut_to_binary("lut_constP_XH_"+std::to_string(max_level)+".bin");
        }else if(const_var == "H")
        {
            eos.createLUT_2D(xmin, xmax, ymin*1E5, ymax*1E5, constZ*1E6, LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level, Update_prop_rho | Update_prop_drhodh | Update_prop_h);
            eos.save_lut_to_vtk("lut_constH_XP_"+std::to_string(max_level)+".vtu");
            eos.save_lut_to_binary("lut_constH_XP_"+std::to_string(max_level)+".bin");
        }
        else
        {
            ERROR("If the first arg is H, the second arg must be one of H, P, X");
        }
    }else
    {
        ERROR("The first arg must be T or H.");
    }

    return 0;
}