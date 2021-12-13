#include "H2ONaCl.H"
#include <iostream>
H2ONaCl::cH2ONaCl eos;

void createTable_constX_TP(int max_level)
{
    int ind = 0;
    std::string dummy;
    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
#if USE_OMP == 1
    eos.set_num_threads(8);
#endif
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =2;
    double TP_min[2] = {1 + 273.15, 5E5}; //T [K], P[Pa]
    double TP_max[2] = {700 + 273.15, 400E5};
    double X_wt = 0.2; //wt% NaCl [0,1]
    int min_level = 4;
    // int max_level = 7;

    eos.createLUT_2D(TP_min, TP_max, X_wt, LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level);
    eos.save_lut_to_vtk("lut_constX_TP.vtu");
    eos.save_lut_to_binary("lut_constX_TP_"+std::to_string(max_level)+".bin");

    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E4;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double T_K = (rand()/(double)RAND_MAX)*(TP_max[0] - TP_min[0]) + TP_min[0];
        double p_Pa = (rand()/(double)RAND_MAX)*(TP_max[1] - TP_min[1]) + TP_min[1];
        eos.lookup(prop_lookup,T_K, p_Pa);
        prop_cal = eos.prop_pTX(p_Pa, T_K, X_wt);
        H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
        if(prop_lookup.Region != phaseRegion_cal)
        {
            // cout<<"Need refine point "<<ind++<<": ";
            // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
        }else
        {
            // double err = prop_cal.Rho-prop_lookup.Rho;
            // if(fabs(err)>0.5)
            // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
        }
    }
    STATUS_time("Searching done", clock() - start);

}
void createTable_constP_XT()
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
#if USE_OMP == 1
    eos.set_num_threads(8);
#endif
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =2;
    double Xmin = 1E-5, Xmax = 0.99999, Tmin = 1 + 273.15, Tmax = 1000 + 273.15;
    double XT_min[2] = {Xmin, Tmin}; //T [K], X[wt%: 0-1]
    double XT_max[2] = {Xmax, Tmax};
    double constP = 25E6; //Pa
    int min_level = 4;
    int max_level = 7;
    eos.createLUT_2D(Xmin, Xmax, Tmin, Tmax, constP, LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level);
    eos.save_lut_to_vtk("lut_constP_XT.vtu");
    eos.save_lut_to_binary("lut_constP_XT_"+std::to_string(max_level)+".bin");

    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E4;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double X_wt = (rand()/(double)RAND_MAX)*(XT_max[0] - XT_min[0]) + XT_min[0];
        double T_K = (rand()/(double)RAND_MAX)*(XT_max[1] - XT_min[1]) + XT_min[1];
        targetLeaf = eos.lookup(prop_lookup, X_wt, T_K);
        prop_cal = eos.prop_pTX(constP, T_K, X_wt);
        H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(constP, T_K, X_wt);
        if(prop_lookup.Region != phaseRegion_cal)
        {
            // cout<<"Need refine point "<<ind++<<": ";
            // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
        }else
        {
            // double err = prop_cal.Rho-prop_lookup.Rho;
            // if(fabs(err)>0.5)
            // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
        }
    }
    STATUS_time("Searching done", clock() - start);

}
void createTable_constP_XH()
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
#if USE_OMP == 1
    eos.set_num_threads(8);
#endif
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =2;
    double Xmin = 1E-5, Xmax = 0.99999; //Tmin = 1 + 273.15, Tmax = 1000 + 273.15;
    double Hmin = 0.1E6, Hmax = 3.5E6;
    double constP = 25E6; //Pa
    int min_level = 4;
    int max_level = 7;
    eos.createLUT_2D(Xmin, Xmax, Hmin, Hmax, constP, LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level);
    eos.save_lut_to_vtk("lut_constP_XH.vtu");
    eos.save_lut_to_binary("lut_constP_XH_"+std::to_string(max_level)+".bin");

    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E4;
    // for (size_t i = 0; i < n_randSample; i++)
    // {
    //     double X_wt = (rand()/(double)RAND_MAX)*(XT_max[0] - XT_min[0]) + XT_min[0];
    //     double T_K = (rand()/(double)RAND_MAX)*(XT_max[1] - XT_min[1]) + XT_min[1];
    //     targetLeaf = eos.lookup(prop_lookup, X_wt, T_K);
    //     prop_cal = eos.prop_pTX(constP, T_K, X_wt);
    //     H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(constP, T_K, X_wt);
    //     if(prop_lookup.Region != phaseRegion_cal)
    //     {
    //         // cout<<"Need refine point "<<ind++<<": ";
    //         // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
    //     }else
    //     {
    //         // double err = prop_cal.Rho-prop_lookup.Rho;
    //         // if(fabs(err)>0.5)
    //         // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
    //     }
    // }
    STATUS_time("Searching done", clock() - start);

}
void createTable_constT_XP(double constT)
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
#if USE_OMP == 1
    eos.set_num_threads(8);
#endif
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =2;
    double Xmin = 1E-5, Xmax = 0.99999, Pmin = 5E5, Pmax = 800E5;
    int min_level = 4;
    int max_level = 7;
    eos.createLUT_2D(Xmin, Xmax, Pmin, Pmax, constT, LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level);
    eos.save_lut_to_vtk("lut_constT_XP.vtu");
    eos.save_lut_to_binary("lut_constT_XP_"+std::to_string(max_level)+".bin");

    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E4;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double X_wt = (rand()/(double)RAND_MAX)*(Xmax - Xmin) + Xmin;
        double p_Pa = (rand()/(double)RAND_MAX)*(Pmax - Pmin) + Pmin;
        targetLeaf = eos.lookup(prop_lookup, X_wt, p_Pa);
        prop_cal = eos.prop_pTX(p_Pa, constT, X_wt);
        H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(p_Pa, constT, X_wt);
        if(prop_lookup.Region != phaseRegion_cal)
        {
            // cout<<"Need refine point "<<ind++<<": ";
            // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
        }else
        {
            // double err = prop_cal.Rho-prop_lookup.Rho;
            // if(fabs(err)>0.5)
            // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
        }
    }
    STATUS_time("Searching done", clock() - start);
}
void createTable_TPX()
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
#if USE_OMP == 1
    eos.set_num_threads(8);
#endif
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =3;
    double Tmin = 1 +273.15, Tmax = 1000+273.15, Xmin = 0.1, Xmax = 0.99999, Pmin = 5E5, Pmax = 2000E5;
    int min_level = 3;
    int max_level = 5;

    eos.createLUT_3D(Tmin, Tmax, Pmin, Pmax, Xmin, Xmax, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level);
    eos.save_lut_to_vtk("lut_TPX_"+std::to_string(max_level)+".vtu");
    eos.save_lut_to_binary("lut_TPX_"+std::to_string(max_level)+".bin");
    cout<<"\n测试读取并写入vtu: "<<endl;
    H2ONaCl::cH2ONaCl eos2;
    eos2.loadLUT("lut_TPX_"+std::to_string(max_level)+".bin");
    eos2.save_lut_to_vtk("lut_TPX_loadwrite.vtu");
    // STATUS("Start search ... ");
    // start = clock();
    // LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    // int n_randSample = 1E4;
    // for (size_t i = 0; i < n_randSample; i++)
    // {
    //     double X_wt = (rand()/(double)RAND_MAX)*(Xmax - Xmin) + Xmin;
    //     double p_Pa = (rand()/(double)RAND_MAX)*(Pmax - Pmin) + Pmin;
    //     // int ind_targetLeaf = forest.searchQuadrant(T_C, p_bar, 0);
    //     // eos.m_lut_PTX_2D->searchQuadrant(targetLeaf, T_K, p_Pa, X_wt);
    //     targetLeaf = eos.searchLUT_2D_PTX(prop_lookup, X_wt, p_Pa);
    //     prop_cal = eos.prop_pTX(p_Pa, constT, X_wt);
    //     H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(p_Pa, constT, X_wt);
    //     // cout<<"T_C: "<<T_K - 273.15<<", p_bar: "<<constP/1E5<<", index: "<<targetLeaf->index<<endl;
    //     // cout<<eos.getPhaseRegionName(targetLeaf->user_data->phaseRegion_cell)<<endl;
    //     // compare to directally calculation
    //     // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
    //     // if(targetLeaf->user_data->phaseRegion_cell != phaseRegion_cal)
    //     if(prop_lookup.Region != phaseRegion_cal)
    //     {
    //         // cout<<"Need refine point "<<ind++<<": ";
    //         // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
    //     }else
    //     {
    //         double err = prop_cal.Rho-prop_lookup.Rho;
    //         if(fabs(err)>0.5)
    //         cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
    //     }
    // }
    // // // forest.searchQuadrant(targetLeaf, 200, 300, 0);
    // // // cout<<targetLeaf->level<<endl;
    // STATUS_time("Searching done", clock() - start);

}
void createTable_HPX()
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
#if USE_OMP == 1
    eos.set_num_threads(8);
#endif
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =3;
    double Hmin = 0.1E6, Hmax = 3.9E6, Xmin = 0.001, Xmax = 1, Pmin = 100E5, Pmax = 2500E5;
    int min_level = 3;
    int max_level = 7;

    eos.createLUT_3D(Hmin, Hmax, Pmin, Pmax, Xmin, Xmax, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level);
    eos.save_lut_to_vtk("lut_HPX_"+std::to_string(max_level)+".vtu");
    eos.save_lut_to_binary("lut_HPX_"+std::to_string(max_level)+".bin");
    // cout<<"\n测试读取并写入vtu: "<<endl;
    // H2ONaCl::cH2ONaCl eos2;
    // eos2.loadLUT("lut_HPX_"+std::to_string(max_level)+".bin");
    // eos2.save_lut_to_vtk("lut_HPX_loadwrite.vtu");
    // STATUS("Start search ... ");
    // start = clock();
    // LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    // int n_randSample = 1E4;
    // for (size_t i = 0; i < n_randSample; i++)
    // {
    //     double X_wt = (rand()/(double)RAND_MAX)*(Xmax - Xmin) + Xmin;
    //     double p_Pa = (rand()/(double)RAND_MAX)*(Pmax - Pmin) + Pmin;
    //     // int ind_targetLeaf = forest.searchQuadrant(T_C, p_bar, 0);
    //     // eos.m_lut_PTX_2D->searchQuadrant(targetLeaf, T_K, p_Pa, X_wt);
    //     targetLeaf = eos.searchLUT_2D_PTX(prop_lookup, X_wt, p_Pa);
    //     prop_cal = eos.prop_pTX(p_Pa, constT, X_wt);
    //     H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(p_Pa, constT, X_wt);
    //     // cout<<"T_C: "<<T_K - 273.15<<", p_bar: "<<constP/1E5<<", index: "<<targetLeaf->index<<endl;
    //     // cout<<eos.getPhaseRegionName(targetLeaf->user_data->phaseRegion_cell)<<endl;
    //     // compare to directally calculation
    //     // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
    //     // if(targetLeaf->user_data->phaseRegion_cell != phaseRegion_cal)
    //     if(prop_lookup.Region != phaseRegion_cal)
    //     {
    //         // cout<<"Need refine point "<<ind++<<": ";
    //         // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
    //     }else
    //     {
    //         double err = prop_cal.Rho-prop_lookup.Rho;
    //         if(fabs(err)>0.5)
    //         cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
    //     }
    // }
    // // // forest.searchQuadrant(targetLeaf, 200, 300, 0);
    // // // cout<<targetLeaf->level<<endl;
    // STATUS_time("Searching done", clock() - start);

}
void load_binary_3d(string filename)
{
    int ind = 0;
    H2ONaCl::cH2ONaCl eos;
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    clock_t start = clock();

    eos.loadLUT(filename);
    cout<<"dim of the bin file: "<<eos.m_dim_lut<<endl;
    // eos.save_lut_to_vtk("lut_TPX.vtu");

    double Tmin = 1 +273.15, Tmax = 1000+273.15, Xmin = 0.1, Xmax = 0.99999, Pmin = 5E5, Pmax = 2000E5;

    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<3,LOOKUPTABLE_FOREST::FIELD_DATA<3> > *targetLeaf = NULL;
    int n_randSample = 1E6;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double T_K = (rand()/(double)RAND_MAX)*(Tmax - Tmin) + Tmin;
        double p_Pa = (rand()/(double)RAND_MAX)*(Pmax - Pmin) + Pmin;
        double X_wt = (rand()/(double)RAND_MAX)*(Xmax - Xmin) + Xmin;
        targetLeaf = eos.lookup(prop_lookup,T_K, p_Pa, X_wt); 
        // prop_cal = eos.prop_pTX(p_Pa, T_K, X_wt);
        // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(p_Pa, T_K, X_wt);
        // if(targetLeaf->user_data->prop_cell.Region != phaseRegion_cal)
        if(targetLeaf->user_data->need_refine)
        {
            ind++;
            // cout<<"Need refine point "<<ind++<<", level: "<<targetLeaf->level<<endl;
            // cout<<", rho: "<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
        }
        else
        {
            // double err = prop_cal.Rho-prop_lookup.Rho;
            // if(fabs(err)>0.5)
            // {
            //     cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
            //     // cout<<"  "<<targetLeaf->user_data->prop_point[0].Rho
            //     //     <<"  "<<targetLeaf->user_data->prop_point[1].Rho
            //     //     <<"  "<<targetLeaf->user_data->prop_point[2].Rho
            //     //     <<"  "<<targetLeaf->user_data->prop_point[3].Rho
            //     //     <<"  "<<targetLeaf->user_data->prop_cell.Rho
            //     //     <<", level: "<<targetLeaf->level<<", refine: "<<targetLeaf->user_data->need_refine
            //     //     <<endl;
            // }
        }
    }
    printf("All %d (%.2f %%) random points close to phase boundary.\n", ind, ind/(double)n_randSample*100);
    STATUS_time("Searching done", clock() - start);
}
void load_binary_2d(string filename)
{
    int ind = 0;
    H2ONaCl::cH2ONaCl eos;
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    clock_t start = clock();

    eos.loadLUT(filename);
    // ((H2ONaCl::LookUpTableForest_2D*)eos.m_pLUT)->print_summary();
    eos.save_lut_to_vtk("lut_load_save.vtu");

    // double Tmin = 1 +273.15, Tmax = 1000+273.15, Xmin = 0.1, Xmax = 0.99999, Pmin = 5E5, Pmax = 2000E5;

    // STATUS("Start search ... ");
    // start = clock();
    // LOOKUPTABLE_FOREST::Quadrant<3,LOOKUPTABLE_FOREST::FIELD_DATA<3> > *targetLeaf = NULL;
    // int n_randSample = 1E6;
    // for (size_t i = 0; i < n_randSample; i++)
    // {
    //     double T_K = (rand()/(double)RAND_MAX)*(Tmax - Tmin) + Tmin;
    //     double p_Pa = (rand()/(double)RAND_MAX)*(Pmax - Pmin) + Pmin;
    //     double X_wt = (rand()/(double)RAND_MAX)*(Xmax - Xmin) + Xmin;
    //     targetLeaf = eos.lookup(prop_lookup,T_K, p_Pa, X_wt); 
    //     // prop_cal = eos.prop_pTX(p_Pa, T_K, X_wt);
    //     // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(p_Pa, T_K, X_wt);
    //     // if(targetLeaf->user_data->prop_cell.Region != phaseRegion_cal)
    //     if(targetLeaf->user_data->need_refine)
    //     {
    //         ind++;
    //         // cout<<"Need refine point "<<ind++<<", level: "<<targetLeaf->level<<endl;
    //         // cout<<", rho: "<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
    //     }
    //     else
    //     {
    //         // double err = prop_cal.Rho-prop_lookup.Rho;
    //         // if(fabs(err)>0.5)
    //         // {
    //         //     cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
    //         //     // cout<<"  "<<targetLeaf->user_data->prop_point[0].Rho
    //         //     //     <<"  "<<targetLeaf->user_data->prop_point[1].Rho
    //         //     //     <<"  "<<targetLeaf->user_data->prop_point[2].Rho
    //         //     //     <<"  "<<targetLeaf->user_data->prop_point[3].Rho
    //         //     //     <<"  "<<targetLeaf->user_data->prop_cell.Rho
    //         //     //     <<", level: "<<targetLeaf->level<<", refine: "<<targetLeaf->user_data->need_refine
    //         //     //     <<endl;
    //         // }
    //     }
    // }
    // printf("All %d (%.2f %%) random points close to phase boundary.\n", ind, ind/(double)n_randSample*100);
    STATUS_time("Searching done", clock() - start);
}

int main(int argc, char** argv)
{
    // cout<<sizeof(LOOKUPTABLE_FOREST::Quadrant<3, LOOKUPTABLE_FOREST::FIELD_DATA<2>>**)<<endl;
    // std::map<LOOKUPTABLE_FOREST::Quad_index, double> test_map;
    // test_map[{1,2,3}] = 30;
    // test_map[{1,2,4}] = 40;
    // test_map[{2,2,4}] = 50;
    // test_map[{2,2,4}] = 70;
    // cout<<test_map[{1,2,3}]<<endl
    //     <<test_map[{1,2,4}]<<endl
    //     <<test_map[{2,2,4}]<<endl
    //     <<endl;

    if(argc<2)return 0;

    int ind = 0;
    std::string dummy;
    int test_index = atoi(argv[1]);
    switch (test_index)
    {
    case 1:
        {
            cout<<"createTable_constX_TP "<<endl;
            int max_level = atoi(argv[2]);
            createTable_constX_TP(max_level);
        }
        break;
    case 2:
        cout<<"createTable_constP_XT();"<<endl;
        // createTable_constP_XT();
        break;
    case 3:
        cout<<"createTable_constP_XH();"<<endl;
        // createTable_constP_XH();
        break;
    case 4: 
        cout<<"createTable_constT_XP(500+273.15);"<<endl;
        // createTable_constT_XP(500+273.15);
        break;
    case 5:
        cout<<"createTable_TPX();"<<endl;
        // createTable_TPX();
        break;
    case 6:
        cout<<"createTable_HPX();"<<endl;
        createTable_HPX();
        break;
    case 7:
        cout<<"load_binary_2d(\"lut_constP_XH_7.bin\");"<<endl;
        break;
    case 8:
        cout<<"load_binary_3d(\"lut_TPX_7.bin\");"<<endl;
        // load_binary_3d("lut_TPX_7.bin");
        break;
    default:
        break;
    }

    // destroy by hand
    // eos.destroyLUT_2D_PTX();

    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
}