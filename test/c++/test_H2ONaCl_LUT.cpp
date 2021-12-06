#include "H2ONaCl.H"
#include <iostream>
H2ONaCl::cH2ONaCl eos;

void createTable_constX_TP()
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
    eos.set_num_threads(8);
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =2;
    double TP_min[2] = {1 + 273.15, 5E5}; //T [K], P[Pa]
    double TP_max[2] = {700 + 273.15, 400E5};
    double X_wt = 0.2; //wt% NaCl [0,1]
    int min_level = 4;
    int max_level = 7;

    eos.createLUT_2D_TPX(TP_min, TP_max, X_wt, LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP, min_level, max_level);
    eos.save_lut_to_vtk("lut_PTX.vtu");
    
    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E4;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double T_K = (rand()/(double)RAND_MAX)*(TP_max[0] - TP_min[0]) + TP_min[0];
        double p_Pa = (rand()/(double)RAND_MAX)*(TP_max[1] - TP_min[1]) + TP_min[1];
        // int ind_targetLeaf = forest.searchQuadrant(T_C, p_bar, 0);
        // eos.m_lut_PTX_2D->searchQuadrant(targetLeaf, T_K, p_Pa, X_wt);
        eos.searchLUT_2D_PTX(prop_lookup,T_K, p_Pa);
        prop_cal = eos.prop_pTX(p_Pa, T_K, X_wt);
        H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
        // cout<<"T_C: "<<T_K - 273.15<<", p_bar: "<<p_Pa/1E5<<", index: "<<targetLeaf->index<<endl;
        // cout<<eos.getPhaseRegionName(targetLeaf->user_data->phaseRegion_cell)<<endl;
        // compare to directally calculation
        // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
        // if(targetLeaf->user_data->phaseRegion_cell != phaseRegion_cal)
        if(prop_lookup.Region != phaseRegion_cal)
        {
            // cout<<"Need refine point "<<ind++<<": ";
            // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
        }else
        {
            double err = prop_cal.Rho-prop_lookup.Rho;
            if(fabs(err)>0.5)
            cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
        }
    }
    // // forest.searchQuadrant(targetLeaf, 200, 300, 0);
    // // cout<<targetLeaf->level<<endl;
    STATUS_time("Searching done", clock() - start);

    // ======= test write to binary =====
    eos.save_lut_to_binary("lut_TPX_"+std::to_string(max_level)+".bin");
}
void createTable_constP_XT()
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
    eos.set_num_threads(8);
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =2;
    double Xmin = 1E-5, Xmax = 0.99999, Tmin = 1 + 273.15, Tmax = 1000 + 273.15;
    double XT_min[2] = {Xmin, Tmin}; //T [K], X[wt%: 0-1]
    double XT_max[2] = {Xmax, Tmax};
    double constP = 25E6; //Pa
    int min_level = 4;
    int max_level = 7;

    // eos.createLUT_2D_TPX(XT_min, XT_max, constP, LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH, min_level, max_level);
    eos.createLUT_2D_TPX(Xmin, Xmax, Tmin, Tmax, constP, LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH, min_level, max_level);
    eos.save_lut_to_vtk("lut_constP_XT.vtu");
    
    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E4;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double X_wt = (rand()/(double)RAND_MAX)*(XT_max[0] - XT_min[0]) + XT_min[0];
        double T_K = (rand()/(double)RAND_MAX)*(XT_max[1] - XT_min[1]) + XT_min[1];
        // int ind_targetLeaf = forest.searchQuadrant(T_C, p_bar, 0);
        // eos.m_lut_PTX_2D->searchQuadrant(targetLeaf, T_K, p_Pa, X_wt);
        targetLeaf = eos.searchLUT_2D_PTX(prop_lookup, X_wt, T_K);
        prop_cal = eos.prop_pTX(constP, T_K, X_wt);
        H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(constP, T_K, X_wt);
        // cout<<"T_C: "<<T_K - 273.15<<", p_bar: "<<constP/1E5<<", index: "<<targetLeaf->index<<endl;
        // cout<<eos.getPhaseRegionName(targetLeaf->user_data->phaseRegion_cell)<<endl;
        // compare to directally calculation
        // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
        // if(targetLeaf->user_data->phaseRegion_cell != phaseRegion_cal)
        if(prop_lookup.Region != phaseRegion_cal)
        {
            // cout<<"Need refine point "<<ind++<<": ";
            // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
        }else
        {
            double err = prop_cal.Rho-prop_lookup.Rho;
            if(fabs(err)>0.5)
            cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
        }
    }
    // // forest.searchQuadrant(targetLeaf, 200, 300, 0);
    // // cout<<targetLeaf->level<<endl;
    STATUS_time("Searching done", clock() - start);

    // // ======= test write to binary =====
    // eos.save_lut_to_binary("lut_TPX_"+std::to_string(max_level)+".bin");
}
void createTable_constT_XP(double constT)
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
    eos.set_num_threads(8);
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =2;
    double Xmin = 1E-5, Xmax = 0.99999, Pmin = 5E5, Pmax = 800E5;
    // double constT = 100+273.15; //K
    int min_level = 4;
    int max_level = 7;

    // eos.createLUT_2D_TPX(XT_min, XT_max, constP, LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH, min_level, max_level);
    eos.createLUT_2D_TPX(Xmin, Xmax, Pmin, Pmax, constT, LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP, min_level, max_level);
    eos.save_lut_to_vtk("lut_constT_XP.vtu");
    
    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E4;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double X_wt = (rand()/(double)RAND_MAX)*(Xmax - Xmin) + Xmin;
        double p_Pa = (rand()/(double)RAND_MAX)*(Pmax - Pmin) + Pmin;
        // int ind_targetLeaf = forest.searchQuadrant(T_C, p_bar, 0);
        // eos.m_lut_PTX_2D->searchQuadrant(targetLeaf, T_K, p_Pa, X_wt);
        targetLeaf = eos.searchLUT_2D_PTX(prop_lookup, X_wt, p_Pa);
        prop_cal = eos.prop_pTX(p_Pa, constT, X_wt);
        H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(p_Pa, constT, X_wt);
        // cout<<"T_C: "<<T_K - 273.15<<", p_bar: "<<constP/1E5<<", index: "<<targetLeaf->index<<endl;
        // cout<<eos.getPhaseRegionName(targetLeaf->user_data->phaseRegion_cell)<<endl;
        // compare to directally calculation
        // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
        // if(targetLeaf->user_data->phaseRegion_cell != phaseRegion_cal)
        if(prop_lookup.Region != phaseRegion_cal)
        {
            // cout<<"Need refine point "<<ind++<<": ";
            // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
        }else
        {
            double err = prop_cal.Rho-prop_lookup.Rho;
            if(fabs(err)>0.5)
            cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
        }
    }
    // // forest.searchQuadrant(targetLeaf, 200, 300, 0);
    // // cout<<targetLeaf->level<<endl;
    STATUS_time("Searching done", clock() - start);

    // // ======= test write to binary =====
    // eos.save_lut_to_binary("lut_TPX_"+std::to_string(max_level)+".bin");
}
void createTable_TPX()
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
    eos.set_num_threads(8);
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =3;
    double Tmin = 1 +273.15, Tmax = 1000+273.15, Xmin = 0.1, Xmax = 0.99999, Pmin = 5E5, Pmax = 2000E5;
    int min_level = 2;
    int max_level = 5;

    eos.createLUT_3D_TPX(Tmin, Tmax, Pmin, Pmax, Xmin, Xmax, min_level, max_level);
    eos.m_lut_PTX_3D->write_to_vtk("lut_TPX.vtu");
    eos.m_lut_PTX_3D->write_to_binary("lut_TPX.bin");
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
void load_binary(string filename)
{
    int ind = 0;
    H2ONaCl::cH2ONaCl eos;
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    clock_t start = clock();

    // const int dim = LOOKUPTABLE_FOREST::get_dim_from_binary(filename); 
    // LOOKUPTABLE_FOREST::LookUpTableForest<dim_in, LOOKUPTABLE_FOREST::FIELD_DATA<dim_in> > forest(filename);
    // load binary from EOS obj
    eos.loadLUT_PTX(filename);
    eos.m_lut_PTX_3D->write_to_vtk("lut_TPX.vtu");
    exit(0);
    double Tmin = 1 +273.15, Tmax = 1000+273.15, Xmin = 1E-5, Xmax = 0.99999, Pmin = 5E5, Pmax = 2000E5;

    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<3,LOOKUPTABLE_FOREST::FIELD_DATA<3> > *targetLeaf = NULL;
    int n_randSample = 1E5;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double T_K = (rand()/(double)RAND_MAX)*(Tmax - Tmin) + Tmin;
        double p_Pa = (rand()/(double)RAND_MAX)*(Pmax - Pmin) + Pmin;
        double X_wt = (rand()/(double)RAND_MAX)*(Xmax - Xmin) + Xmin;
        // int ind_targetLeaf = forest.searchQuadrant(T_C, p_bar, 0);
        eos.m_lut_PTX_3D->searchQuadrant(targetLeaf, T_K, p_Pa, X_wt);
        // targetLeaf = eos.searchLUT_2D_PTX(prop_lookup,T_K, p_Pa); 
        // prop_cal = eos.prop_pTX(p_Pa, T_K, X_wt);
        H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(p_Pa, T_K, X_wt);
        // cout<<"T_C: "<<T_K - 273.15<<", p_bar: "<<p_Pa/1E5<<", index: "<<targetLeaf->index<<endl;
        // cout<<eos.getPhaseRegionName(targetLeaf->user_data->phaseRegion_cell)<<endl;
        // compare to directally calculation
        // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
        // if(targetLeaf->user_data->phaseRegion_cell != phaseRegion_cal)
        if(targetLeaf->user_data->prop_cell.Region != phaseRegion_cal)
        {
            cout<<"Need refine point "<<ind++<<endl;
            // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<endl;
        }
        // else
        // {
        //     double err = prop_cal.Rho-prop_lookup.Rho;
        //     if(fabs(err)>0.5)
        //     {
        //         cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<err<<endl;
        //         cout<<"  "<<targetLeaf->user_data->prop_point[0].Rho
        //             <<"  "<<targetLeaf->user_data->prop_point[1].Rho
        //             <<"  "<<targetLeaf->user_data->prop_point[2].Rho
        //             <<"  "<<targetLeaf->user_data->prop_point[3].Rho
        //             <<"  "<<targetLeaf->user_data->prop_cell.Rho
        //             <<", level: "<<targetLeaf->level<<", refine: "<<targetLeaf->user_data->need_refine
        //             <<endl;
        //     }
        // }
    }
    // // forest.searchQuadrant(targetLeaf, 200, 300, 0);
    // // cout<<targetLeaf->level<<endl;
    STATUS_time("Searching done", clock() - start);
}
int main()
{
    // 1. 
    // createTable_constX_TP();
    // createTable_constP_XT();
    // createTable_constT_XP(500+273.15);
    createTable_TPX();

    // 2. 
    load_binary("lut_TPX.bin");

    // destroy by hand
    // eos.destroyLUT_2D_PTX();

    // std::cout << "Enter to continue..." << std::endl;
    // std::getline(std::cin, dummy);
}