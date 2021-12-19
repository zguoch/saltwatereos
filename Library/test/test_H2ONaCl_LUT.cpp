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

    eos.createLUT_2D(TP_min, TP_max, X_wt, LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level, Update_prop_rho | Update_prop_drhodh | Update_prop_h);
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
void createTable_TPX(int max_level)
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
    // int max_level = 5;

    eos.createLUT_3D(Tmin, Tmax, Pmin, Pmax, Xmin, Xmax, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level, Update_prop_rho | Update_prop_drhodh | Update_prop_h);
    eos.save_lut_to_vtk("lut_TPX_"+std::to_string(max_level)+".vtu");
    eos.save_lut_to_binary("lut_TPX_"+std::to_string(max_level)+".bin");
    // cout<<"\n测试读取并写入vtu: "<<endl;
    // H2ONaCl::cH2ONaCl eos2;
    // eos2.loadLUT("lut_TPX_"+std::to_string(max_level)+".bin");
    // eos2.save_lut_to_vtk("lut_TPX_loadwrite.vtu");
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
void createTable_HPX(int max_level)
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
    // int max_level = 7;

    eos.createLUT_3D(Hmin, Hmax, Pmin, Pmax, Xmin, Xmax, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level, Update_prop_rho | Update_prop_h | Update_prop_drhodh);
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
    const int dim = 3;
    eos.loadLUT(filename);
    cout<<"dim of the bin file: "<<eos.m_dim_lut<<endl;
    // eos.save_lut_to_vtk("lut_TPX.vtu");
    H2ONaCl::LookUpTableForest_3D* pLUT = (H2ONaCl::LookUpTableForest_3D*)eos.m_pLUT;

    double TorHmin = pLUT->m_xyz_min[0], TorHmax = pLUT->m_xyz_max[0], Pmin = pLUT->m_xyz_min[1], Pmax = pLUT->m_xyz_max[1];
    double Xmin = pLUT->m_xyz_min[2], Xmax = pLUT->m_xyz_max[2];

    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    double* props = new double[pLUT->m_map_props.size()];
    int ind_rho = 0;
    int n_randSample = 10;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double TorH = (rand()/(double)RAND_MAX)*(TorHmax - TorHmin) + TorHmin;
        double p_Pa = (rand()/(double)RAND_MAX)*(Pmax - Pmin) + Pmin;
        double X_wt = (rand()/(double)RAND_MAX)*(Xmax - Xmin) + Xmin;
        targetLeaf = eos.lookup(props, TorH, p_Pa, X_wt, false); 
        if(pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)
        {
            prop_cal = eos.prop_pTX(p_Pa, TorH, X_wt);
        }else
        {
            prop_cal = eos.prop_pHX(p_Pa, TorH, X_wt);
        }
        // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(p_Pa, T_K, X_wt);
        // if(targetLeaf->user_data->prop_cell.Region != phaseRegion_cal)
        if(targetLeaf->user_data->need_refine)
        {
            ind++;
            cout<<"Need refine point "<<ind++<<", level: "<<targetLeaf->level;
            cout<<", rho: "<<prop_cal.Rho<<"  "<<props[ind_rho]<<endl;
        }
        else
        {
            cout<<"Rho: "<<prop_cal.Rho<<"  "<<props[ind_rho]<<endl;
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

    delete[] props;

    STATUS_time("Searching done", clock() - start);
}
void load_binary_2d(string filename)
{
    int ind = 0;
    H2ONaCl::cH2ONaCl eos;
    H2ONaCl::PROP_H2ONaCl prop_cal;
    clock_t start = clock();
    const int dim = 2;
    int dim_file = LOOKUPTABLE_FOREST::get_dim_from_binary(filename);
    if(dim!=dim_file)ERROR("The input LUT file is not 2D.");
    eos.loadLUT(filename);
    H2ONaCl::LookUpTableForest_2D* pLUT = (H2ONaCl::LookUpTableForest_2D*)eos.m_pLUT;
    // eos.save_lut_to_vtk("lut_load_save.vtu");
    double Xmin = pLUT->m_xyz_min[0], Xmax = pLUT->m_xyz_max[0], Ymin = pLUT->m_xyz_min[1], Ymax = pLUT->m_xyz_max[1];
    double constZ = pLUT->m_constZ;
    int index_rho = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_rho));
    int index_h = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_h));
    int index_T = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_T));

    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E3;
    double* props = new double[pLUT->m_map_props.size()];
    int ind_rho = 0;
    int nx = 100, ny = 150;
    double dx = (Xmax - Xmin)/(nx - 1);
    double dy = (Ymax - Ymin)/(ny - 1);
    double x, y;
    ofstream fpout_xx("XX.txt");
    ofstream fpout_yy("YY.txt");
    ofstream fpout_zz("ZZ.txt");
    ofstream fpout_zz_lut("ZZ_lut.txt");
    bool isCal = true;
    for (int i = 0; i < ny; i++)
    {
        y = Ymin + i*dy;
        for (int j = 0; j < nx; j++)
        {
            x = Xmin + j*dx;
            if(isCal)
            {
                switch (pLUT->m_TorH)
                {
                case LOOKUPTABLE_FOREST::EOS_ENERGY_T:
                    {
                        switch (pLUT->m_const_which_var)
                        {
                        case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
                            prop_cal = eos.prop_pTX(y, x, constZ);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
                            prop_cal = eos.prop_pTX(constZ, y, x);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
                            prop_cal = eos.prop_pTX(y, constZ, x);
                            break;
                        default:
                            break;
                        }
                    }
                    break;
                case LOOKUPTABLE_FOREST::EOS_ENERGY_H:
                    {
                        switch (pLUT->m_const_which_var)
                        {
                        case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
                            prop_cal = eos.prop_pHX(y, x, constZ);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
                            prop_cal = eos.prop_pHX(constZ, y, x);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
                            prop_cal = eos.prop_pHX(y, constZ, x);
                            break;
                        default:
                            break;
                        }
                    }
                    break;
                default:
                    break;
                }
                fpout_zz<<prop_cal.Rho<<" ";
            }
            
            targetLeaf = eos.lookup(props, x, y, false);
            fpout_xx<<x<<" ";
            fpout_yy<<y<<" ";
            fpout_zz_lut<<props[index_rho]<<" ";
        }
        fpout_xx<<endl;
        fpout_yy<<endl;
        fpout_zz<<endl;
        fpout_zz_lut<<endl;
    }
    fpout_xx.close();
    fpout_yy.close();
    fpout_zz.close();
    fpout_zz_lut.close();
    
    // for (size_t i = 0; i < n_randSample; i++)
    // {
    //     double T_K = (rand()/(double)RAND_MAX)*(Tmax - Tmin) + Tmin;
    //     double p_Pa = (rand()/(double)RAND_MAX)*(Pmax - Pmin) + Pmin;
    //     // double X_wt = (rand()/(double)RAND_MAX)*(Xmax - Xmin) + Xmin;
    //     // targetLeaf = eos.lookup(prop_lookup,T_K, p_Pa); 
    //     targetLeaf = eos.lookup(props, T_K, p_Pa, true);
    //     prop_cal = eos.prop_pTX(p_Pa, T_K, x_const);
    //     // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion_pTX(p_Pa, T_K, X_wt);
    //     // if(targetLeaf->user_data->prop_cell.Region != phaseRegion_cal)
    //     // if(targetLeaf->user_data->need_refine)
    //     {
    //         // ind++;
    //         // cout<<"Need refine point "<<ind++<<", level: "<<(int)targetLeaf->level<<endl;
    //         // // cout<<", rho: "<<props[distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_rho))]<<endl;
    //         // cout<<"  Rho: "<<props[0] //distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_rho))
    //         //     <<" "
    //         //     <<prop_cal.Rho
    //         //     <<endl;
    //     }
    //     // else
    //     {
    //         // cout<<"Rho: "<<props[ind_rho]
    //         //     <<" "
    //         //     <<prop_cal.Rho
    //         //     <<endl;
    //         // double err = prop_cal.Rho-props[ind_rho];
    //         // if(fabs(err)>pLUT->m_RMSD_RefineCriterion.Rho)
    //         // {
    //         //     cout<<prop_cal.Rho<<"  "<<props[ind_rho]<<", err: "<<err<<endl;
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

    delete[] props;

    STATUS_time("Searching done", clock() - start);
}
void help(char** argv)
{
    cout<<"Help information ... "<<endl;
    cout<<argv[0]<<" 1 [max_level]: createTable_constX_TP"<<endl;
    cout<<argv[0]<<" 2 [max_level]: createTable_constP_XT"<<endl;
    cout<<argv[0]<<" 3 [max_level]: createTable_constP_XH"<<endl;
    cout<<argv[0]<<" 4 [max_level]: createTable_constT_XP"<<endl;
    cout<<argv[0]<<" 5 [max_level]: createTable_TPX"<<endl;
    cout<<argv[0]<<" 6 [max_level]: createTable_HPX"<<endl;
    cout<<argv[0]<<" 7 [my.bin]: load_binary_2d consX_TP"<<endl;
    cout<<argv[0]<<" 8 [my.bin]: load_binary_3d"<<endl;

    exit(0);
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
    if(argc==1)help(argv);

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
        {
            int max_level = atoi(argv[2]);
            cout<<"createTable_TPX();"<<endl;
            createTable_TPX(max_level);
        }
        break;
    case 6:
        {
            int max_level = atoi(argv[2]);
            cout<<"createTable_HPX();"<<endl;
            createTable_HPX(max_level);
        }
        break;
    case 7:
        cout<<"Test loading 2D LUT: "<<argv[2]<<endl;
        load_binary_2d(argv[2]);
        break;
    case 8:
        cout<<"Test loading 3D LUT: "<<argv[2]<<endl;
        load_binary_3d(argv[2]);
        break;
    default:
        break;
    }

    // destroy by hand
    // eos.destroyLUT_2D_PTX();

    std::cout << "Enter to continue..." << std::endl;
    std::getline(std::cin, dummy);
}