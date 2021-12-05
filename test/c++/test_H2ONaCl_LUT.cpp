#include "H2ONaCl.H"
#include <iostream>
H2ONaCl::cH2ONaCl eos;

void createTable()
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    const int dim =2;
    double TP_min[2] = {1 + 273.15, 5E5}; //T [K], P[Pa]
    double TP_max[2] = {700 + 273.15, 400E5};
    double X_wt = 0.2; //wt% NaCl [0,1]
    int min_level = 4;
    int max_level = 7;

    eos.createLUT_2D_PTX("constX", TP_min, TP_max, X_wt, min_level, max_level);
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
void load_binary(string filename)
{
    int ind = 0;
    clock_t start = clock();
    const int dim = 2; 
    // LOOKUPTABLE_FOREST::LookUpTableForest<dim_in, LOOKUPTABLE_FOREST::FIELD_DATA<dim_in> > forest(filename);
    // load binary from EOS obj
    H2ONaCl::cH2ONaCl eos;
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    eos.loadLUT_PTX(filename);
    
    double TP_min[2] = {1 + 273.15, 5E5}; //T [K], P[Pa]
    double TP_max[2] = {700 + 273.15, 400E5};
    double X_wt = 0.2; //wt% NaCl [0,1]

    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E7;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double T_K = (rand()/(double)RAND_MAX)*(TP_max[0] - TP_min[0]) + TP_min[0];
        double p_Pa = (rand()/(double)RAND_MAX)*(TP_max[1] - TP_min[1]) + TP_min[1];
        // int ind_targetLeaf = forest.searchQuadrant(T_C, p_bar, 0);
        // eos.m_lut_PTX_2D->searchQuadrant(targetLeaf, T_K, p_Pa, X_wt);
        eos.searchLUT_2D_PTX(prop_lookup,T_K, p_Pa);
        // prop_cal = eos.prop_pTX(p_Pa, T_K, X_wt);
        // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
        // cout<<"T_C: "<<T_K - 273.15<<", p_bar: "<<p_Pa/1E5<<", index: "<<targetLeaf->index<<endl;
        // cout<<eos.getPhaseRegionName(targetLeaf->user_data->phaseRegion_cell)<<endl;
        // compare to directally calculation
        // H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
        // if(targetLeaf->user_data->phaseRegion_cell != phaseRegion_cal)
        // if(prop_lookup.Region != phaseRegion_cal)
        // {
        //     cout<<"Need refine point "<<ind++<<": ";
        //     // cout<<prop_cal.Rho<<"  ";
        //     cout<<prop_lookup.Rho<<endl;
        // }
        // else
        // {
        //     // cout<<prop_cal.Rho<<"  "<<prop_lookup.Rho<<", err: "<<prop_cal.Rho-prop_lookup.Rho<<endl;
        // }
    }
    // // forest.searchQuadrant(targetLeaf, 200, 300, 0);
    // // cout<<targetLeaf->level<<endl;
    STATUS_time("Searching done", clock() - start);
}
int main()
{
    // 1. 
    createTable();

    // 2. 
    // load_binary("lut_TPX_13.bin");

    // destroy by hand
    // eos.destroyLUT_2D_PTX();

    // std::cout << "Enter to continue..." << std::endl;
    // std::getline(std::cin, dummy);
}