#include "H2ONaCl.H"
#include <iostream>
H2ONaCl::cH2ONaCl eos;

int main()
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
    const int dim =2;
    double TP_min[2] = {1 + 273.15, 5E5}; //T [K], P[Pa]
    double TP_max[2] = {700 + 273.15, 400E5};
    double X_wt = 0.032; //wt% NaCl [0,1]
    int min_level = 4;
    int max_level = 7;

    eos.createLUT_2D_PTX("constX", TP_min, TP_max, X_wt, min_level, max_level, "lut_PTX.vtu");
    
    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E4;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double T_K = (rand()/(double)RAND_MAX)*(TP_max[0] - TP_min[0]) + TP_min[0];
        double p_Pa = (rand()/(double)RAND_MAX)*(TP_max[1] - TP_min[1]) + TP_min[1];
        // int ind_targetLeaf = forest.searchQuadrant(T_C, p_bar, 0);
        eos.m_lut_PTX->searchQuadrant(targetLeaf, T_K, p_Pa, X_wt);
        // cout<<"T_C: "<<T_K - 273.15<<", p_bar: "<<p_Pa/1E5<<", index: "<<targetLeaf->index<<endl;
        // cout<<eos.getPhaseRegionName(targetLeaf->user_data->phaseRegion_cell)<<endl;
        // compare to directally calculation
        H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
        if(targetLeaf->user_data->phaseRegion_cell != phaseRegion_cal)
        {
            cout<<"Need refine point "<<ind++<<": "<<targetLeaf->user_data->need_refine<<endl;
        }
    }
    // // forest.searchQuadrant(targetLeaf, 200, 300, 0);
    // // cout<<targetLeaf->level<<endl;
    STATUS_time("Searching done", clock() - start);

    // eos.destroyLUT_2D_PTX();
    // eos.m_lut_PTX->destroy();

    // std::cout << "Enter to continue..." << std::endl;
    // std::getline(std::cin, dummy);
}