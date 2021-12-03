#include "LookUpTableForest.h"
using namespace LOOKUPTABLE_FOREST;
H2ONaCl::cH2ONaCl eos;

int phaseRegion_sin_func(double x, double y)
{
    double value = sin(x);
    int region = 0;
    if(value==y)
    {
        region = 0;
    }else if(y > value)
    {
        region = 1;
    }else
    {
        region = -1;
    }
    return region;
}

H2ONaCl::PhaseRegion phaseRegion_H2ONaCl_constantX(double T_C, double p_bar, double X = 3.2)
{
    double xv, xl;
    H2ONaCl::PhaseRegion phaseregion=eos.findPhaseRegion(T_C, p_bar, eos.Wt2Mol(X/100.0),xl,xv);

    return phaseregion;
}

void calProp_H2ONaCl_consX(H2ONaCl::PROP_H2ONaCl& prop, double T_C, double p_bar, double X = 3.2)
{
    prop = eos.prop_pTX(p_bar*1E5, T_C+273.15, X/100.0);
}

template <int dim, typename USER_DATA>
bool refine_fn(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level)
{
    if(quad->isHasChildren) return true; //if a quad has children, of course it need refine, but we don't need do anything at here, just return true.
    
    USER_DATA       *data = (USER_DATA *) quad->user_data;
    bool need_refine_phaseBoundary  = false;
    bool need_refine_Rho            = false;
    bool need_refine_H              = false;
    double physical_length[3];
    forest->get_quadrant_physical_length(quad->level, physical_length);
    const int num_sample_x =2; //最简单的情况就是只取xmin, xmax作为采样点判断这些采样点的函数计算返回值(flat)是否全部相等.但是有时候会有漏掉的情况，所以可以考虑在这里加密采样
    const int num_sample_y =2;
    double dx_qua = physical_length[0] / (num_sample_x - 1.0);
    double dy_qua = physical_length[1] / (num_sample_y - 1.0);
    double x_qua, x_ref, y_qua, y_ref;
    double xyz_tmp[3];
    H2ONaCl::PhaseRegion regionIndex[num_sample_x*num_sample_y];
    for (int iy = 0; iy < num_sample_y; iy++)
    {
        y_qua = quad->xyz[1] + dy_qua*iy;
        for (int ix = 0; ix < num_sample_x; ix++)
        {
            x_qua = quad->xyz[0] + dx_qua*ix;
            regionIndex[iy*num_sample_x + ix] = phaseRegion_H2ONaCl_constantX(x_qua, y_qua);
        }
    }
    // ========== 1. refinement check for phase index ============
    bool isSame_phaseIndex = true;
    for (int i = 1; i < num_sample_x*num_sample_y; i++)
    {
        isSame_phaseIndex = (isSame_phaseIndex && (regionIndex[0] == regionIndex[i]));
    }
    if(!isSame_phaseIndex) //if phase indices of all sample points are equal, do not refine; otherwise do refine
    {
        need_refine_phaseBoundary = true;
    }else
    {
        need_refine_phaseBoundary = false;
    }
    // ========================================================
    // calculate properties: four vertices and one midpoint
    calProp_H2ONaCl_consX(data->prop_point[0], quad->xyz[0],                      quad->xyz[1]);                         //xmin,ymin
    calProp_H2ONaCl_consX(data->prop_point[1], quad->xyz[0] + physical_length[0], quad->xyz[1]);                         //xmax,ymin
    calProp_H2ONaCl_consX(data->prop_point[2], quad->xyz[0],                      quad->xyz[1] + physical_length[1]);    //xmin,ymax
    calProp_H2ONaCl_consX(data->prop_point[3], quad->xyz[0] + physical_length[0], quad->xyz[1] + physical_length[1]);    //xmax,ymax
    calProp_H2ONaCl_consX(data->prop_cell, quad->xyz[0] + physical_length[0]/2.0, quad->xyz[1] + physical_length[1]/2.0); //xc,yc
    data->phaseRegion_point[0] = regionIndex[0]; //phase index
    data->phaseRegion_point[1] = regionIndex[num_sample_x-1];
    data->phaseRegion_point[2] = regionIndex[num_sample_x*num_sample_y-num_sample_x];
    data->phaseRegion_point[3] = regionIndex[num_sample_x*num_sample_y-1];
    data->phaseRegion_cell     = phaseRegion_H2ONaCl_constantX(quad->xyz[0] + physical_length[0]/2.0, quad->xyz[1] + physical_length[1]/2.0);
    // set some special indicator if cell need refine
    if(need_refine_phaseBoundary)
    {
        data->phaseRegion_cell = H2ONaCl::MixPhaseRegion;
        data->need_refine = NeedRefine_PhaseBoundary;
    }else
    {
        data->need_refine = NeedRefine_NoNeed;
    }

    // ========== 2. refinement check for Rho ===================== \todo maybe use another criterion
    double mean_Rho = data->prop_cell.Rho;
    for(int i=0;i<forest->m_num_children;i++)mean_Rho += data->prop_point[i].Rho;
    mean_Rho = mean_Rho / (forest->m_num_children + 1); // vertices data + one midpoint data
    double RMSD_Rho = pow(data->prop_cell.Rho - mean_Rho, 2.0);
    for(int i=0;i<forest->m_num_children;i++)RMSD_Rho += pow(data->prop_point[i].Rho - mean_Rho, 2.0);
    RMSD_Rho = sqrt(RMSD_Rho/(forest->m_num_children + 1));
    if(RMSD_Rho > forest->RMSD_Rho_min)
    {
        need_refine_Rho = true;
        if(data->need_refine == NeedRefine_NoNeed) data->need_refine = NeedRefine_Rho;
    }
    // // ========== 3. refinement check for H enthalpy ===================== \todo maybe use another criterion
    // double mean_H = data->prop_cell.H;
    // for(int i=0;i<forest->m_num_children;i++)mean_H += data->prop_point[i].H;
    // mean_H = mean_H / (forest->m_num_children + 1); // vertices data + one midpoint data
    // double RMSD_H = pow(data->prop_cell.H - mean_H, 2.0);
    // for(int i=0;i<forest->m_num_children;i++)RMSD_H += pow(data->prop_point[i].H - mean_H, 2.0);
    // RMSD_H = sqrt(RMSD_H/(forest->m_num_children + 1));
    // if(RMSD_H > forest->RMSD_H_min)
    // {
    //     need_refine_H = true;
    //     if(data->need_refine == NeedRefine_NoNeed) data->need_refine = NeedRefine_H;
    // }

    // ============ return refine indicator ===========
    if(quad->level > forest->m_max_level)return false;
    if(need_refine_phaseBoundary)return true;
    // if(need_refine_Rho) return true;
    // if(need_refine_H) return true;

    return false;
}

template <int dim, typename USER_DATA>
bool refine_uniform(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level)
{
    if(quad->level <= 3)return true;

    return false;
}

int main()
{
    clock_t start, end;
    double xyzmin[3] = {2,5, 0};
    double xyzmax[3] = {700, 400, 1};

    int max_level = 2;
    const int dim =2;
    LOOKUPTABLE_FOREST::LookUpTableForest<dim, LOOKUPTABLE_FOREST::FIELD_DATA<dim> > forest(xyzmin, xyzmax, max_level, sizeof(LOOKUPTABLE_FOREST::FIELD_DATA<dim>));
    // refine 
    start = clock();
    // forest.refine(refine_uniform);
    forest.refine(refine_fn);
    cout<<"Refinement end: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<" s"<<endl;
    
    start = clock();
    forest.write_to_vtk("lookuptable.vtu");
    cout<<"Write to vtk end: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<" s"<<endl;
    
    cout<<"search start ..."<<endl;
    start = clock();
    Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E7;
    // for (size_t i = 0; i < n_randSample; i++)
    // {
    //     double T_C = (rand()/(double)RAND_MAX)*(xyzmax[0] - xyzmin[0]) + xyzmin[0];
    //     double p_bar = (rand()/(double)RAND_MAX)*(xyzmax[1] - xyzmin[1]) + xyzmin[1];
    //     // int ind_targetLeaf = forest.searchQuadrant(T_C, p_bar, 0);
    //     forest.searchQuadrant(targetLeaf, T_C, p_bar, 0);
    //     // cout<<"T_C: "<<T_C<<", p_bar: "<<p_bar<<", index: "<<ind_targetLeaf<<endl;
    //     // cout<<eos.getPhaseRegionName(targetLeaf->user_data->phaseRegion_cell)<<endl;
    // }
    // forest.searchQuadrant(targetLeaf, 200, 300, 0);
    // cout<<targetLeaf->level<<endl;
    cout<<"Search "<<n_randSample<<" points: "<<(double)(clock()-start)/CLOCKS_PER_SEC<<" s"<<endl;
    
    // // std::string dummy;
    // // std::cout << "Enter to continue..." << std::endl;
    // // std::getline(std::cin, dummy);
    return 0;
}