/**
 * @file H2ONaCl_LUT_RefineFuncI.H
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief In order to easily manage the code file, collect all the refinement function for the H2ONaCl lookup table(LUT) to this file.
 * These refinement function can not be implement in the .cpp file, so just include this file in \file H2ONaCl.cpp
 * @version 0.1
 * @date 2021-12-03
 * 
 * @copyright Copyright (c) 2021
 * 
 */

/**
 * @brief Table refine function for the 2D case in T-P space with constant salinity.
 * 
 * @tparam dim 
 * @tparam USER_DATA 
 * @param forest 
 * @param quad 
 * @param max_level 
 * @return true 
 * @return false 
 */
template <int dim, typename USER_DATA>
bool RefineFunc_PTX_consX(LOOKUPTABLE_FOREST::LookUpTableForest<dim,USER_DATA>* forest, LOOKUPTABLE_FOREST::Quadrant<dim,USER_DATA>* quad, int max_level)
{
    // cout<<"level: "<<quad->level<<endl;
    if(quad->isHasChildren) return true; //if a quad has children, of course it need refine, but we don't need do anything at here, just return true.
    
    H2ONaCl::cH2ONaCl* eosPointer =(H2ONaCl::cH2ONaCl*)(forest->m_eosPointer); //read only! please DO NOT use this pointer to change any data in the EOS object!!! although it can change the member data. 
    double xv, xl;
    USER_DATA       *data = (USER_DATA *) quad->user_data;
    bool need_refine_phaseBoundary  = false;
    bool need_refine_Rho            = false;
    bool need_refine_H              = false;
    double physical_length[3];
    forest->get_quadrant_physical_length(quad->level, physical_length);
    const int num_sample_x =2; //if want to make more safe check, could increase this number to generate more sample points, for most of cases, the default value 2 is enough.
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
            // cout<<"T_C: "<<x_qua-273.15<<", p_bar: "<<y_qua/1E5<<", x_wt: "<<forest->m_constZ*100<<endl;
            regionIndex[iy*num_sample_x + ix] = eosPointer->findPhaseRegion(x_qua - 273.15, y_qua /1E5, eosPointer->Wt2Mol(forest->m_constZ),xl,xv); //note that the unit of T,P,X is deg.C, bar, and mol fraction in this function
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
    data->prop_point[0] = eosPointer->prop_pTX(quad->xyz[1],                            quad->xyz[0],                           forest->m_constZ);      //xmin,ymin
    data->prop_point[1] = eosPointer->prop_pTX(quad->xyz[1],                            quad->xyz[0] + physical_length[0],      forest->m_constZ);      //xmax,ymin
    data->prop_point[2] = eosPointer->prop_pTX(quad->xyz[1] + physical_length[1],       quad->xyz[0],                           forest->m_constZ);      //xmin,ymax
    data->prop_point[3] = eosPointer->prop_pTX(quad->xyz[1] + physical_length[1],       quad->xyz[0] + physical_length[0],      forest->m_constZ);      //xmax,ymax
    data->prop_cell     = eosPointer->prop_pTX(quad->xyz[1] + physical_length[1]/2.0,   quad->xyz[0] + physical_length[0]/2.0,  forest->m_constZ);      //xc,yc
    data->phaseRegion_point[0] = regionIndex[0]; //phase index
    data->phaseRegion_point[1] = regionIndex[num_sample_x-1];
    data->phaseRegion_point[2] = regionIndex[num_sample_x*num_sample_y-num_sample_x];
    data->phaseRegion_point[3] = regionIndex[num_sample_x*num_sample_y-1];
    data->phaseRegion_cell     = eosPointer->findPhaseRegion(quad->xyz[0] + physical_length[0]/2.0 - 273.15, (quad->xyz[1] + physical_length[1]/2.0)/1E5, eosPointer->Wt2Mol(forest->m_constZ),xl,xv); //note that the unit of T,P,X is deg.C, bar, and mol fraction in this function
    // set some special indicator if cell need refine
    if(need_refine_phaseBoundary)
    {
        data->phaseRegion_cell = H2ONaCl::MixPhaseRegion;
        data->need_refine = LOOKUPTABLE_FOREST::NeedRefine_PhaseBoundary;
    }else
    {
        data->need_refine = LOOKUPTABLE_FOREST::NeedRefine_NoNeed;
    }

    // ========== 2. refinement check for Rho ===================== \todo maybe use another criterion
    double mean_Rho = data->prop_cell.Rho;
    for(int i=0;i<forest->m_num_children;i++)mean_Rho += data->prop_point[i].Rho;
    mean_Rho = mean_Rho / (forest->m_num_children + 1); // vertices data + one midpoint data
    double RMSD_Rho = pow(data->prop_cell.Rho - mean_Rho, 2.0);
    for(int i=0;i<forest->m_num_children;i++)RMSD_Rho += pow(data->prop_point[i].Rho - mean_Rho, 2.0);
    RMSD_Rho = sqrt(RMSD_Rho/(forest->m_num_children + 1));
    if(RMSD_Rho > forest->m_RMSD_RefineCriterion.Rho)
    {
        need_refine_Rho = true;
        if(data->need_refine == LOOKUPTABLE_FOREST::NeedRefine_NoNeed) data->need_refine = LOOKUPTABLE_FOREST::NeedRefine_Rho;
    }
    // // ========== 3. refinement check for H enthalpy ===================== \todo maybe use another criterion
    // double mean_H = data->prop_cell.H;
    // for(int i=0;i<forest->m_num_children;i++)mean_H += data->prop_point[i].H;
    // mean_H = mean_H / (forest->m_num_children + 1); // vertices data + one midpoint data
    // double RMSD_H = pow(data->prop_cell.H - mean_H, 2.0);
    // for(int i=0;i<forest->m_num_children;i++)RMSD_H += pow(data->prop_point[i].H - mean_H, 2.0);
    // RMSD_H = sqrt(RMSD_H/(forest->m_num_children + 1));
    // if(RMSD_H > forest->m_RMSD_RefineCriterion.H)
    // {
    //     need_refine_H = true;
    //     if(data->need_refine == LOOKUPTABLE_FOREST::NeedRefine_NoNeed) data->need_refine = LOOKUPTABLE_FOREST::NeedRefine_H;
    // }

    // ============ return refine indicator ===========
    if(quad->level > forest->m_max_level)return false;
    if(need_refine_phaseBoundary)return true;
    if(need_refine_Rho) return true;
    // if(need_refine_H) return true;

    return false;
}


/**
 * @brief Uniform refinement. Used to define the minimum level mesh.
 * 
 * @tparam dim 
 * @tparam USER_DATA 
 * @param forest 
 * @param quad 
 * @param max_level 
 * @return true 
 * @return false 
 */
template <int dim, typename USER_DATA>
static bool refine_uniform(LOOKUPTABLE_FOREST::LookUpTableForest<dim,USER_DATA>* forest, LOOKUPTABLE_FOREST::Quadrant<dim,USER_DATA>* quad, int max_level)
{
    if(quad->level < forest->m_min_level)return true;

    return false;
};