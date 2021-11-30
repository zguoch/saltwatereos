/**
 * @file step5.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief 尝试对H2ONaCl的二位相图进行自适应网格细化，并且输出phase id。X=0.032, T in [1, 1000] C, P in [5, 500] bar
 * @version 0.1
 * @date 2021-11-30
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <p4est_vtk.h> //用于输出vtk格式
#include <p4est_extended.h> //用于直接调用p4est_new_ext函数创建tree
#include "H2ONaCl.H"  //加入H2ONaCl动态库头文件
H2ONaCl::cH2ONaCl eos;

#include <iostream>
using namespace std;
#define MAX_Level_Refine 10
#define PI 3.141592653

struct FIELD_DATA
{
    bool need_refine; 
    int phaseRegion_point[P4EST_CHILDREN], phaseRegion_cell;
    H2ONaCl::PROP_H2ONaCl prop_point[P4EST_CHILDREN]; // properties at vertiex
    H2ONaCl::PROP_H2ONaCl prop_cell; // properties at midpoint as cell value (for vtk output)
    bool need_refine_Rho;
};

int phaseRegion_quadr_func(double x_ref, double y_ref, double xmin = 0, double xmax = 1, double ymin =0, double ymax = 0.25)
{
    // convert coordinate x of reference system to physical system
    double lenx = xmax - xmin;
    double leny = ymax - ymin;
    double x = lenx * x_ref + xmin;
    double y = leny * y_ref + ymin;
    double value = pow(x-0.5, 2.0);
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
int phaseRegion_sin_func(double x_ref, double y_ref, double xmin = 0, double xmax = 2*PI, double ymin =-1, double ymax = 1)
{
    // convert coordinate x of reference system to physical system
    double lenx = xmax - xmin;
    double leny = ymax - ymin;
    double x = lenx * x_ref + xmin;
    double y = leny * y_ref + ymin;
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
int phaseRegion_H2ONaCl_constantX(double x_ref, double y_ref, double X = 3.2, double Tmin = 1, double Tmax = 700, double pmin = 5, double pmax = 400)
{
    // T: x axis; p: y axis
    // convert coordinate x of reference system to physical system
    double lenx = Tmax - Tmin;
    double leny = pmax - pmin;
    double T_C = lenx * x_ref + Tmin;
    double p_bar = leny * y_ref + pmin;
    // calculate phase region index
    double xv, xl;
    H2ONaCl::PhaseRegion phaseregion=eos.findPhaseRegion(T_C, p_bar, eos.Wt2Mol(X/100.0),xl,xv);
    // cout<<phaseregion<<": "<<eos.m_phaseRegion_name[phaseregion].c_str()<<endl;
    // printf("T: %f C, P: %f bar, X: %f wt%% NaCl %f mol fraction\nphase region: %s xl: %f wt%% %f mol fraction, xv: %f wt%% %f mol fraction \n", T, P, X,eos.Wt2Mol(X/100.0),eos.m_phaseRegion_name[phaseregion].c_str(), eos.Mol2Wt(xl)*100.0, xl, eos.Mol2Wt(xv)*100.0, xv);

    return phaseregion;
}
void calProp_H2ONaCl(H2ONaCl::PROP_H2ONaCl& prop,double x_ref, double y_ref, double X = 3.2, double Tmin = 1, double Tmax = 700, double pmin = 5, double pmax = 400)
{
    double lenx = Tmax - Tmin;
    double leny = pmax - pmin;
    double T_C = lenx * x_ref + Tmin;
    double p_bar = leny * y_ref + pmin;
    prop = eos.prop_pTX(p_bar*1E5, T_C+273.15, X/100.0);
}
static int refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    FIELD_DATA       *data = (FIELD_DATA *) q->p.user_data;
    p4est_qcoord_t      half_length = P4EST_QUADRANT_LEN (q->level) / 2;
    p4est_qcoord_t      length = P4EST_QUADRANT_LEN (q->level);
    const int num_sample_x =2; //最简单的情况就是只取xmin, xmax作为采样点判断这些采样点的函数计算返回值(flat)是否全部相等.但是有时候会有漏掉的情况，所以可以考虑在这里加密采样
    const int num_sample_y =2;
    double dx_qua = P4EST_QUADRANT_LEN (q->level) / (num_sample_x - 1.0);
    double dy_qua = P4EST_QUADRANT_LEN (q->level) / (num_sample_y - 1.0);
    double x_qua, x_ref, y_qua, y_ref;
    double xyz_tmp[3];
    int regionIndex[num_sample_x*num_sample_y];
    for (int iy = 0; iy < num_sample_y; iy++)
    {
        y_qua = q->y + dy_qua*iy;
        for (int ix = 0; ix < num_sample_x; ix++)
        {
            x_qua = q->x + dx_qua*ix;
            p4est_qcoord_to_vertex (p4est->connectivity, which_tree, x_qua, y_qua, xyz_tmp);
            regionIndex[iy*num_sample_x + ix] = phaseRegion_H2ONaCl_constantX(xyz_tmp[0], xyz_tmp[1]);
        }
    }
    // ========== 1. refinement check for phase index ============
    bool isSame_phaseIndex = true;
    for (int i = 1; i < num_sample_x*num_sample_y; i++)
    {
        isSame_phaseIndex = (isSame_phaseIndex && (regionIndex[0] == regionIndex[i]));
    }
    if(!isSame_phaseIndex) //如果采样点的phase index全相等，则不refine；否则refine
    {
        data->need_refine = true;
    }else
    {
        data->need_refine = false;
    }
    // ========================================================
    // calculate properties: four vertices and one midpoint
    double xyzmin[3], xyzmax[3],xyzc[3];
    // cal min,max,center coordinate
    p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x, q->y, xyzmin);
    p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x + length, q->y + length, xyzmax);
    p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x + half_length, q->y + half_length, xyzc);
    // calculate props and give it to user data
    calProp_H2ONaCl(data->prop_point[0], xyzmin[0], xyzmin[1]); //xmin,ymin
    calProp_H2ONaCl(data->prop_point[1], xyzmax[0], xyzmin[1]); //xmax,ymin
    calProp_H2ONaCl(data->prop_point[2], xyzmin[0], xyzmax[1]); //xmin,ymax
    calProp_H2ONaCl(data->prop_point[3], xyzmax[0], xyzmax[1]); //xmax,ymax
    calProp_H2ONaCl(data->prop_cell, xyzc[0], xyzc[1]); //xc,yc
    data->phaseRegion_point[0] = regionIndex[0]; //phase index
    data->phaseRegion_point[1] = regionIndex[num_sample_x-1];
    data->phaseRegion_point[2] = regionIndex[num_sample_x*num_sample_y-num_sample_x];
    data->phaseRegion_point[3] = regionIndex[num_sample_x*num_sample_y-1];
    data->phaseRegion_cell     = phaseRegion_H2ONaCl_constantX(xyzc[0], xyzc[1]);
    // ========== 1. refinement check for Rho =====================
    // 如果梯度大于某个值则细化
    double threshold = 50;
    double mean_Rho = data->prop_cell.Rho;
    for(int i=0;i<4;i++)mean_Rho += data->prop_point[i].Rho;
    mean_Rho = mean_Rho / 5.0;
    double MSE_Rho = pow(data->prop_cell.Rho - mean_Rho, 2.0);
    for(int i=0;i<4;i++)MSE_Rho += pow(data->prop_point[i].Rho - mean_Rho, 2.0);
    MSE_Rho = MSE_Rho/5.0;
    
    // double dRho1 = fabs(data->prop_point[0].Rho - data->prop_point[1].Rho);
    // double dRho2 = fabs(data->prop_point[2].Rho - data->prop_point[3].Rho);
    // double dRho3 = fabs(data->prop_point[0].Rho - data->prop_point[2].Rho);
    // double dRho4 = fabs(data->prop_point[1].Rho - data->prop_point[3].Rho);
    // if(dRho1>threshold || dRho2>threshold || dRho3>threshold || dRho4>threshold) 
    // {
    //     data->need_refine_Rho = true;
        
    // }
    if(MSE_Rho > 10)
    {
        data->need_refine_Rho = true;
    }else{
        data->need_refine_Rho = false;
    }
    // ============================================================

    if(q->level > MAX_Level_Refine) return 0; //maximum-level control
    if(data->need_refine) return 1;
    if(data->need_refine_Rho) return 1;

    return 0;
}
static void fillPointData(p4est_iter_volume_info_t * info, void *user_data)
{
  sc_array_t         **user_data_ptr = (sc_array_t **)user_data;
  // point field
  sc_array_t          *phaseIndex = (sc_array_t *) (user_data_ptr[0]);
  sc_array_t          *Rho = (sc_array_t *) (user_data_ptr[1]);
  sc_array_t          *coord_vector = (sc_array_t *) (user_data_ptr[2]);
  // cell field
  sc_array_t          *coord_phaseIndex_cell = (sc_array_t *) (user_data_ptr[3]);
  sc_array_t          *coord_y_cell = (sc_array_t *) (user_data_ptr[4]);
  sc_array_t          *coord_vector_cell = (sc_array_t *) (user_data_ptr[5]);

  p4est_t            *p4est = info->p4est;
  p4est_quadrant_t   *q = info->quad;
  FIELD_DATA       *data = (FIELD_DATA *) q->p.user_data;
  p4est_locidx_t      arrayoffset, vectoroffsect;
  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid;  /* this is the index of q *within its tree's numbering*.  We want to convert it its index for all the quadrants on this process, which we do below */
  double              this_c;
  double              *this_c_ptr;
  double              *this_phaseIndex_ptr, *this_Rho_ptr, *this_coord_vector_ptr;
  double              *this_phaseIndex_cell_ptr, *this_y_cell_ptr, *this_coord_vector_cell_ptr;
  int                 i, j;
  double              xyz[3];
  double              midpoint[3];
  p4est_qcoord_t      length = P4EST_QUADRANT_LEN (q->level);
  p4est_qcoord_t      half_length = length/2;

  // cell data
  p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x + half_length, q->y + half_length, midpoint);
  this_phaseIndex_cell_ptr = (double *) sc_array_index (coord_phaseIndex_cell, local_id);
  this_y_cell_ptr = (double *) sc_array_index (coord_y_cell, local_id);
  this_coord_vector_cell_ptr = (double *) sc_array_index (coord_vector_cell, local_id*3);
  this_phaseIndex_cell_ptr[0] = data->phaseRegion_cell;//phaseRegion_H2ONaCl_constantX(midpoint[0], midpoint[1]); //midpoint[0];
  if(data->need_refine) this_phaseIndex_cell_ptr[0] = 999; //the smallest cell but still need refine, for this quadrant maybe just call function directally calculate everything rather than interpolate from lookuptable
  this_y_cell_ptr[0] = midpoint[1];
  this_coord_vector_cell_ptr[0] = midpoint[0];
  this_coord_vector_cell_ptr[1] = midpoint[1];
  this_coord_vector_cell_ptr[2] = midpoint[2];
  
  // point data
  arrayoffset = P4EST_CHILDREN * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in u_interp */
  vectoroffsect = arrayoffset*3;
  for (i = 0; i < P4EST_CHILDREN; i++) 
  {
    // this_c = data->c;
    p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x + ((i%2) ? length : 0), q->y + ((i/2) ? length : 0), xyz);
    // this_c = xyz[0]; // give it x coordinate
    this_phaseIndex_ptr = (double *) sc_array_index (phaseIndex, arrayoffset + i);
    this_Rho_ptr = (double *) sc_array_index (Rho, arrayoffset + i);
    this_coord_vector_ptr = (double *) sc_array_index (coord_vector, vectoroffsect + i*3);
    this_phaseIndex_ptr[0] = data->phaseRegion_point[i]; 
    this_Rho_ptr[0] = data->prop_point[i].Rho;//xyz[1];
    this_coord_vector_ptr[0] = xyz[0];
    this_coord_vector_ptr[1] = xyz[1];
    this_coord_vector_ptr[2] = xyz[2];
  }
}
static void write2vtk (p4est_t * p4est)
{
  char                filename[BUFSIZ] = "";
  int                 retval;
  p4est_locidx_t      numquads;
  // point field
  sc_array_t         *phaseIndex;  //scalar
  sc_array_t         *Rho;  //scalar
  sc_array_t         *coord_vector;  //vector field
  // cell field
  sc_array_t         *coord_phaseIndex_cell;  //scalar
  sc_array_t         *coord_y_cell;  //scalar
  sc_array_t         *coord_vector_cell;  //vector field
  p4est_vtk_context_t *context;

  snprintf (filename, BUFSIZ, "step6");
  numquads = p4est->local_num_quadrants; //cell number
  cout<<"Number of quads: "<<numquads<<endl;
  // have to use sc library even if you don't want to use MPI
  phaseIndex = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN); //create array of the point data: corner(vertex) number always equal to children number (dim-related)
  Rho = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN);
  coord_vector = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN * 3); //vector field must have 3 component, explained in p4est_vtk_write_cell_dataf
  coord_phaseIndex_cell = sc_array_new_size (sizeof (double), numquads); //create array of the point data: corner(vertex) number always equal to children number (dim-related)
  coord_y_cell = sc_array_new_size (sizeof (double), numquads);
  coord_vector_cell = sc_array_new_size (sizeof (double), numquads * 3);
  // package the pointer of field array to a pointer array and send the pointer array to p4est_iterate
  sc_array_t *fieldArray_ptr[]={phaseIndex, Rho, coord_vector, coord_phaseIndex_cell, coord_y_cell, coord_vector_cell}; //two scalar field and 1 vector field
  // fill data, option 1: use the iterator; option 2: loop every tree and quadrant of a tree
  p4est_iterate (p4est, NULL,   /* we don't need any ghost quadrants for this loop */
                 (void *) fieldArray_ptr,     /* pass in u_interp so that we can fill it */
                 fillPointData,    /* callback function that interpolates from the cell center to the cell corners, defined above */
                 NULL,          /* there is no callback for the faces between quadrants */
                 NULL);         /* there is no callback for the corners between quadrants */

  /* create VTK output context and set its parameters */
  context = p4est_vtk_context_new (p4est, filename);
  p4est_vtk_context_set_scale (context, 1);  /* quadrant at almost full scale */
//   p4est_geometry_t *geom;
//   p4est_vtk_context_set_geom (context, geom);
  /* begin writing the output files */
  context = p4est_vtk_write_header (context);
  SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing vtk header");
  context = p4est_vtk_write_cell_dataf (context, 1, 1,  /* do write the refinement level of each quadrant */
                                        0,      /* do write the mpi process id of each quadrant */
                                        0,      /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                        2,      /* there is no custom cell scalar data. */
                                        1,      /* there is no custom cell vector data. */
                                        "phaseIndex", coord_phaseIndex_cell,
                                        "midpoint_y", coord_y_cell,
                                        "coord_vector", coord_vector_cell,
                                        context);       /* mark the end of the variable cell data. */
  SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing cell data");

  /* write one scalar field: the solution value */
  context = p4est_vtk_write_point_dataf (context, 
        2, 1,
        "phaseIndex", phaseIndex, 
        "Rho",Rho, 
        "coord_vec", 
        coord_vector, 
        context);
  SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing point data");

  retval = p4est_vtk_write_footer (context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

  sc_array_destroy (phaseIndex); //remember destroy the new created array
  sc_array_destroy (Rho);
  sc_array_destroy (coord_vector);
  sc_array_destroy (coord_phaseIndex_cell); 
  sc_array_destroy (coord_y_cell);
  sc_array_destroy (coord_vector_cell);
}
int main()
{
    int level_init = 6; // will create a initial mesh with (2^dim)^(level_init) quadrants
    int  recursive;
    FIELD_DATA field_data;
    p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare ();
    p4est_t *p4est = p4est_new_ext (0, conn, 0, level_init, 1, sizeof(field_data), NULL, NULL);

    // refine 
    recursive = 1;
    p4est_refine (p4est, recursive, refine_fn, NULL);
    
    // write quadrant to vtu file
    // p4est_vtk_write_file (p4est, NULL, "step6");
    write2vtk(p4est); //Use customized output function
    // cout<<"global first qua level: "<<p4est->global_first_position->level<<endl;
    /* Destroy the p4est and the connectivity structure. */
    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);
    return 0;
}