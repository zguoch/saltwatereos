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
#define MAX_Level_Refine 12
#define PI 3.141592653

struct FIELD_DATA
{
  bool need_refine; 
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
static int refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    FIELD_DATA       *data = (FIELD_DATA *) q->p.user_data;
    // quadrant的中心点坐标
    double c[3], min[3], max[3];
    double xmin, xmax, ymin, ymax, xc, yc;
    p4est_qcoord_t      half_length = P4EST_QUADRANT_LEN (q->level) / 2;
    // p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x + half_length, q->y + half_length, c);
    p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x + half_length*2, q->y + half_length*2, max);
    p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x, q->y, min);
    // 计算二次函数的值: ((x-0.5)/0.5)^2
    xmin = min[0], xmax = max[0];
    ymin = min[1], ymax = max[1];
    // 计算此quadrant的四个顶点的phase index， 如果全部相等且不为0（0表示在函数曲线的边界上，其实这个判断只是逻辑上需要）则需要细化
    int region1 = phaseRegion_H2ONaCl_constantX(xmin, ymin);
    int region2 = phaseRegion_H2ONaCl_constantX(xmax, ymin);
    int region3 = phaseRegion_H2ONaCl_constantX(xmin, ymax);
    int region4 = phaseRegion_H2ONaCl_constantX(xmax, ymax);
    if( !((region1 == region2) && (region1 == region3) && (region1 == region4)) )
    {
        data->need_refine = true;
    }else
    {
        data->need_refine = false;
    }
    
    if(q->level > MAX_Level_Refine) return 0;
    if(data->need_refine) return 1;

    return 0;
}
static void fillPointData(p4est_iter_volume_info_t * info, void *user_data)
{
  sc_array_t         **user_data_ptr = (sc_array_t **)user_data;
  // point field
  sc_array_t          *phaseIndex = (sc_array_t *) (user_data_ptr[0]);
  sc_array_t          *coord_y = (sc_array_t *) (user_data_ptr[1]);
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
  double              *this_phaseIndex_ptr, *this_y_ptr, *this_coord_vector_ptr;
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
  this_phaseIndex_cell_ptr[0] = phaseRegion_H2ONaCl_constantX(midpoint[0], midpoint[1]); //midpoint[0];
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
    this_y_ptr = (double *) sc_array_index (coord_y, arrayoffset + i);
    this_coord_vector_ptr = (double *) sc_array_index (coord_vector, vectoroffsect + i*3);
    this_phaseIndex_ptr[0] = phaseRegion_H2ONaCl_constantX(xyz[0], xyz[1]); //xyz[0]; //This also works: *this_c_ptr = xyz[0];
    this_y_ptr[0] = xyz[1];
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
  sc_array_t         *coord_y;  //scalar
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
  coord_y = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN);
  coord_vector = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN * 3); //vector field must have 3 component, explained in p4est_vtk_write_cell_dataf
  coord_phaseIndex_cell = sc_array_new_size (sizeof (double), numquads); //create array of the point data: corner(vertex) number always equal to children number (dim-related)
  coord_y_cell = sc_array_new_size (sizeof (double), numquads);
  coord_vector_cell = sc_array_new_size (sizeof (double), numquads * 3);
  // package the pointer of field array to a pointer array and send the pointer array to p4est_iterate
  sc_array_t *fieldArray_ptr[]={phaseIndex, coord_y, coord_vector, coord_phaseIndex_cell, coord_y_cell, coord_vector_cell}; //two scalar field and 1 vector field
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
        "y",coord_y, 
        "coord_vec", 
        coord_vector, 
        context);
  SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing point data");

  retval = p4est_vtk_write_footer (context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

  sc_array_destroy (phaseIndex); //remember destroy the new created array
  sc_array_destroy (coord_y);
  sc_array_destroy (coord_vector);
  sc_array_destroy (coord_phaseIndex_cell); 
  sc_array_destroy (coord_y_cell);
  sc_array_destroy (coord_vector_cell);
}
int main()
{
    int level_init = 3; // will create a initial mesh with (2^dim)^(level_init) quadrants
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