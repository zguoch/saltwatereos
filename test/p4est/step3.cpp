/**
 * @file step3.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief 将每个quadrant的中心点坐标作为cell data输出
 * @version 0.1
 * @date 2021-11-29
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <p4est_vtk.h> //用于输出vtk格式
#include <p4est_extended.h> //用于直接调用p4est_new_ext函数创建tree
#include <iostream>
using namespace std;

// data package
struct USER_DATA
{
  double              c; 
};
// get the midpoint of a quadrant
static void get_midpoint (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q, double xyz[3])
{
  p4est_qcoord_t      half_length = P4EST_QUADRANT_LEN (q->level) / 2;
  p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x + half_length, q->y + half_length, xyz);
}
// get coordinate of the lower-left corner of quadrent
static void get_coord (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q, double xyz[3])
{
  p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x, q->y, xyz); //must use p4est_qcoord_to_vertex transform coordinate to vertex space (range 0-1), otherwise q->x and q->y are very large value.
}
static void process (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    // double              midpoint[3];
    // get_midpoint(p4est, which_tree, q, midpoint);
    // double              coord[3];
    // get_coord(p4est, which_tree, q, coord);
    // // test output
    // cout<<"q->level: "<<int(q->level)<<" which tree: "<<(int)(which_tree)<<endl;
    // cout<<"  midpoint: "; for(int i=0;i<3;i++)cout<<midpoint[i]<<" ";
    // cout<<"  coord: "; for(int i=0;i<3;i++)cout<<coord[i]<<" ";
    // cout<<endl;
}
static void fillPointData(p4est_iter_volume_info_t * info, void *user_data)
{
  sc_array_t         **user_data_ptr = (sc_array_t **)user_data;
  // point field
  sc_array_t          *coord_x = (sc_array_t *) (user_data_ptr[0]);
  sc_array_t          *coord_y = (sc_array_t *) (user_data_ptr[1]);
  sc_array_t          *coord_vector = (sc_array_t *) (user_data_ptr[2]);
  // cell field
  sc_array_t          *coord_x_cell = (sc_array_t *) (user_data_ptr[3]);
  sc_array_t          *coord_y_cell = (sc_array_t *) (user_data_ptr[4]);
  sc_array_t          *coord_vector_cell = (sc_array_t *) (user_data_ptr[5]);

  p4est_t            *p4est = info->p4est;
  p4est_quadrant_t   *q = info->quad;
  USER_DATA       *data = (USER_DATA *) q->p.user_data;
  p4est_locidx_t      arrayoffset, vectoroffsect;
  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid;  /* this is the index of q *within its tree's numbering*.  We want to convert it its index for all the quadrants on this process, which we do below */
  double              this_c;
  double              *this_c_ptr;
  double              *this_x_ptr, *this_y_ptr, *this_coord_vector_ptr;
  double              *this_x_cell_ptr, *this_y_cell_ptr, *this_coord_vector_cell_ptr;
  int                 i, j;
  double              xyz[3];
  double              midpoint[3];
  p4est_qcoord_t      length = P4EST_QUADRANT_LEN (q->level);
  p4est_qcoord_t      half_length = length/2;

  // cell data
  p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x + half_length, q->y + half_length, midpoint);
  this_x_cell_ptr = (double *) sc_array_index (coord_x_cell, local_id);
  this_y_cell_ptr = (double *) sc_array_index (coord_y_cell, local_id);
  this_coord_vector_cell_ptr = (double *) sc_array_index (coord_vector_cell, local_id*3);
  this_x_cell_ptr[0] = midpoint[0];
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
    this_x_ptr = (double *) sc_array_index (coord_x, arrayoffset + i);
    this_y_ptr = (double *) sc_array_index (coord_y, arrayoffset + i);
    this_coord_vector_ptr = (double *) sc_array_index (coord_vector, vectoroffsect + i*3);
    this_x_ptr[0] = xyz[0]; //This also works: *this_c_ptr = xyz[0];
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
  sc_array_t         *coord_x;  //scalar
  sc_array_t         *coord_y;  //scalar
  sc_array_t         *coord_vector;  //vector field
  // cell field
  sc_array_t         *coord_x_cell;  //scalar
  sc_array_t         *coord_y_cell;  //scalar
  sc_array_t         *coord_vector_cell;  //vector field
  p4est_vtk_context_t *context;

  snprintf (filename, BUFSIZ, "step3");
  numquads = p4est->local_num_quadrants; //cell number
  cout<<"Number of quads: "<<numquads<<endl;
  // have to use sc library even if you don't want to use MPI
  coord_x = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN); //create array of the point data: corner(vertex) number always equal to children number (dim-related)
  coord_y = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN);
  coord_vector = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN * 3); //vector field must have 3 component, explained in p4est_vtk_write_cell_dataf
  coord_x_cell = sc_array_new_size (sizeof (double), numquads); //create array of the point data: corner(vertex) number always equal to children number (dim-related)
  coord_y_cell = sc_array_new_size (sizeof (double), numquads);
  coord_vector_cell = sc_array_new_size (sizeof (double), numquads * 3);
  // package the pointer of field array to a pointer array and send the pointer array to p4est_iterate
  sc_array_t *fieldArray_ptr[]={coord_x, coord_y, coord_vector, coord_x_cell, coord_y_cell, coord_vector_cell}; //two scalar field and 1 vector field
  // fill data, option 1: use the iterator; option 2: loop every tree and quadrant of a tree
  p4est_iterate (p4est, NULL,   /* we don't need any ghost quadrants for this loop */
                 (void *) fieldArray_ptr,     /* pass in u_interp so that we can fill it */
                 fillPointData,    /* callback function that interpolates from the cell center to the cell corners, defined above */
                 NULL,          /* there is no callback for the faces between quadrants */
                 NULL);         /* there is no callback for the corners between quadrants */

  /* create VTK output context and set its parameters */
  context = p4est_vtk_context_new (p4est, filename);
  p4est_vtk_context_set_scale (context, 0.99);  /* quadrant at almost full scale */
  /* begin writing the output files */
  context = p4est_vtk_write_header (context);
  SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing vtk header");
  context = p4est_vtk_write_cell_dataf (context, 1, 1,  /* do write the refinement level of each quadrant */
                                        0,      /* do write the mpi process id of each quadrant */
                                        0,      /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                        2,      /* there is no custom cell scalar data. */
                                        1,      /* there is no custom cell vector data. */
                                        "midpoint_x", coord_x_cell,
                                        "midpoint_y", coord_y_cell,
                                        "coord_vector", coord_vector_cell,
                                        context);       /* mark the end of the variable cell data. */
  SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing cell data");

  /* write one scalar field: the solution value */
  context = p4est_vtk_write_point_dataf (context, 2, 1,"x", coord_x, "y",coord_y, "coord_vec", coord_vector, context);
  SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing point data");

  retval = p4est_vtk_write_footer (context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

  sc_array_destroy (coord_x); //remember destroy the new created array
  sc_array_destroy (coord_y);
  sc_array_destroy (coord_vector);
  sc_array_destroy (coord_x_cell); 
  sc_array_destroy (coord_y_cell);
  sc_array_destroy (coord_vector_cell);
}
int main()
{
    int level_init = 4; // will create a initial mesh with (2^dim)^(level_init) quadrants
    USER_DATA user_data;
    p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare ();
    p4est_t *p4est = p4est_new_ext (0, conn, 0, level_init, 1, sizeof(user_data), process, (void *) (&user_data));

    // write quadrant to vtu file
    // p4est_vtk_write_file (p4est, NULL, "step3");
    write2vtk(p4est); //Use customized output function

    /* Destroy the p4est and the connectivity structure. */
    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);
    return 0;
}