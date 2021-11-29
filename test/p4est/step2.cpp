/**
 * @file step2.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief 获取每个quadrent的节点坐标
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
    double              midpoint[3];
    get_midpoint(p4est, which_tree, q, midpoint);
    double              coord[3];
    get_coord(p4est, which_tree, q, coord);
    // test output
    cout<<"q->level: "<<int(q->level)<<" which tree: "<<(int)(which_tree)<<endl;
    cout<<"  midpoint: "; for(int i=0;i<3;i++)cout<<midpoint[i]<<" ";
    cout<<"  coord: "; for(int i=0;i<3;i++)cout<<coord[i]<<" ";
    cout<<endl;
}
int main()
{
    int level_init = 2; // will create a initial mesh with (2^dim)^(level_init) quadrants
    p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare ();
    p4est_t *p4est = p4est_new_ext (0, conn, 0, level_init, 0, 0, process, NULL);

    // write quadrant to vtu file
    p4est_vtk_write_file (p4est, NULL, "step2");
    
    /* Destroy the p4est and the connectivity structure. */
    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);
    return 0;
}