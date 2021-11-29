/**
 * @file step1.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief 在 x: [0,2pi], y: [-1,1]的平面上对二次曲线refine，曲线穿过的cell给1，没穿过的cell给0
 * @version 0.1
 * @date 2021-11-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <p4est_vtk.h> //用于输出vtk格式
#include <p4est_extended.h> //用于直接调用p4est_new_ext函数创建tree
#include <iostream>
using namespace std;
#define MAX_Level_Refine 9

double quadr_func(double x)
{
    return pow(x-0.5, 2.0)*4;
}
static int refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    if(q->level > MAX_Level_Refine)return 0;

    // quadrant的中心点坐标
    double c[3], min[3], max[3];
    double xmin, xmax, ymin, ymax, xc, yc;
    p4est_qcoord_t      half_length = P4EST_QUADRANT_LEN (q->level) / 2;
    p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x + half_length, q->y + half_length, c);
    p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x + half_length*2, q->y + half_length*2, max);
    p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x, q->y, min);
    // 计算二次函数的值: ((x-0.5)/0.5)^2
    xmin = min[0], xmax = max[0];
    ymin = min[1], ymax = max[1];
    xc = c[0]; yc = c[1];
    double value_xmin = quadr_func(xmin);
    double value_xmax = quadr_func(xmax);
    double value_xc = quadr_func(xc);
    double min_value = value_xmin, max_value = value_xmax;
    if(min_value > max_value){ min_value = value_xmax; max_value = value_xmin; }
    cout<<"xmin: "<<xmin<<", xmax: "<<xmax<<", ymin: "<<ymin<<", ymax: "<<ymax<<", minvalue: "<<min_value<<", maxValue: "<<max_value<<endl;
    // 如果函数点落到此quadrant里面则细化: 这是个比较困难的问题，如何判断函数曲线穿过该cell，需要寻找一个更聪明的办法解决之！
    if( ((value_xmin >= ymin) && (value_xmin <= ymax)) || ((value_xmax >= ymin) && (value_xmax <= ymax)) || ((value_xc >= ymin) && (value_xc <= ymax)))return 1;
    if(min_value<=ymin && max_value >= ymax ) return 1;

    return 0;
}
int main()
{
    int level_init = 1; // will create a initial mesh with (2^dim)^(level_init) quadrants
    int  recursive;
    p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare ();
    p4est_t *p4est = p4est_new_ext (0, conn, 0, level_init, 1, 0, NULL, NULL);

    // refine 
    recursive = 1;
    p4est_refine (p4est, recursive, refine_fn, NULL);

    p4est_vtk_write_file (p4est, NULL, "step4");
    
    /* Destroy the p4est and the connectivity structure. */
    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);
    return 0;
}