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
#define PI 3.141592653

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
int phaseRegion_sin_func(double x_ref, double y_ref, double xmin = 0, double xmax = 2*PI, double ymin =-1, double ymax = 2)
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

static int refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{
    bool need_refine = false;
    p4est_qcoord_t      half_length = P4EST_QUADRANT_LEN (q->level) / 2;
    p4est_qcoord_t      length = P4EST_QUADRANT_LEN (q->level);
    const int num_sample_x =10; //最简单的情况就是只取xmin, xmax作为采样点判断这些采样点的函数计算返回值(flat)是否全部相等.但是有时候会有漏掉的情况，所以可以考虑在这里加密采样
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
            regionIndex[iy*num_sample_x + ix] = phaseRegion_sin_func(xyz_tmp[0], xyz_tmp[1]);
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
        need_refine = true;
    }else
    {
        need_refine = false;
    }

    if(q->level > MAX_Level_Refine) return 0; //maximum-level control
    if(need_refine) return 1;

    return 0;
}
int main()
{
    int level_init = 0; // will create a initial mesh with (2^dim)^(level_init) quadrants
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