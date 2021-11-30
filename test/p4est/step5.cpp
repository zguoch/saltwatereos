/**
 * @file step5.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief 尝试对H2ONaCl的二位相图进行自适应网格细化。X=0.032, T in [1, 1000] C, P in [5, 500] bar
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
    if(q->level > MAX_Level_Refine)return 0;

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
        return 1;
    }
    // phaseRegion_H2ONaCl_constantX(xmin, ymin);

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

    p4est_vtk_write_file (p4est, NULL, "step5");
    
    /* Destroy the p4est and the connectivity structure. */
    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);
    return 0;
}