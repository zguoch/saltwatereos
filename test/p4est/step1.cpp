/**
 * @file step1.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Create a uniform mesh with 4-level tree. This is a good start.
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

int main()
{
    int level_init = 2; // will create a initial mesh with (2^dim)^(level_init) quadrants
    p4est_connectivity_t *conn = p4est_connectivity_new_unitsquare ();
    p4est_t *p4est = p4est_new_ext (0, conn, 0, level_init, 1, 0, NULL, NULL);

    p4est_vtk_write_file (p4est, NULL, "step1");
    
    /* Destroy the p4est and the connectivity structure. */
    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);
    return 0;
}