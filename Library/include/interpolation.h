/**
 * @file interpolation.h
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Interpolation method definition
 * @version 0.1
 * @date 2021-12-05
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <cmath>
#ifndef INTERPOLATION
#define INTERPOLATION

// linear
//       ^
// f(x2) |..........*
//       |         /.
//       |        / .
// p=f(x)|.......?
//       |      /.  .
//       |     / .  .
//       |    /  .  .
// f(x1) |...*   .  .
//       |   .   .  .
//       -------------------->
//           x1  x  x2


// bilinear
//       ^
//       |    Q12     R2   Q22
//    y2 |...*.......o....*
//       |   .       .    .
//       |   .       .P   .
//     y |...........?.....
//       |   .       .    .
//       |   .       .    .
//       |   .Q11    .R1  .Q21
//    y1 |...*.......o....*
//       |   .       .    .
//       |   .       .    .
//       ---------------------->
//           x1      x    x2

namespace INTERPOLATION
{
    double linear(double x1, double f_x1, double x2, double f_x2, double x)
    {
        double result = (x - x1)/(x2-x1)*f_x2  + (x2-x)/(x2-x1)*f_x1;
        return result;
    }

    template<int dim>
    void bilinear(const double *xyz_min, const double *length, const double* values_at_vertices, const double* xyz, double& result)
    {
        double R1 = linear(xyz_min[0], values_at_vertices[0], xyz_min[0]+length[0], values_at_vertices[1], xyz[0]);
        double R2 = linear(xyz_min[0], values_at_vertices[2], xyz_min[0]+length[0], values_at_vertices[3], xyz[0]);
        result    = linear(xyz_min[1], R1, xyz_min[1] + length[1], R2, xyz[1]);
    }
}

#endif