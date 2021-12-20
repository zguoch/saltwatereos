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

/**
 * @brief Get the coeff bilinear object
 * 
 * @tparam dim 
 * @param xyz_min 
 * @param length 
 * @param xyz 
 * @param coeff Size of the coeff is [dim][2]
 */
    template<int dim>
    void get_coeff_bilinear(const double *xyz_min, const double *length, const double* xyz,double coeff[dim][2])
    {
        // double coeff[dim][2];
        coeff[0][0]     = (xyz[0]               - xyz_min[0])/length[0];
        coeff[0][1]     = (xyz_min[0] + length[0]   - xyz[0])/length[0];
        coeff[1][0]     = (xyz[1]               - xyz_min[1])/length[1];
        coeff[1][1]     = (xyz_min[1] + length[1]   - xyz[1])/length[1];
        if(dim==3)
        {
            coeff[2][0] = (xyz[2]               - xyz_min[2])/length[2];
            coeff[2][1] = (xyz_min[2] + length[2]   - xyz[2])/length[2];
        }
    }

    template<int dim>
    void bilinear_cal(double coeff[dim][2], const double* values_at_vertices, double& result)
    {
        double Y1_Z1    = coeff[0][0]*values_at_vertices[1]  + coeff[0][1]*values_at_vertices[0]; 
        double Y2_Z1    = coeff[0][0]*values_at_vertices[3]  + coeff[0][1]*values_at_vertices[2]; 
        double Z1       = coeff[1][0]*Y2_Z1  + coeff[1][1]*Y1_Z1; 
        if (dim==3)
        {
            double Y1_Z2 = coeff[0][0]*values_at_vertices[5]  + coeff[0][1]*values_at_vertices[4]; 
            double Y2_Z2 = coeff[0][0]*values_at_vertices[7]  + coeff[0][1]*values_at_vertices[6]; 
            double Z2    = coeff[1][0]*Y2_Z2  + coeff[1][1]*Y1_Z2; 
            result = coeff[2][0]*Z2  + coeff[2][1]*Z1; 
        }else
        {
            result = Z1;
        }
    }

    template<int dim>
    void bilinear(const double *xyz_min, const double *length, const double* values_at_vertices, const double* xyz, double& result)
    {
        double coeff[dim][2];
        get_coeff_bilinear<dim> (xyz_min, length, xyz, coeff);
        bilinear_cal<dim>(coeff, values_at_vertices, result);
    }
}

#endif