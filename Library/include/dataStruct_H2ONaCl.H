/**
 * @file dataStruct_H2ONaCl.H
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Definition of data structure/type for H2ONaCl::cH2ONaCl
 * @version 1.0
 * @date 2020-11-16
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#ifndef DATASTRUCT_H2ONACL_H
#define DATASTRUCT_H2ONACL_H

#include "stdfunc.H"
#include "dataStruct.H"

namespace H2ONaCl
{
    /// Output file format
    enum fmtOutPutFile {
        /// <a href=https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf> VTK </a>
        fmt_vtk, 
        /// ASCII file with multi-columns
        fmt_dat
        };
    // enum CELL_TYPES_VTK{VTK_EMPTY, VTK_VERTEX, VTK_POLY_VERTEX, VTK_LINE, VTK_POLY_LINE, VTK_TRIANGLE, VTK_TRIANGLE_STRIP,
    //                     VTK_POLYGON, VTK_PIXEL, VTK_QUAD, VTK_TETRA, VTK_VOXEL, VTK_HEXAHEDRON, VTK_WEDGE, VTK_PYRAMID};
    // define index of region
    enum PhaseRegion {MixPhaseRegion=-1, SinglePhase_L, TwoPhase_L_V_X0, SinglePhase_V, 
                    TwoPhase_L_H, TwoPhase_V_H, ThreePhase_V_L_H, TwoPhase_V_L_L, TwoPhase_V_L_V, UnknownPhaseRegion}; //MixPhaseRegion is used in AMR lookup table
                    
    typedef std::map<int,std::string> MAP_PHASE_REGION;

    struct PROP_H2ONaCl
    {
        PhaseRegion Region; 
        double T, P, X_wt, H, Cp, Rho, Mu; /**< temperature(C), P(Pa), X_wt(wt%),bulk enthalpy [J/kg], bulk specific heat [J/kg/K], bulk density [kg/m3], viscosity [Pa s] */
        double Rho_l, Rho_v, Rho_h; //density [kg/m3]
        double H_l, H_v, H_h; //enthalpy [J/kg]
        double Cp_l, Cp_v, Cp_h;
        double S_l, S_v, S_h; //saturation [-]
        double X_l, X_v; // volume fraction of NaCl in vaper and liquid, it is a composition fraction. H2O + NaCl
        double Mu_l, Mu_v;//viscosity [Pa s]
        // derivatives
        double dRhodH;
    };
    
    struct MP_STRUCT
    {
        double b1, b1t, b1tt;
        double b2, b2t, b2tt;
    };
    struct ID_STRUCT
    {
        double f, ft, ftt;
    };
    struct TWOPHASEPROP_STRUCT
    {
        double f, p, s, g, u, h, dpd, dpt, cv, x;
    };
    struct BS_STRUCT
    {
        double x, f;
        double fd, fdd, ft, ftd, ftt;
    };
    struct RS_STRUCT
    {
        double f, ft, ftd;
        double fd, fdd, ftt;
    };
    struct Cr_STRUCT
    {
        double g[9][6], k[4], l[4], gg[4], t[4], d[4], a[4], b[4];
    };
    struct f_STRUCT
    {
        double f[11];
        double sum_f10;
    };
    struct TABLE4
    {
        double c[14], cA[11], d[11];
    };
}

#endif