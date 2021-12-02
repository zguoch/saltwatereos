#ifndef QUADTREE
#define QUADTREE
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
#include <cmath>
#include "H2ONaCl.H" 
H2ONaCl::cH2ONaCl eos;

#define MAX_FOREST_LEVEL 29

template <int dim, typename USER_DATA>
struct Quadrant
{
    double              xyz[3]; //real coordinate of the lower left corner of a quadrant
    int                 level;
    bool                isHasChildren;
    Quadrant            *children[1<<dim]; //2^dim
    USER_DATA           *user_data;
};

enum NeedRefine {NeedRefine_NoNeed, NeedRefine_PB_L2V, NeedRefine_PB_V2VL_L, NeedRefine_PB_V2VL_V, NeedRefine_PB_others, NeedRefine_Rho};
template <int dim>
struct FIELD_DATA
{
    // point data field
    int phaseRegion_point[1<<dim];
    H2ONaCl::PROP_H2ONaCl prop_point[1<<dim]; // properties at vertiex
    // cell data field
    NeedRefine need_refine; // indicator of what kind of the need-refined quad position (phase boundary), or what kind of properties need to refine
    int phaseRegion_cell;
    H2ONaCl::PROP_H2ONaCl prop_cell; // properties at midpoint as cell value (for vtk output)
};


// typedef template <int dim, typename USER_DATA> bool (*refine_func)(Quadrant<dim,USER_DATA>* quad);

/**
 * @brief Pass dimension and data type to the class
 * 
 * @tparam dim 
 * @tparam USER_DATA 
 */
template <int dim, typename USER_DATA> 
class cForest
{

private:
    int     m_num_children;
    size_t  m_data_size;
    double m_xyz_min[3];
    double m_xyz_max[3];
    double length_scale[3]; /**< The reference space is a square or a cube with length=2^{MAX_FOREST_LEVEL}, so the length scale in x,y,z axis is calculated as, e.g. length_scale[0] = (m_xyz_max[0] - m_xyz_min[0])/length, so the length of a quadrant is len_quad = 2^{MAX_FOREST_LEVEL - level}, so its real length in x-axis is len_quad*length_scale[0] */
    Quadrant<dim,USER_DATA> m_root;
    void init_Root(Quadrant<dim,USER_DATA>& quad);
    void release_quadrant_data(Quadrant<dim,USER_DATA>* quad);
    void release_children(Quadrant<dim,USER_DATA>* quad);
    void getLeaves(vector<Quadrant<dim,USER_DATA>* >& leaves, Quadrant<dim,USER_DATA>* quad);
    void refine(Quadrant<dim,USER_DATA>* quad, bool (*is_refine)(cForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level));
    
public:
    int     m_max_level;
    void get_quadrant_physical_length(int level, double physical_length[3]);
    void refine(bool (*is_refine)(cForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level));
    void write_to_vtk(string filename, bool write_data=true, bool isNormalizeXYZ=true);
    cForest(double xyz_min[3], double xyz_max[3], int max_level, size_t data_size);
    ~cForest();
};


#endif