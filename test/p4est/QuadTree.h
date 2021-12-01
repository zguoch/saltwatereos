#ifndef QUADTREE
#define QUADTREE
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
#include <cmath>

template <int dim, typename USER_DATA>
struct Quadrant
{
    int         xyz[3]; 
    int         level;
    bool        isHasChildren;
    Quadrant    *children[1<<dim]; //2^dim
    USER_DATA        *user_data;
};

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
    int     m_max_level;
    double m_xyz_min[3];
    double m_xyz_max[3];
    Quadrant<dim,USER_DATA> m_root;
    void init_Root(Quadrant<dim,USER_DATA>& quad);
    void release_quadrant_data(Quadrant<dim,USER_DATA>* quad);
    void write_to_vtk(const Quadrant<dim,USER_DATA>* quad);
    void getLeaves(vector<Quadrant<dim,USER_DATA>* >& leaves, Quadrant<dim,USER_DATA>* quad);
    void refine_uniform(Quadrant<dim,USER_DATA>* quad, int max_level);
public:
    void refine_uniform(int max_level);
    void write_to_vtk(string filename);
    cForest(double xyz_min[3], double xyz_max[3], int max_level, size_t data_size);
    ~cForest();
};


#endif