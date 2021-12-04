#ifndef QUADTREE
#define QUADTREE
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
#include <cmath>
#include "H2ONaCl.H" 

namespace LOOKUPTABLE_FOREST
{
    #define MAX_FOREST_LEVEL 29

    template <int dim, typename USER_DATA>
    struct Quadrant
    {
        double              xyz[dim]; //real coordinate of the lower left corner of a quadrant
        int                 level;
        bool                isHasChildren;
        Quadrant            *children[1<<dim]; //2^dim
        USER_DATA           *user_data;
        // DEBUG
        int index = -1;
    };

    enum NeedRefine {NeedRefine_NoNeed, NeedRefine_PhaseBoundary, NeedRefine_Rho, NeedRefine_H};
    /**
     * @brief Property refinement criterion, minimum RMSD of a quadran, if the RMSD of a property in a quadran grater than this criterion, it will be refined.
     * 
     */
    struct RMSD_RefineCriterion
    {
        double Rho;
        double H;
    };
    

    template <int dim>
    struct FIELD_DATA
    {
        // point data field
        H2ONaCl::PhaseRegion phaseRegion_point[1<<dim];
        H2ONaCl::PROP_H2ONaCl prop_point[1<<dim]; // properties at vertiex
        // cell data field
        NeedRefine need_refine; // indicator of what kind of the need-refined quad position (phase boundary), or what kind of properties need to refine
        H2ONaCl::PhaseRegion phaseRegion_cell;
        H2ONaCl::PROP_H2ONaCl prop_cell; // properties at midpoint as cell value (for vtk output)
    };

    inline int get_dim_from_binary(string filename)
    {
        FILE* fpin = NULL;
        fpin = fopen(filename.c_str(), "rb");
        if(!fpin)ERROR("Open file failed: "+filename);
        int dim0;
        fread(&dim0, sizeof(dim0), 1, fpin);
        fclose(fpin);
        return dim0;
    }
    /**
     * @brief Pass dimension and data type to the class
     * 
     * @tparam dim 
     * @tparam USER_DATA 
     */
    template <int dim, typename USER_DATA> 
    class LookUpTableForest
    {
    private:
        size_t  m_data_size;
        double m_xyz_min[dim];
        double m_xyz_max[dim];
        double m_length_scale[dim]; /**< The reference space is a square or a cube with length=2^{MAX_FOREST_LEVEL}, so the length scale in x,y,z axis is calculated as, e.g. length_scale[0] = (m_xyz_max[0] - m_xyz_min[0])/length, so the length of a quadrant is len_quad = 2^{MAX_FOREST_LEVEL - level}, so its real length in x-axis is len_quad*length_scale[0] */
        Quadrant<dim,USER_DATA> m_root;
        void init_Root(Quadrant<dim,USER_DATA>& quad);
        void release_quadrant_data(Quadrant<dim,USER_DATA>* quad);
        void release_children(Quadrant<dim,USER_DATA>* quad);
        void getLeaves(vector<Quadrant<dim,USER_DATA>* >& leaves, Quadrant<dim,USER_DATA>* quad);
        void refine(Quadrant<dim,USER_DATA>* quad, bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level));
        void write_vtk_cellData(ofstream* fout, string type, string name, string format);
        void searchQuadrant(Quadrant<dim,USER_DATA>* quad_source, Quadrant<dim,USER_DATA> *&quad_target, double x_ref, double y_ref, double z_ref);
        void init(double xyz_min[dim], double xyz_max[dim], int max_level, size_t data_size, void* eosPointer);
        void write_forest(FILE* fpout, Quadrant<dim,USER_DATA>* quad, bool is_write_data);
        void read_forest(FILE* fpin, Quadrant<dim,USER_DATA>* quad, bool is_read_data);
    public:
        void    *m_eosPointer;      //pass pointer of EOS object (e.g., the pointer of a object of cH2ONaCl class) to the forest through construct function, this will give access of EOS stuff in the refine call back function, e.g., calculate phase index and properties
        double  m_constZ;         // only valid when dim==2, i.e., 2D case, the constant value of third dimension, e.g. in T-P space with constant X.
        int     m_min_level;
        int     m_max_level;
        int     m_num_children;
        RMSD_RefineCriterion m_RMSD_RefineCriterion;
        inline void set_min_level(int min_level){m_min_level = min_level;};
        Quadrant<dim,USER_DATA>* get_root(){return &m_root;};
        int searchQuadrant(double x, double y, double z);
        void searchQuadrant(Quadrant<dim,USER_DATA> *&targetLeaf, double x, double y, double z);
        void get_quadrant_physical_length(int level, double physical_length[dim]);
        void refine(bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level));
        void write_to_vtk(string filename, bool write_data=true, bool isNormalizeXYZ=true);
        void write_to_binary(string filename, bool is_write_data=true);
        void read_from_binary(string filename, bool is_read_data=true);
        /**
         * @brief Construct a new Look Up Table Forest object. This is always used to create a 3D table
         * xyz would be corresponding to TPX or PHX. Note that the unit of T is K, unit of P is Pa, unit of X is wt% NaCl (e.g., seawater is 0.032), unit of H is J/kg. The same as H2ONaCl::cH2ONaCl::prop_pTX and The same as H2ONaCl::cH2ONaCl::prop_pHX
         * @param xyz_min 
         * @param xyz_max 
         * @param max_level 
         * @param data_size 
         * @param pointer 
         */
        LookUpTableForest(double xyz_min[dim], double xyz_max[dim], int max_level, void* eosPointer=NULL); //3D case
        /**
         * @brief Construct a new Look Up Table Forest object. This is always used to create a 2D table
         * xyz would be corresponding to TPX or PHX. Note that the unit of T is K, unit of P is Pa, unit of X is wt% NaCl (e.g., seawater is 0.032), unit of H is J/kg. The same as H2ONaCl::cH2ONaCl::prop_pTX and The same as H2ONaCl::cH2ONaCl::prop_pHX
         * @param xyz_min 
         * @param xyz_max 
         * @param constZ constant value of the third dimension, if want to create table in T-P space, constZ will be a constant value of salinity.
         * @param max_level 
         * @param data_size 
         * @param pointer 
         */
        LookUpTableForest(double xyz_min[dim], double xyz_max[dim], double constZ, int max_level, void* eosPointer=NULL); //2D case
        LookUpTableForest(string filename, void* pointer=NULL); //load from exist binary file
        void destory();
        ~LookUpTableForest();
    };

   #include "LookUpTableForestI.h" //have to organize the template class header file and sourfile together
    
}


#endif