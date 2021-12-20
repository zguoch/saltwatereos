#ifndef LOOKUPTABLEFOREST_H
#define LOOKUPTABLEFOREST_H
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
using namespace std;
#include <cmath>
// #include "H2ONaCl.H" 
#if USE_OMP == 1
    #include <omp.h>
#endif
namespace LOOKUPTABLE_FOREST
{
    #define MAX_FOREST_LEVEL 29

    struct Quad_index
    {
        int i, j, k;
        bool operator< (const Quad_index &ijk) const
        {
            return i < ijk.i || (i==ijk.i && j<ijk.j) || (i==ijk.i && j==ijk.j && k<ijk.k);
        }
    };
    
    template <int dim, typename USER_DATA>
    struct Quadrant
    {
        double              xyz[dim]; //real coordinate of the lower left corner of a quadrant
        Quad_index          ijk; //index of the LowerLeft corner in the reference space [2^MAX_FOREST_LEVEL, 2^MAX_FOREST_LEVEL, 2^MAX_FOREST_LEVEL]
        int       level;
        Quadrant*           parent;
        bool                isHasChildren;
        Quadrant            *children[1<<dim];// = NULL; //[1<<dim]; //2^dim: use dynamic array to save memory. Use dynamic array will cause bugs in read_forest, fix it later.
        USER_DATA           *user_data = NULL;
        // double              *pointData = NULL; //store all point data of a leaf quad, calculate it after refining process.
        // DEBUG
        // int index = -1;
    };
    
    /**
     * @brief For 2D case, define which variable is constant and the variable order of xy.
     * 
     */
    enum CONST_WHICH_VAR 
    {
        CONST_NO_VAR_TorHPX,     /**< No constant variable, 3D case. The x, y, z represents T/H, P, and X, respectively. T or H is specified by EOS_SPACE*/
        CONST_TorH_VAR_XP,      /**< Constant temperature T or specific enthalpy H, x represents salinity X and y represents pressure P. T or H is specified by EOS_SPACE */
        CONST_P_VAR_XTorH,      /**< Constant pressure P, x represents salinity X and y represents temperature T or specific enthalpy H. T or H is specified by EOS_SPACE */
        CONST_X_VAR_TorHP       /**< Constant salinity X, x represents temperature T or specific enthalpy H, and y represents pressure P. T or H is specified by EOS_SPACE */
    }; //only used for 2D case, CONST_NO means 3D
    
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

    // add data struct definition for different system, e.g., H2ONaCl. Actually we can move this data type definition to H2ONaCl.H, because LookUpTableForest class never care about this data type definite, it just accept whatever data type through template argument. But for the TCL API, if move this to other place, it will cause some compiling errors. So keep it here before finding better solution.
    
    /**
     * @brief Use which variable to express energy
     * 
     */
    enum EOS_ENERGY {
        EOS_ENERGY_T, /**< TPX space */
        EOS_ENERGY_H  /**< HPX space */
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
        long int m_num_quads;
        int m_num_leaves;
        size_t  m_data_size;
        double m_length_scale[dim]; /**< The reference space is a square or a cube with length=2^{MAX_FOREST_LEVEL}, so the length scale in x,y,z axis is calculated as, e.g. length_scale[0] = (m_xyz_max[0] - m_xyz_min[0])/length, so the length of a quadrant is len_quad = 2^{MAX_FOREST_LEVEL - level}, so its real length in x-axis is len_quad*length_scale[0] */
        Quadrant<dim,USER_DATA> m_root;
        void init_Root(Quadrant<dim,USER_DATA>& quad);
        void release_quadrant_data(Quadrant<dim,USER_DATA>* quad);
        void release_children(Quadrant<dim,USER_DATA>* quad);
        void getLeaves(vector<Quadrant<dim,USER_DATA>* >& leaves, long int& quad_counts, Quadrant<dim,USER_DATA>* quad);
        void refine(Quadrant<dim,USER_DATA>* quad, bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level));
        void release_map2data();
        void write_vtk_cellData(ofstream* fout, string type, string name, string format);
        void searchQuadrant(Quadrant<dim,USER_DATA>* quad_source, Quadrant<dim,USER_DATA> *&quad_target, double x_ref, double y_ref, double z_ref);
        void init(double xyz_min[dim], double xyz_max[dim], int max_level, size_t data_size, void* eosPointer);
        void write_forest(FILE* fpout, Quadrant<dim,USER_DATA>* quad, int order_child, bool is_write_data);
        void read_forest(FILE* fpin, Quadrant<dim,USER_DATA>* quad, int order_child);
        double* get_lowerleft_xyz(Quadrant<dim,USER_DATA>* quad);
        Quad_index get_lowerleft_ijk(Quadrant<dim,USER_DATA>* quad);
        void cal_xyz_quad(double* xyz_lower_left, int order_child, Quadrant<dim,USER_DATA>* quad);
        void cal_ijk_quad(Quad_index ijk_lower_left, int order_child, Quadrant<dim,USER_DATA>* quad);
        void construct_map2dat();
        void read_props_from_binary(string filename_forest);
        void read_forest_from_binary(string filename, bool read_only_header=false);
    public:
        void    *m_eosPointer;      //pass pointer of EOS object (e.g., the pointer of a object of cH2ONaCl class) to the forest through construct function, this will give access of EOS stuff in the refine call back function, e.g., calculate phase index and properties
        double  m_constZ;         // only valid when dim==2, i.e., 2D case, the constant value of third dimension, e.g. in T-P space with constant X.
        int     m_min_level;
        int     m_max_level;
        double m_xyz_min[dim];
        double m_xyz_max[dim];
        int     m_num_children;
        int     m_num_node_per_quad; //how many nodes will be stored in a quad: only for data storage
        int     m_num_props; //How many properties will be stored in the node
        std::map<int, propInfo> m_map_props;
        std::map<Quad_index, double*> m_map_ijk2data; //std::map<Quad_index, double> one pair of (i,j,k) to one array of data
        // int     m_index_TorH, m_index_P, m_index_X; //specify the index of variable T/H, P, X in the xyz array.
        CONST_WHICH_VAR m_const_which_var; 
        EOS_ENERGY m_TorH; 
        // double  m_physical_length_quad[MAX_FOREST_LEVEL][dim]; //Optimization: store the length of quad in each dimension as a member data of the forest, therefore don't need to calculate length of quad, just access this 2D array according to the quad level. 
        RMSD_RefineCriterion m_RMSD_RefineCriterion;
        inline void set_min_level(int min_level){m_min_level = min_level;};
        Quadrant<dim,USER_DATA>* get_root(){return &m_root;};
        // int searchQuadrant(double x, double y, double z);
        void searchQuadrant(Quadrant<dim,USER_DATA> *&targetLeaf, double x, double y, double z);
        void get_quadrant_physical_length(int level, double physical_length[dim]);
        void refine(bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level));
        void get_ijk_nodes_quadrant(Quadrant<dim,USER_DATA>* quad, int num_nodes_per_quad, Quad_index* ijk);
        void assemble_data(void (*cal_prop)(LookUpTableForest<dim,USER_DATA>* forest, std::map<Quad_index, double*>& map_ijk2data));
        void ijk2xyz(const Quad_index* ijk, double& x, double& y, double& z);
        void write_to_vtk_v1(string filename, bool write_data=true, bool isNormalizeXYZ=true);
        void write_to_vtk(string filename, bool write_data=true, bool isNormalizeXYZ=true);
        void write_to_binary(string filename, bool is_write_data=true);
        void print_summary();
        /**
         * @brief Construct a new Look Up Table Forest object. This is always used to create a 3D table
         * xyz would be corresponding to TPX or PHX. Note that the unit of T is K, unit of P is Pa, unit of X is wt% NaCl (e.g., seawater is 0.032), unit of H is J/kg. The same as H2ONaCl::cH2ONaCl::prop_pTX and The same as H2ONaCl::cH2ONaCl::prop_pHX.
         * For 3D case, the order of the variable MUST BE [T/H, p, X], T or H is specify by argument eos_space
         * 
         * @param xyz_min 
         * @param xyz_max 
         * @param max_level 
         * @param eosPointer 
         */
        LookUpTableForest(double xyz_min[dim], double xyz_max[dim], EOS_ENERGY TorH, int max_level, std::map<int, propInfo> name_props, void* eosPointer=NULL); //3D case
        /**
         * @brief Construct a new Look Up Table Forest object. This is always used to create a 2D table
         * xyz would be corresponding to TPX or PHX. Note that the unit of T is K, unit of P is Pa, unit of X is wt% NaCl (e.g., seawater is 0.032), unit of H is J/kg. 
         * The same as H2ONaCl::cH2ONaCl::prop_pTX and The same as H2ONaCl::cH2ONaCl::prop_pHX.
         * For 2D case, (1) if constant variable is X, xy order MUST BE T/H, P; (2) if constant variable is T/H, xy order MUST BE X, P; (3) if constant variable is P, xy order MUST BE X, T/H
         * @param xy_min 
         * @param xy_max 
         * @param constZ 
         * @param const_which_var 
         * @param max_level 
         * @param eosPointer 
         */
        LookUpTableForest(double xy_min[dim], double xy_max[dim], double constZ, CONST_WHICH_VAR const_which_var, EOS_ENERGY TorH, int max_level, std::map<int, propInfo> name_props, void* eosPointer=NULL); //2D case
        LookUpTableForest(string filename_forest, void* pointer=NULL); //load from exist binary file
        void destory();
        ~LookUpTableForest();
    };

//    #include "LookUpTableForestI.h" //have to organize the template class header file and sourfile together
    
}


#endif