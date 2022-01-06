
#include "LookUpTableForest.h"

namespace LOOKUPTABLE_FOREST
{

    template <int dim, typename USER_DATA> 
    LookUpTableForest<dim,USER_DATA>::LookUpTableForest(double xyz_min[dim], double xyz_max[dim], EOS_ENERGY TorH, int max_level, std::map<int, propInfo> name_props, void* eosPointer)
    :
    m_constZ(xyz_min[dim-1]), //make this compatible with 2D case in the refine function.
    m_const_which_var(CONST_NO_VAR_TorHPX),
    m_TorH(TorH),
    m_map_props(name_props)
    {
        if(dim!=3)ERROR("This construct function only support dim=3, if you want do 2D table, please specify constZ and const_which_var! Note that there is no 1D support!");
        m_props_unique_points_leaves.num_props = name_props.size();
        init(xyz_min, xyz_max, max_level, sizeof(USER_DATA), eosPointer);
    }

    template <int dim, typename USER_DATA> 
    LookUpTableForest<dim,USER_DATA>::LookUpTableForest(double xyz_min[dim], double xyz_max[dim], double constZ, CONST_WHICH_VAR const_which_var, EOS_ENERGY TorH, int max_level, std::map<int, propInfo> name_props, void* eosPointer)
    :
    m_constZ(constZ),
    m_const_which_var(const_which_var),
    m_TorH(TorH),
    m_map_props(name_props)
    {
        if(dim!=2)ERROR("This construct function only support dim=2, if you want do 3D table, please get rid of constZ and const_which_var! Note that there is no 1D support!");
        m_props_unique_points_leaves.num_props = name_props.size();
        init(xyz_min, xyz_max, max_level, sizeof(USER_DATA), eosPointer);
    }

    template <int dim, typename USER_DATA> 
    LookUpTableForest<dim,USER_DATA>::LookUpTableForest(string filename_forest, void* eosPointer)
    {
        m_eosPointer = eosPointer;
        m_num_children = 1<<dim;
        m_num_node_per_quad = m_num_children; //use 4 nodes for 2d and 8 nodes for 3D at this moment, there is no necessary use more points!!! 
        m_data_size = sizeof(USER_DATA);
        init_Root(m_root);
        
        // read from binary file
        if(eosPointer==NULL) //if the eosPointer is NULL, only read header for lutInfo app
        {
            read_forest_from_binary(filename_forest, true);
        }else
        {
            // read forest
            bool isHasPointIndexFile = false;
            isHasPointIndexFile = read_forest_from_binary(filename_forest);
            // WAIT("after loading forest");
            if(isHasPointIndexFile) // if the point index exist, the point index will be filled in each leaf quad, then need to create properties dynamic array
            {
                m_props_unique_points_leaves.create();
            }else    // if the point index file: *.pi doesn't exist, construct point index for properties at unique point of leaves
            {
                construct_props_leaves(NULL);
            }
            
            // WAIT("after constructing map2data");
            // read data from binary
            read_props_from_binary(filename_forest); 
            // WAIT("after loading property data");
        }
        // print info
        print_summary();
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::init(double xyz_min[dim], double xyz_max[dim], int max_level, size_t data_size, void* eosPointer)
    {
        m_eosPointer = eosPointer;
        m_num_children = 1<<dim;
        m_num_node_per_quad = m_num_children; //use 4 nodes for 2d and 8 nodes for 3D at this moment
        m_data_size = data_size;
        m_min_level = 0;
        m_max_level = max_level;
        m_RMSD_RefineCriterion.Rho  = 0.01; // 1%
        m_RMSD_RefineCriterion.H    = 0.01; // 1%

        double length_forest = (1<<MAX_FOREST_LEVEL);
        for (size_t i = 0; i < dim; i++)
        {
            m_xyz_max[i] = xyz_max[i];
            m_xyz_min[i] = xyz_min[i];
            m_length_scale[i] = (m_xyz_max[i] - m_xyz_min[i])/length_forest;
        }
        // made a simple test, it doesn't speed up in this way.
        // init length of quad at every level
        // for (int i = 0; i < MAX_FOREST_LEVEL; i++)
        // {
        //     for (int j = 0; j < dim; j++)
        //     {
        //         m_physical_length_quad[i][j] = (1<<(MAX_FOREST_LEVEL-i))*m_length_scale[j];
        //     }
        // }

        init_Root(m_root);
    }

    template <int dim, typename USER_DATA> 
    LookUpTableForest<dim,USER_DATA>::~LookUpTableForest()
    {
        destory();
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::destory()
    {
        // cout<<"destroy forest"<<endl;
        release_quadrant_data(&m_root);
        while (m_root.isHasChildren)
        {
            release_leaves(&m_root); //release leaf quad every round
        }
        // release properties data
        m_props_unique_points_leaves.clear();
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::release_leaves(Quadrant<dim,USER_DATA>* quad)
    {
        if(quad==NULL)return; 
        static int quadnums = 0;
        if(quad->isHasChildren)
        {
            // bool release_all = false;
            for (int i = 0; i < m_num_children; i++)
            {
                if(quad->qData.nonleaf->children[i]->isHasChildren)
                {
                    release_leaves(quad->qData.nonleaf->children[i]);
                }else
                {
                    // delete child's qData.nonleaf first: here the child could also be a leaf quad, but we still can use qData.nonleaf to access its qData.leaf, because they (are union pointer) has the same pointer value.
                    delete quad->qData.nonleaf->children[i]->qData.nonleaf;
                    quad->qData.nonleaf->children[i]->qData.nonleaf = NULL;
                    delete quad->qData.nonleaf->children[i];
                    quad->qData.nonleaf->children[i] = NULL;
                    quad->isHasChildren = false;
                    // quadnums++;
                    // cout<<quad->level<<" "<<quadnums<<endl;
                }
                // release_all = (release_all || quad->qData.nonleaf->children[i]);
            }
            // if(!release_all)
            // {
            //     // cout<<"level: "<<quad->level<<endl;
            //     delete[] quad->qData.nonleaf->children;
            //     quad->qData.nonleaf->children = NULL;
            // }
        }else
        {
            return;
        }
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::release_quadrant_data(Quadrant<dim,USER_DATA>* quad)
    {
        if(quad==NULL)return;
        // if(quad->pointData)
        // {
        //     delete[] quad->pointData;
        //     quad->pointData = NULL;
        // }

        if(quad->isHasChildren)
        {
            for (int i = 0; i < m_num_children; i++)
            {
                release_quadrant_data(quad->qData.nonleaf->children[i]);
            }
        }else if (quad->qData.leaf->user_data)
        {
            delete quad->qData.leaf->user_data;
            quad->qData.leaf->user_data = NULL;
            // delete quad->qData.leaf; !!! DON'T need to delete leaf quad here, it will be deleted in the release_leaves function, because qData.leaf and qData.nonleaf can point to the same memory block.
            // That means the value of qData.leaf and qData.nonleaf always the same, because they are all pointer type, but they point to different type of data.
        }else
        {
            return;
        }
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::init_Root(Quadrant<dim,USER_DATA>& quad)
    {
        // quad.xyz[0] = m_xyz_min[0];
        // quad.xyz[1] = m_xyz_min[1];
        // dim == 3 ? quad.xyz[2] = m_xyz_min[2] : 0;
        // quad.ijk={0,0,0};
        quad.level  = 0;
        quad.isHasChildren = false;
        quad.qData.leaf = new LeafQuad<dim, USER_DATA>;
        // quad.qData.leaf->parent = NULL;
        
        // quad.children = NULL;
        if(m_data_size!=0) quad.qData.leaf->user_data = new USER_DATA;  //only allocate memory if it is a leaf, this will be released when a quadrent is refined.
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::write_forest(FILE* fpout_forest, FILE* fpout_point_index, Quadrant<dim,USER_DATA>* quad, int order_child, bool isWriteData)
    {
        fwrite(&quad->level, sizeof(unsigned char), 1, fpout_forest); //write level
        fwrite(&quad->isHasChildren, sizeof(bool), 1, fpout_forest); //write isHasChildren
        // cout<<"write, level: "<<quad->level<<", index: "<<quad->index<<endl;
        if(quad->isHasChildren)
        {
            for (int i = 0; i < m_num_children; i++)
            {
                write_forest(fpout_forest, fpout_point_index, quad->qData.nonleaf->children[i], i, isWriteData);
            }
        }else
        {
            // write cell data of leaf quad
            fwrite(quad->qData.leaf->user_data, sizeof(USER_DATA), 1, fpout_forest);
            fwrite(quad->qData.leaf->index_props, sizeof(int_pointIndex), m_num_children, fpout_point_index);
        }
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::write_point_index(FILE* fpout_point_index, Quadrant<dim,USER_DATA>* quad)
    {
        if(quad->isHasChildren)
        {
            for (int i = 0; i < m_num_children; i++)
            {
                write_point_index(fpout_point_index, quad->qData.nonleaf->children[i]);
            }
        }else
        {
            fwrite(quad->qData.leaf->index_props, sizeof(int_pointIndex), m_num_children, fpout_point_index);
        }
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::write_point_index(string filename_forest)
    {
        STATUS("Writting point index ... ");
        string filename_point_index = filename_forest + "." + ExtensionName_PointIndexFile;
        FILE* fpout_point_index = NULL;
        fpout_point_index = fopen(filename_point_index.c_str(), "wb");
        if(fpout_point_index == NULL)ERROR("Open file failed: "+filename_point_index);

        write_point_index(fpout_point_index, &m_root);

        fclose(fpout_point_index);
        STATUS("Write point index done: " + filename_forest + "."+ExtensionName_PointIndexFile);
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::write_to_binary(string filename, bool isWriteData)
    {
        STATUS("Write lookup table forest to binary file ...");
        int dim0 = dim;
        FILE* fpout_forest = NULL;
        FILE* fpout_point_index = NULL;
        fpout_forest = fopen(filename.c_str(), "wb");
        if(fpout_forest == NULL)ERROR("Open file failed: "+filename);
        string filename_point_index = filename + "." + ExtensionName_PointIndexFile;
        fpout_point_index = fopen(filename_point_index.c_str(), "wb");
        if(fpout_point_index == NULL)ERROR("Open file failed: "+filename_point_index);
        // write header
        fwrite(&dim0,       sizeof(int),    1,      fpout_forest);
        fwrite(&m_TorH,     sizeof(EOS_ENERGY),    1,      fpout_forest);
        fwrite(&m_const_which_var, sizeof(CONST_WHICH_VAR),    1,      fpout_forest);
        fwrite(m_xyz_min,   sizeof(double), dim,    fpout_forest);
        fwrite(m_xyz_max,   sizeof(double), dim,    fpout_forest);
        fwrite(&m_constZ,    sizeof(double), dim,    fpout_forest);
        fwrite(m_length_scale,  sizeof(double), dim, fpout_forest);
        fwrite(&m_min_level,    sizeof(int), 1, fpout_forest);
        fwrite(&m_max_level,    sizeof(int), 1, fpout_forest);
        fwrite(&m_num_node_per_quad, sizeof(int), 1, fpout_forest);
        // leaves and total quads info 
        fwrite(&m_num_quads, sizeof(long int), 1, fpout_forest);
        fwrite(&m_num_leaves, sizeof(int), 1, fpout_forest);
        fwrite(&m_props_unique_points_leaves.num_points, sizeof(m_props_unique_points_leaves.num_points), 1, fpout_forest);
        fwrite(&m_num_need_refine, sizeof(m_num_need_refine), 1, fpout_forest);
        // write props info
        int num_props = m_map_props.size();
        fwrite(&num_props, sizeof(int), 1, fpout_forest);
        for(auto &m : m_map_props)
        {
            fwrite(&m.first, sizeof(int), 1, fpout_forest);
            // short name
            fwrite(m.second.shortName, sizeof(char), STR_LENGTH_PROPINFO, fpout_forest);
            // long name
            fwrite(m.second.longName, sizeof(char), STR_LENGTH_PROPINFO, fpout_forest);
            // unit
            fwrite(m.second.unit, sizeof(char), STR_LENGTH_PROPINFO, fpout_forest);
        }
        fwrite(&m_RMSD_RefineCriterion, sizeof(RMSD_RefineCriterion), 1, fpout_forest);
        // recursion write forest and data
        write_forest(fpout_forest, fpout_point_index, &m_root, 0, isWriteData);
        // close file
        fclose(fpout_forest);
        fclose(fpout_point_index);
        STATUS("Writting lookup table forest to binary file done.");

        // write properties to binary file
        STATUS("Writting properties data to binary file ...");
        // \todo 应该给forst文件和所有匹配的prop数据文件头上写一个UUID用于匹配识别
        int ind_prop = 0;
        for(auto &map_props : m_map_props)
        {
            string filename_prop = filename+"."+map_props.second.shortName;
            STATUS_color(to_string(ind_prop)+" "+map_props.second.longName+": "+filename_prop, COLOR_BLUE);
            FILE* fpout_prop = NULL;
            fpout_prop = fopen(filename_prop.c_str(), "wb");
            if(fpout_prop == NULL)ERROR("Open file failed: "+filename_prop);
            for(int_pointIndex i = 0; i<m_props_unique_points_leaves.num_points; i++)
            {
                fwrite(&m_props_unique_points_leaves.data[i][ind_prop], sizeof(double), 1, fpout_prop);
            }
            fclose(fpout_prop);
            ind_prop++;
        }
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::ijk2xyz(const Quad_index* ijk, double& x, double& y, double& z)
    {
        x = m_xyz_min[0] + ijk->i*m_length_scale[0];
        y = m_xyz_min[1] + ijk->j*m_length_scale[1];
        if(dim==3) z = m_xyz_min[2] + ijk->k*m_length_scale[2];
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::read_forest(FILE* fpin_forest, FILE* fpin_point_index, Quadrant<dim,USER_DATA>* quad, int order_child)
    {
        fread(&quad->level, sizeof(quad->level), 1, fpin_forest); //write level
        fread(&quad->isHasChildren, sizeof(bool), 1, fpin_forest); //write isHasChildren
        // cout<<"read, level: "<<quad->level<<", has child: "<<quad->isHasChildren<<endl;
        
        if(quad->isHasChildren)
        {
            if(quad->qData.leaf->user_data)delete quad->qData.leaf->user_data;
            delete quad->qData.leaf;
            quad->qData.nonleaf = new NonLeafQuad<dim, USER_DATA>;
            for (int i = 0; i < m_num_children; i++)
            {
                quad->qData.nonleaf->children[i] = new Quadrant<dim,USER_DATA>;
                quad->qData.nonleaf->children[i]->qData.leaf = new LeafQuad<dim,USER_DATA>;
                quad->qData.nonleaf->children[i]->qData.leaf->parent = quad;
                read_forest(fpin_forest, fpin_point_index, quad->qData.nonleaf->children[i], i); 
            }
        }else
        {
            // load data
            quad->qData.leaf->user_data = new USER_DATA;
            fread(quad->qData.leaf->user_data, sizeof(USER_DATA), 1, fpin_forest);
            if(fpin_point_index)fread(quad->qData.leaf->index_props, sizeof(int_pointIndex), m_num_children, fpin_point_index);
        }
    }

    template <int dim, typename USER_DATA> 
    void LookUpTableForest<dim,USER_DATA>::read_props_from_binary(string filename_forest)
    {
        STATUS("Read lookup table properties from binary file ...");
        
        int ind_prop = 0;
        for(auto &map_props : m_map_props)
        {
            string filename_prop = filename_forest+"."+map_props.second.shortName;
            STATUS_color(to_string(ind_prop)+" "+map_props.second.longName+": "+filename_prop, COLOR_BLUE);
            FILE* fpin = NULL;
            fpin = fopen(filename_prop.c_str(), "rb");
            if(!fpin)
            {
                ERROR("Open file failed: "+filename_prop);
            }
            for(int_pointIndex i=0; i<m_props_unique_points_leaves.num_points; i++)
            {
                fread(&m_props_unique_points_leaves.data[i][ind_prop], sizeof(double), 1, fpin);
            }
            fclose(fpin);
            ind_prop++;
        }
    }

    template <int dim, typename USER_DATA> 
    bool LookUpTableForest<dim,USER_DATA>::read_forest_from_binary(string filename_forest, bool read_only_header)
    {
        STATUS("Check lookup table forest header ...");
        FILE* fpin = NULL;
        FILE* fpin_point_index = NULL;
        fpin = fopen(filename_forest.c_str(), "rb");
        if(!fpin)ERROR("Open file failed: "+filename_forest);
        // open point index file if it exist
        string filename_point_index = filename_forest + "." + ExtensionName_PointIndexFile;
        fpin_point_index = fopen(filename_point_index.c_str(), "rb");
        if(!fpin_point_index)WARNING("Open point index file failed, the point index will be reconstructed. "+filename_forest);

        int dim0;
        fread(&dim0, sizeof(dim0), 1, fpin);
        if(dim0 != dim)
        {
            cout<<"-- Dimension in the file is "<<dim0<<", but the temperate argument <dim> is "<<dim<<endl;
            ERROR("Dimension is not consistent, maybe change the template argument <dim>");
            fclose(fpin);
        }
        // read header
        fread(&m_TorH,              sizeof(EOS_ENERGY),         1,      fpin);
        fread(&m_const_which_var,   sizeof(CONST_WHICH_VAR),    1,      fpin);
        fread(m_xyz_min,            sizeof(double),             dim,    fpin);
        fread(m_xyz_max,            sizeof(double),             dim,    fpin);
        fread(&m_constZ,            sizeof(double),             dim,    fpin);
        fread(m_length_scale,       sizeof(double),             dim,    fpin);
        fread(&m_min_level,         sizeof(int),                1,      fpin);
        fread(&m_max_level,         sizeof(int),                1,      fpin);
        fread(&m_num_node_per_quad, sizeof(int),                1,      fpin);
        // leaves and total quads info
        fread(&m_num_quads, sizeof(long int), 1, fpin);
        fread(&m_num_leaves, sizeof(int), 1, fpin);
        fread(&m_props_unique_points_leaves.num_points, sizeof(m_props_unique_points_leaves.num_points), 1, fpin);
        fread(&m_num_need_refine, sizeof(m_num_need_refine), 1, fpin);
        // read props info
        int tmp_num_props = 0;
        fread(&tmp_num_props, sizeof(int), 1, fpin);
        for (int i = 0; i < tmp_num_props; i++)
        {
            int ind_prop;
            fread(&ind_prop, sizeof(int), 1, fpin);
            // short name
            fread(m_map_props[ind_prop].shortName, sizeof(char), STR_LENGTH_PROPINFO, fpin);
            // long name
            fread(m_map_props[ind_prop].longName, sizeof(char), STR_LENGTH_PROPINFO, fpin);
            // unit
            fread(m_map_props[ind_prop].unit, sizeof(char), STR_LENGTH_PROPINFO, fpin);
            // check if property file exist
            string filename_prop = filename_forest+"."+m_map_props[ind_prop].shortName;
            struct stat buffer;
            if(stat(filename_prop.c_str(), &buffer) != 0)
            {
                WARNING("Property: " + string(m_map_props[ind_prop].longName)+COLOR_RED+" (" + "Not Found" + COLOR_DEFAULT + "): "+filename_prop);
                m_map_props.erase(ind_prop);
            }
        }
        // only return existed properties
        m_props_unique_points_leaves.num_props = m_map_props.size();
        STATUS("Read lookup table forest from binary file ...");
        fread(&m_RMSD_RefineCriterion, sizeof(RMSD_RefineCriterion), 1, fpin);
        // recursion read forest and data
        if(!read_only_header)read_forest(fpin, fpin_point_index, &m_root, 0); //child order of the root it self is 0
        // close file
        fclose(fpin);
        if(fpin_point_index)fclose(fpin_point_index);
        STATUS("Reading lookup table forest done");
        return (bool)fpin_point_index;
    }

    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::get_ijk_nodes_quadrant(Quadrant<dim,USER_DATA>* quad, const Quad_index* ijk_quad, int num_nodes_per_quad, Quad_index* ijk)
    {
        int length_quad = 1<<(MAX_FOREST_LEVEL - quad->level);
        switch (num_nodes_per_quad)
        {
        case 1<<dim:
            {
                // initialize as lower left corner ijk
                for (int i = 0; i < 4; i++)
                {
                    ijk[i] = *ijk_quad;
                }
                ijk[1].i += length_quad;
                ijk[2].j += length_quad;
                ijk[3].i += length_quad; ijk[3].j += length_quad;
                if(dim==3)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        ijk[i+4] = ijk[i];
                        ijk[i+4].k += length_quad;
                    }
                }
            }
            break;
        default:
            // assert(num_nodes_per_quad!=1<<dim);
            ERROR("Number of nodes per quad only supports 2^dim so far");
            break;
        }
    }

    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::refine(Quadrant<dim,USER_DATA>* quad, double xmin_quad, double ymin_quad, double zmin_quad, bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, double xmin_quad, double ymin_quad, double zmin_quad, int max_level))
    {
        // WAIT("inside refine");
        // printf("refine on thread %d, min_level: %d, max_level: %d, level: %d\n", omp_get_thread_num(), m_min_level, m_max_level, quad->level);
        // 
        // NOTE that do not use for loop in the recursion, it will slow down the calculation
        bool do_refine = is_refine(this, quad, xmin_quad, ymin_quad, zmin_quad, m_max_level); //refine calculation function need xyz coordinate

        if(do_refine)
        {
            int length_child = 1<<(MAX_FOREST_LEVEL - quad->level -1); //Note that must -1, because this is the child length
            // if no children, create children
            if(!quad->isHasChildren)
            {
                // make children and release user_data of parent, because the non-leaf quad doesn't need user data
                // z index = 0

                // quad->qData.nonleaf->children = new Quadrant<dim,USER_DATA>*[1<<dim];
                // WAIT("Create children");
                // backup parent leaf quad data and create new nonleaf data because the quad has become a parent
                LeafQuad<dim, USER_DATA>* leaf_quad_backup = quad->qData.leaf;
                quad->qData.nonleaf = new NonLeafQuad<dim, USER_DATA>;
                // 1st child: lower left
                quad->qData.nonleaf->children[0] = new Quadrant<dim,USER_DATA>;
                quad->qData.nonleaf->children[0]->qData.leaf = new LeafQuad<dim, USER_DATA>;
                quad->qData.nonleaf->children[0]->level  = quad->level+1; //only calculate once
                quad->qData.nonleaf->children[0]->qData.leaf->parent = quad;
                quad->qData.nonleaf->children[0]->isHasChildren = false;
                // 2nd child: lower right
                quad->qData.nonleaf->children[1] = new Quadrant<dim,USER_DATA>;
                quad->qData.nonleaf->children[1]->qData.leaf = new LeafQuad<dim, USER_DATA>;
                quad->qData.nonleaf->children[1]->level  = quad->qData.nonleaf->children[0]->level;
                quad->qData.nonleaf->children[1]->qData.leaf->parent = quad;
                quad->qData.nonleaf->children[1]->isHasChildren = false;
                // 3th child: upper left
                quad->qData.nonleaf->children[2] = new Quadrant<dim,USER_DATA>;
                quad->qData.nonleaf->children[2]->qData.leaf = new LeafQuad<dim, USER_DATA>;
                quad->qData.nonleaf->children[2]->level  = quad->qData.nonleaf->children[0]->level;
                quad->qData.nonleaf->children[2]->qData.leaf->parent = quad;
                quad->qData.nonleaf->children[2]->isHasChildren = false;
                // 4th child: upper right
                quad->qData.nonleaf->children[3] = new Quadrant<dim,USER_DATA>;
                quad->qData.nonleaf->children[3]->qData.leaf = new LeafQuad<dim, USER_DATA>;
                quad->qData.nonleaf->children[3]->level  = quad->qData.nonleaf->children[0]->level;
                quad->qData.nonleaf->children[3]->qData.leaf->parent = quad;
                quad->qData.nonleaf->children[3]->isHasChildren = false;
                if(m_data_size!=0) // only check once
                {
                    quad->qData.nonleaf->children[0]->qData.leaf->user_data = new USER_DATA;
                    quad->qData.nonleaf->children[1]->qData.leaf->user_data = new USER_DATA;
                    quad->qData.nonleaf->children[2]->qData.leaf->user_data = new USER_DATA;
                    quad->qData.nonleaf->children[3]->qData.leaf->user_data = new USER_DATA;
                }
                // 3D case: z index =1
                if(dim == 3)
                {
                    // ------- the last four children ------
                    // 5th child: lower left
                    quad->qData.nonleaf->children[4] = new Quadrant<dim,USER_DATA>;
                    quad->qData.nonleaf->children[4]->qData.leaf = new LeafQuad<dim, USER_DATA>;
                    quad->qData.nonleaf->children[4]->level  = quad->qData.nonleaf->children[0]->level;
                    quad->qData.nonleaf->children[4]->qData.leaf->parent = quad;
                    quad->qData.nonleaf->children[4]->isHasChildren = false;
                    // 6th child: lower right
                    quad->qData.nonleaf->children[5] = new Quadrant<dim,USER_DATA>;
                    quad->qData.nonleaf->children[5]->qData.leaf = new LeafQuad<dim, USER_DATA>;
                    quad->qData.nonleaf->children[5]->level  = quad->qData.nonleaf->children[0]->level;
                    quad->qData.nonleaf->children[5]->qData.leaf->parent = quad;
                    quad->qData.nonleaf->children[5]->isHasChildren = false;
                    // 7th child: upper left
                    quad->qData.nonleaf->children[6] = new Quadrant<dim,USER_DATA>;
                    quad->qData.nonleaf->children[6]->qData.leaf = new LeafQuad<dim, USER_DATA>;
                    quad->qData.nonleaf->children[6]->level  = quad->qData.nonleaf->children[0]->level;
                    quad->qData.nonleaf->children[6]->qData.leaf->parent = quad;
                    quad->qData.nonleaf->children[6]->isHasChildren = false;
                    // 8th child: upper right
                    quad->qData.nonleaf->children[7] = new Quadrant<dim,USER_DATA>;
                    quad->qData.nonleaf->children[7]->qData.leaf = new LeafQuad<dim, USER_DATA>;
                    quad->qData.nonleaf->children[7]->level  = quad->qData.nonleaf->children[0]->level;
                    quad->qData.nonleaf->children[7]->qData.leaf->parent = quad;
                    quad->qData.nonleaf->children[7]->isHasChildren = false;
                    if(m_data_size!=0) // only check once
                    {
                        quad->qData.nonleaf->children[4]->qData.leaf->user_data = new USER_DATA; 
                        quad->qData.nonleaf->children[5]->qData.leaf->user_data = new USER_DATA; 
                        quad->qData.nonleaf->children[6]->qData.leaf->user_data = new USER_DATA; 
                        quad->qData.nonleaf->children[7]->qData.leaf->user_data = new USER_DATA; 
                    }
                }
                // delete parent data
                delete leaf_quad_backup->user_data;
                leaf_quad_backup->user_data = NULL;
                delete leaf_quad_backup; //delete leafquad data itself, because we don't need it any more.
                quad->isHasChildren = true;
                // WAIT("Create children done");
            }

            // for (int i = 0; i < m_num_children; i++)
            // {
            #if USE_OMP == 1
                #pragma omp task shared(is_refine) //firstprivate(quad, m_max_level) //
            #endif
                refine(quad->qData.nonleaf->children[0], xmin_quad, ymin_quad, zmin_quad, is_refine);

            #if USE_OMP == 1
                #pragma omp task shared(is_refine) //firstprivate(quad, m_max_level)
            #endif
                refine(quad->qData.nonleaf->children[1], xmin_quad + length_child*m_length_scale[0], ymin_quad, zmin_quad, is_refine);

            #if USE_OMP == 1
                #pragma omp task shared(is_refine) //firstprivate(quad, m_max_level)
            #endif
                refine(quad->qData.nonleaf->children[2], xmin_quad, ymin_quad + length_child*m_length_scale[1], zmin_quad, is_refine);

            #if USE_OMP == 1
                #pragma omp task shared(is_refine) //firstprivate(quad, m_max_level)
            #endif
                refine(quad->qData.nonleaf->children[3], xmin_quad + length_child*m_length_scale[0], ymin_quad + length_child*m_length_scale[1], zmin_quad, is_refine);

                if(dim==3)
                {
                    zmin_quad += length_child*m_length_scale[2];
                #if USE_OMP == 1
                    #pragma omp task shared(is_refine)
                #endif
                    refine(quad->qData.nonleaf->children[4], xmin_quad, ymin_quad, zmin_quad, is_refine);

                #if USE_OMP == 1
                    #pragma omp task shared(is_refine)
                #endif
                    refine(quad->qData.nonleaf->children[5], xmin_quad + length_child*m_length_scale[0], ymin_quad, zmin_quad, is_refine);

                #if USE_OMP == 1
                    #pragma omp task shared(is_refine)
                #endif
                    refine(quad->qData.nonleaf->children[6], xmin_quad, ymin_quad + length_child*m_length_scale[1], zmin_quad, is_refine);

                #if USE_OMP == 1
                    #pragma omp task shared(is_refine)
                #endif
                    refine(quad->qData.nonleaf->children[7], xmin_quad + length_child*m_length_scale[0], ymin_quad + length_child*m_length_scale[1], zmin_quad, is_refine);
                }
            #if USE_OMP == 1
                #pragma omp taskwait
            #endif
            // }
        }
    }

    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::refine(bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, double xmin_quad, double ymin_quad, double zmin_quad, int max_level))
    {
        // WAIT("Do refine");
        refine(&m_root, m_xyz_min[0], m_xyz_min[1], dim==3 ? m_xyz_min[2] : 0, is_refine);
        // WAIT("Refine done");
    }

    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::get_unique_points_leaves(std::map<Quad_index, int_pointIndex>& map_unique_points, int& num_leaves, long int& num_quads, int& num_need_refine, Quadrant<dim,USER_DATA>* quad, Quad_index ijk_quad, unsigned int length_quad)
    {
        num_quads++;
        if(quad->isHasChildren)
        {
            //update the length of children
            length_quad = length_quad>>1; //divided by 2 //1<<(MAX_FOREST_LEVEL - quad->level - 1);

            get_unique_points_leaves(map_unique_points, num_leaves, num_quads, num_need_refine, quad->qData.nonleaf->children[0], ijk_quad, length_quad);

            ijk_quad.i += length_quad;
            get_unique_points_leaves(map_unique_points, num_leaves, num_quads, num_need_refine, quad->qData.nonleaf->children[1], ijk_quad, length_quad);

            ijk_quad.i -= length_quad; ijk_quad.j += length_quad;
            get_unique_points_leaves(map_unique_points, num_leaves, num_quads, num_need_refine, quad->qData.nonleaf->children[2], ijk_quad, length_quad);

            ijk_quad.i += length_quad;
            get_unique_points_leaves(map_unique_points, num_leaves, num_quads, num_need_refine, quad->qData.nonleaf->children[3], ijk_quad, length_quad);
            if(dim==3)
            {
                ijk_quad.k += length_quad;

                ijk_quad.i -= length_quad; ijk_quad.j -= length_quad;
                get_unique_points_leaves(map_unique_points, num_leaves, num_quads, num_need_refine, quad->qData.nonleaf->children[4], ijk_quad, length_quad);

                ijk_quad.i += length_quad;
                get_unique_points_leaves(map_unique_points, num_leaves, num_quads, num_need_refine, quad->qData.nonleaf->children[5], ijk_quad, length_quad);

                ijk_quad.i -= length_quad; ijk_quad.j += length_quad;
                get_unique_points_leaves(map_unique_points, num_leaves, num_quads, num_need_refine, quad->qData.nonleaf->children[6], ijk_quad, length_quad);

                ijk_quad.i += length_quad;
                get_unique_points_leaves(map_unique_points, num_leaves, num_quads, num_need_refine, quad->qData.nonleaf->children[7], ijk_quad, length_quad);
            }
        }else
        {
            num_leaves++;
            Quad_index *ijk_nodes_quad = new Quad_index[m_num_node_per_quad];
            get_ijk_nodes_quadrant(quad, &ijk_quad, m_num_node_per_quad, ijk_nodes_quad);
            for (int i_node = 0; i_node < m_num_node_per_quad; i_node++)
            {
                if(map_unique_points.count(ijk_nodes_quad[i_node]) == 0)
                {
                    // cout<<ijk_nodes_quad[i_node].i<<" "<<ijk_nodes_quad[i_node].j<<" "<<ijk_nodes_quad[i_node].k<<endl;
                    map_unique_points[ijk_nodes_quad[i_node]] = 1;
                }
            }
            delete[] ijk_nodes_quad;
            // check need-refine quads
            if(quad->qData.leaf->user_data->need_refine)num_need_refine++;
        }
    }

    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::pass_props_pointer_leaves(std::map<Quad_index, int_pointIndex>& map_unique_points, Quadrant<dim,USER_DATA>* quad, Quad_index ijk_quad, unsigned int length_quad)
    {
        if(quad->isHasChildren)
        {
            //update the length of children
            length_quad = length_quad>>1; //divided by 2 //1<<(MAX_FOREST_LEVEL - quad->level - 1);

            pass_props_pointer_leaves(map_unique_points, quad->qData.nonleaf->children[0], ijk_quad, length_quad);

            ijk_quad.i += length_quad;
            pass_props_pointer_leaves(map_unique_points, quad->qData.nonleaf->children[1], ijk_quad, length_quad);

            ijk_quad.i -= length_quad; ijk_quad.j += length_quad;
            pass_props_pointer_leaves(map_unique_points, quad->qData.nonleaf->children[2], ijk_quad, length_quad);

            ijk_quad.i += length_quad;
            pass_props_pointer_leaves(map_unique_points, quad->qData.nonleaf->children[3], ijk_quad, length_quad);
            if(dim==3)
            {
                ijk_quad.k += length_quad;

                ijk_quad.i -= length_quad; ijk_quad.j -= length_quad;
                pass_props_pointer_leaves(map_unique_points, quad->qData.nonleaf->children[4], ijk_quad, length_quad);

                ijk_quad.i += length_quad;
                pass_props_pointer_leaves(map_unique_points, quad->qData.nonleaf->children[5], ijk_quad, length_quad);

                ijk_quad.i -= length_quad; ijk_quad.j += length_quad;
                pass_props_pointer_leaves(map_unique_points, quad->qData.nonleaf->children[6], ijk_quad, length_quad);

                ijk_quad.i += length_quad;
                pass_props_pointer_leaves(map_unique_points, quad->qData.nonleaf->children[7], ijk_quad, length_quad);
            }
        }else
        {
            Quad_index *ijk_nodes_quad = new Quad_index[m_num_node_per_quad];
            get_ijk_nodes_quadrant(quad, &ijk_quad, m_num_node_per_quad, ijk_nodes_quad);
            for (int i_node = 0; i_node < m_num_node_per_quad; i_node++)
            {
                // cout<<ijk_nodes_quad[i_node].i<<" "<<ijk_nodes_quad[i_node].j<<" "<<ijk_nodes_quad[i_node].k<<endl;
                quad->qData.leaf->index_props[i_node] = map_unique_points[ijk_nodes_quad[i_node]];
                // cout<<map_unique_points[ijk_nodes_quad[i_node]]<<endl;
            }
            delete[] ijk_nodes_quad;
        }
    }

    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::construct_props_leaves(void (*cal_prop)(LookUpTableForest<dim,USER_DATA>* forest, std::map<Quad_index, unsigned int>& map_ijk2data, double** data))
    {
        if(m_props_unique_points_leaves.num_props==0)return; //if there is not any property included, skip it!
        
        STATUS("Construct properties on unique points of leaves, it will take a while ...");
        // 0. clear old array of properties
        m_props_unique_points_leaves.clear();

        // 1. get point ijk map
        // WAIT("get point ijk map");
        std::map<Quad_index, int_pointIndex> map_unique_points;
        Quad_index ijk_quad;
        ijk_quad.i = 0;
        ijk_quad.j = 0;
        ijk_quad.k = 0;
        m_num_leaves = 0; m_num_quads = 0; m_num_need_refine = 0;
        get_unique_points_leaves(map_unique_points, m_num_leaves, m_num_quads, m_num_need_refine, &m_root, ijk_quad, 1<<(MAX_FOREST_LEVEL));

        // 2. create m_props_unique_points_leaves.data
        // WAIT("update map value");
        ASSERT(m_props_unique_points_leaves.data==NULL, "m_props_unique_points_leaves不为空，在new之前需要释放之前的内存");
        m_props_unique_points_leaves.num_points = map_unique_points.size();
        m_props_unique_points_leaves.data = new double*[m_props_unique_points_leaves.num_points];
        int_pointIndex ind = 0;
        for(auto &m : map_unique_points)
        {
            m.second = ind;
            m_props_unique_points_leaves.data[ind] = new double[m_props_unique_points_leaves.num_props];
            ind++;
        }

        // 3. pass index of property array pointer to each leaf quad
        // WAIT("pass index of property array pointer to each leaf quad");
        pass_props_pointer_leaves(map_unique_points, &m_root, ijk_quad, 1<<(MAX_FOREST_LEVEL));
        // cout<<"unique points: "<<map_unique_points.size()<<endl;

        // 4. calculate properties and fill to the data array
        // WAIT("calculate properties and fill to the data array");
        // if cal_prop is NULL, the construct_props_leaves will be used to construct property data and index for binary file reading
        if(cal_prop)cal_prop(this, map_unique_points, m_props_unique_points_leaves.data); 
        // for (int_pointIndex i = 0; i < m_props_unique_points_leaves.num_points; i++)
        // {
        //     cout<<m_props_unique_points_leaves.data[i][0]<<endl;
        // }
        // WAIT("Test output props value");
    }
    
    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::getLeaves(vector<Quadrant<dim,USER_DATA>* >& leaves, long int& quad_counts, Quadrant<dim,USER_DATA>* quad)
    {
        quad_counts++;

        if(quad->isHasChildren)
        {
            for (int i = 0; i < m_num_children; i++)
            {
                getLeaves(leaves, quad_counts, quad->qData.nonleaf->children[i]);
            }
        }else
        {
            // quad->index = leaves.size();
            leaves.push_back(quad);
        }
    }
    // ASCII version
    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::write_to_vtk(string filename, bool write_data, bool isNormalizeXYZ)
    {
        clock_t start = clock();
        STATUS("Write to vtu file starting ...");
        Quadrant<dim,USER_DATA> *targetLeaf = NULL;
        vector<Quadrant<dim,USER_DATA>* > leaves;
        long int quad_counts = 0;
        getLeaves(leaves, quad_counts, &m_root);
        // WAIT(to_string(leaves.size()));
        // get valid leaves, get rid of leaves with  H2ONaCl::UnknownPhaseRegion
        vector<size_t> index_valid_leaves;
        for (size_t i = 0; i < leaves.size(); i++)
        {
            if(leaves[i]->qData.leaf->user_data->phaseRegion_cell != H2ONaCl::UnknownPhaseRegion)
            {
                index_valid_leaves.push_back(i);
            }
        }
        int num_points_per_cell = 1<<dim;
        int num_cells = index_valid_leaves.size();
        int num_points = m_props_unique_points_leaves.num_points;
        int ind =0;
        double length_cell[dim]; //physical length
        double length_ref_cell = 0;
        Quad_index *ijk_nodes_quad = new Quad_index[m_num_node_per_quad]; // \todo 如果使用二阶插值，则需要更多节点，需要通过cellType进行判断：比如二维情况九点quad，那么需要限制max_level必须小于MAX_FOREST_LEVEL-2，不过这个好办，在构造函数里面判断一下进行安全检查就行

        ofstream fout(filename);
        int CELL_TYPE = dim==2 ? 8 : 11;
        // write head
        // WAIT("before writting vtu");
        // STATUS("write header");
        cout<<"    head ";
        fout<<"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<<endl;
        fout<<"  <UnstructuredGrid>"<<endl;
        fout<<"    <Piece NumberOfPoints=\""<<num_points<<"\" NumberOfCells=\""<<num_cells<<"\">"<<endl;
        // write point data 
        // STATUS("write point data");
        cout<<"pd ";
        fout<<"      <PointData>"<<endl;
        ind = 0;
        for(auto &m : m_map_props)
        {
            fout<<"        <DataArray type=\"Float32\" Name=\""<<m.second.longName<<"\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
            for (unsigned int i = 0; i < m_props_unique_points_leaves.num_points; i++)
            {
                fout<<" "<<m_props_unique_points_leaves.data[i][ind];
            }
            
            fout<<"\n        </DataArray>"<<endl;
            ind++;
        }
        fout<<"      </PointData>"<<endl;
        // write cell data 
        // STATUS("write cell data");
        cout<<"cd ";
        fout<<"      <CellData>"<<endl;
        // ---------- 1. phase index
        fout<<"        <DataArray type=\"Int32\" Name=\"phaseIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
        for (size_t i = 0; i < index_valid_leaves.size(); i++){ fout<<" "<<leaves[index_valid_leaves[i]]->qData.leaf->user_data->phaseRegion_cell;}
        fout<<"\n        </DataArray>"<<endl;
        // ---------- 2. need refine
        fout<<"        <DataArray type=\"Int32\" Name=\"needRefine\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
        for (size_t i = 0; i < index_valid_leaves.size(); i++){ fout<<" "<<leaves[index_valid_leaves[i]]->qData.leaf->user_data->need_refine;}
        fout<<"\n        </DataArray>"<<endl;
        fout<<"      </CellData>"<<endl;
        // write points
        // STATUS("write points (xyz)");
        cout<<"xyz ";
        fout<<"      <Points>"<<endl;
        fout<<"        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0.008838834896540746\" RangeMax=\"1.4053746938907843\">"<<endl;
        double scale=0; //scale = 0.001, keep some space between each quadrant
        double physical_length[dim]; //={m_xyz_max[0] - m_xyz_min[0], m_xyz_max[1] - m_xyz_min[1], m_xyz_max[2] - m_xyz_min[2]};
        for (int i = 0; i < dim; i++){ physical_length[i] = m_xyz_max[i] - m_xyz_min[i];}
        std::map<Quad_index, int_pointIndex> map_unique_points;
        Quad_index ijk_quad;
        ijk_quad.i = 0;
        ijk_quad.j = 0;
        ijk_quad.k = 0;
        m_num_leaves = 0; m_num_quads = 0; m_num_need_refine = 0;
        get_unique_points_leaves(map_unique_points, m_num_leaves, m_num_quads, m_num_need_refine, &m_root, ijk_quad, 1<<(MAX_FOREST_LEVEL));
        for(auto &ijk2data : map_unique_points)
        {
            fout<<"         "<<ijk2data.first.i<<" "<<ijk2data.first.j<<" "<<ijk2data.first.k<<endl;
        }
        fout<<"        </DataArray>"<<endl;
        fout<<"      </Points>"<<endl;
        // write cells
        // STATUS("write cells");
        cout<<"cells ";
        fout<<"      <Cells>"<<endl;
        fout<<"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<num_cells-1<<"\">"<<endl;
        for (int i = 0; i < num_cells; i++)
        {
            fout<<"         ";
            for (int i_node = 0; i_node < m_num_node_per_quad; i_node++)
            {
                fout<<leaves[index_valid_leaves[i]]->qData.leaf->index_props[i_node]<<" "; 
            }
            fout<<endl;
        }
        
        
        fout<<endl;
        fout<<"        </DataArray>"<<endl;
        fout<<"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"4\" RangeMax=\"64\">"<<endl;
        fout<<"        ";
        for (int i = 0; i < num_cells; i++)fout<<" "<<(i+1)*num_points_per_cell;
        fout<<endl;
        fout<<"        </DataArray>"<<endl;
        fout<<"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"8\" RangeMax=\"8\">"<<endl;
        fout<<"        ";
        for (int i = 0; i < num_cells; i++)fout<<" "<<CELL_TYPE; //dim dependent
        fout<<endl;
        fout<<"        </DataArray>"<<endl;
        fout<<"      </Cells>"<<endl;

        // write footer
        fout<<"    </Piece>"<<endl;
        fout<<"  </UnstructuredGrid>"<<endl;
        fout<<"</VTKFile>"<<endl;
        fout.close();
        cout<<"done"<<endl;

        delete[] ijk_nodes_quad;
        STATUS_time("Write to vtu file done: "+filename, clock()-start);
    }


    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::print_summary()
    {
        int factor = 1<<m_max_level;
        cout<<"======= Summary of the LookUp Table forest ======="<<endl;
        cout<<"Dimension: "<<dim<<" in ";
        if(m_TorH == EOS_ENERGY_T)cout<<"TPX space"<<endl;
        else if(m_TorH == EOS_ENERGY_H)cout<<"HPX space"<<endl;
        switch (m_const_which_var)
        {
        case CONST_P_VAR_XTorH:
            {
                cout<<"Constant P = "<<m_constZ/1E5<<" bar"<<endl;
                cout<<"X range ["<<m_xyz_min[0]*100<<", "<<m_xyz_max[0]*100<<"] wt.% NaCl, max resolution: "<<(m_xyz_max[0] - m_xyz_min[0])/factor*100<<" wt.% NaCl"<<endl;
                if(m_TorH == EOS_ENERGY_T)cout<<"T range ["<<m_xyz_min[1]-273.15<<", "<<m_xyz_max[1]-273.15<<"] deg.C, max resolution: "<<(m_xyz_max[1] - m_xyz_min[1])/factor<<" deg.C"<<endl;
                else if(m_TorH == EOS_ENERGY_H)cout<<"H range ["<<m_xyz_min[1]/1E6<<", "<<m_xyz_max[1]/1E6<<"] MJ/kg, max resolution: "<<(m_xyz_max[1] - m_xyz_min[1])/factor/1000<<" kJ/kg"<<endl;
            }
            break;
        case CONST_TorH_VAR_XP:
            {
                if(m_TorH == EOS_ENERGY_T)cout<<"Constant T = "<<m_constZ - 273.15<<" deg.C"<<endl;
                else if(m_TorH == EOS_ENERGY_H)cout<<"Constant H = "<<m_constZ/1E6<<" MJ/kg"<<endl;
                cout<<"X range ["<<m_xyz_min[0]*100<<", "<<m_xyz_max[0]*100<<"] wt.% NaCl, max resolution: "<<(m_xyz_max[0] - m_xyz_min[0])/factor*100<<" wt.% NaCl"<<endl;
                cout<<"P range ["<<m_xyz_min[1]/1E5<<", "<<m_xyz_max[1]/1E5<<"] bar, max resolution: "<<(m_xyz_max[1] - m_xyz_min[1])/factor/1E5<<" bar"<<endl;
            }
            break;
        case CONST_X_VAR_TorHP:
            {
                cout<<"Constant X = "<<m_constZ*100<<" wt.% NaCl"<<endl;
                if(m_TorH == EOS_ENERGY_T)cout<<"T range ["<<m_xyz_min[0]-273.15<<", "<<m_xyz_max[0]-273.15<<"] deg.C, max resolution: "<<(m_xyz_max[0] - m_xyz_min[0])/factor<<" deg.C"<<endl;
                else if(m_TorH == EOS_ENERGY_H)cout<<"H range ["<<m_xyz_min[0]/1E6<<", "<<m_xyz_max[0]/1E6<<"] MJ/kg, max resolution: "<<(m_xyz_max[0] - m_xyz_min[0])/factor/1000<<" kJ/kg"<<endl;
                cout<<"P range ["<<m_xyz_min[1]/1E5<<", "<<m_xyz_max[1]/1E5<<"] bar, max resolution: "<<(m_xyz_max[1] - m_xyz_min[1])/factor/1E5<<" bar"<<endl;
            }
            break;
        case CONST_NO_VAR_TorHPX:
            {
                if(m_TorH == EOS_ENERGY_T)cout<<"T range ["<<m_xyz_min[0]-273.15<<", "<<m_xyz_max[0]-273.15<<"] deg.C, max resolution: "<<(m_xyz_max[0] - m_xyz_min[0])/factor<<" deg.C"<<endl;
                else if(m_TorH == EOS_ENERGY_H)cout<<"H range ["<<m_xyz_min[0]/1E6<<", "<<m_xyz_max[0]/1E6<<"] MJ/kg, max resolution: "<<(m_xyz_max[0] - m_xyz_min[0])/factor/1000<<" kJ/kg"<<endl;
                cout<<"P range ["<<m_xyz_min[1]/1E5<<", "<<m_xyz_max[1]/1E5<<"] bar, max resolution: "<<(m_xyz_max[1] - m_xyz_min[1])/factor/1E5<<" bar"<<endl;
                cout<<"X range ["<<m_xyz_min[2]*100<<", "<<m_xyz_max[2]*100<<"] wt.% NaCl, max resolution: "<<(m_xyz_max[2] - m_xyz_min[2])/factor*100<<" wt.% NaCl"<<endl;
            }
            break;
        default:
            break;
        }
        cout<<"Min level: "<<m_min_level<<", max level: "<<m_max_level<<endl;
        cout<<"All "<<m_num_leaves<<" leaves, "<<m_num_quads-m_num_leaves<<" non-leaf quads. "<<m_props_unique_points_leaves.num_points<<" unique points on leaves."<<endl;
        cout<<"All "<<m_num_need_refine<<" leaves need refine. "
            <<"Next refinement will add "<<m_num_need_refine*(m_num_children - 1)<<" leaves."
            <<endl;
        // props info
        cout<<"Include "<<m_props_unique_points_leaves.num_props<<" properties on each node."<<endl;
        int ind =0;
        for(auto & m : m_map_props)
        {
            cout<<"  "<<ind<<": "<<COLOR_GREEN<<m.second.longName<<COLOR_DEFAULT<<": "<<m.second.shortName<<m.second.unit<<endl;
            ind++;
        }
        double byte_forest_leaves = sizeof(LeafQuad<dim, USER_DATA>) * m_num_leaves;
        double byte_forest_nonleaves = sizeof(NonLeafQuad<dim, USER_DATA>) * (m_num_quads - m_num_leaves);
        double byte_quads = sizeof(Quadrant<dim, USER_DATA>) * m_num_quads;
        double byte_per_property = sizeof(double)*m_props_unique_points_leaves.num_points;
        double byte_total = (byte_forest_leaves + byte_forest_nonleaves + byte_quads + byte_per_property*m_props_unique_points_leaves.num_props);
        cout<<"Memory estimate. Total: "
            <<byte2string(byte_total)<<"\n"
            <<"  Leaves: "<<byte2string(byte_forest_leaves)<<"; Nonleaves: "
            <<byte2string(byte_forest_nonleaves)<<"\n"
            <<"  Quads: "<<byte2string(byte_quads)<<"; Properties: "
            <<byte2string(byte_per_property)<<"/property."
            <<endl;
        cout<<"================== Summary end ==================="<<endl;
    }

    template <int dim, typename USER_DATA>
    string LookUpTableForest<dim,USER_DATA>::byte2string(double bytes)
    {
        std::stringstream ss;
        if(bytes<1024)
        {
            ss<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<bytes<<" bytes";
            return ss.str();
        }else if(bytes>=1024 && bytes<1024*1024)
        {
            ss<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<bytes/1024<<" kb";
            return ss.str();
        }else if(bytes>=1024*1024 && bytes<1024*1024*1024)
        {
            ss<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<bytes/1024/1024<<" Mb";
            return ss.str();
        }else if(bytes>=1024*1024*1024)
        {
            ss<<std::setiosflags(std::ios::fixed)<<std::setprecision(1)<<bytes/1024/1024/1024<<" Gb";
            return ss.str();
        }
        return " ";
    }

    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::get_quadrant_physical_length(int level, double physical_length[dim])
    {
        int length_ref = 1<<(MAX_FOREST_LEVEL - level);
        for (int i = 0; i < dim; i++)
        {
            physical_length[i] = m_length_scale[i] * length_ref;
        }
    }

    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::searchQuadrant(Quadrant<dim,USER_DATA> *quad_source, Quadrant<dim,USER_DATA> *&quad_target, double* xyz_min_target, double x_ref, double y_ref, double z_ref)
    {
        // NOTE that do not use for loop in the recursion, it will slow down the calculation
        if(!quad_source->isHasChildren)
        {
            quad_target = quad_source;
            return;
        }else
        {
            int length_child = 1<<(MAX_FOREST_LEVEL - quad_source->qData.nonleaf->children[0]->level); //any of a child to access the child level number, e.g., 0
            int childID_x = x_ref/length_child;
            int childID_y = y_ref/length_child;
            int childID_z = dim == 3 ? z_ref/length_child : 0;
            // some time if the input point at the upper right corner point, the childID_x==2, or childID_y==2, so need to make it to 1
            childID_x = childID_x >=2 ? 1 : childID_x;
            childID_y = childID_y >=2 ? 1 : childID_y;
            childID_z = childID_z >=2 ? 1 : childID_z;
            // cout<<"  childID_x: "<<childID_x<<", childID_y: "<<childID_y<<", level: "<<quad_source->level<<endl;
            if(childID_x==1)xyz_min_target[0] += length_child; //here only calculate the reference coordinate, multiply scale after searching
            if(childID_y==1)xyz_min_target[1] += length_child;
            if(dim==3)if(childID_z==1)xyz_min_target[2] += length_child;
            searchQuadrant(quad_source->qData.nonleaf->children[childID_z*4 + childID_y*2 + childID_x], quad_target, xyz_min_target,
                childID_x == 1 ? x_ref - length_child : x_ref, 
                childID_y == 1 ? y_ref - length_child : y_ref, 
                childID_z == 1 ? z_ref - length_child : z_ref
                );
        }
    }

    template <int dim, typename USER_DATA>
    void LookUpTableForest<dim,USER_DATA>::searchQuadrant(Quadrant<dim,USER_DATA> *&targetLeaf, double* xyz_min_target, double x, double y, double z)
    {
        // NOTE that do not use for loop in the recursion, it will slow down the calculation
        double x_ref = (x - m_xyz_min[0])/m_length_scale[0];
        double y_ref = (y - m_xyz_min[1])/m_length_scale[1];
        double z_ref = dim == 3 ? (z - m_xyz_min[2])/m_length_scale[2] : 0;
        // initialized by reference coordinate and calculate reference xyz as well inside the searchQuadrant
        xyz_min_target[0] = 0;
        xyz_min_target[1] = 0;
        if(dim==3)xyz_min_target[2] = 0;
        // WAIT("开始搜索");
        searchQuadrant(&m_root, targetLeaf, xyz_min_target, x_ref, y_ref, z_ref);
        // now calculate the physical xyz of the left corner of the target quad
        xyz_min_target[0] = xyz_min_target[0]*m_length_scale[0] + m_xyz_min[0];
        xyz_min_target[1] = xyz_min_target[1]*m_length_scale[1] + m_xyz_min[1];
        if(dim==3)xyz_min_target[2] = xyz_min_target[2]*m_length_scale[2] + m_xyz_min[2];
        // cout<<targetLeaf->level<<endl;
        // WAIT("搜索结束");
    }

    // template <int dim, typename USER_DATA>
    // int LookUpTableForest<dim,USER_DATA>::searchQuadrant(double x, double y, double z)
    // {
    //     double x_ref = (x - m_xyz_min[0])/m_length_scale[0];
    //     double y_ref = (y - m_xyz_min[1])/m_length_scale[1];
    //     double z_ref = dim == 3 ? (z - m_xyz_min[2])/m_length_scale[2] : 0;
    //     Quadrant<dim,USER_DATA>* targetLeaf;
    //     searchQuadrant(&m_root, targetLeaf, x_ref, y_ref, z_ref);
    //     return targetLeaf->index;
    // }
}