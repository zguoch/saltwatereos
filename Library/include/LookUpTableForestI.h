
template <int dim, typename USER_DATA> 
LookUpTableForest<dim,USER_DATA>::LookUpTableForest(double xyz_min[dim], double xyz_max[dim], EOS_ENERGY TorH, int max_level, std::map<int, propInfo> name_props, void* eosPointer)
:
m_constZ(xyz_min[dim-1]), //make this compatible with 2D case in the refine function.
m_const_which_var(CONST_NO_VAR_TorHPX),
m_TorH(TorH),
m_num_props(name_props.size()),
m_map_props(name_props)
{
    if(dim!=3)ERROR("This construct function only support dim=3, if you want do 2D table, please specify constZ and const_which_var! Note that there is no 1D support!");
    init(xyz_min, xyz_max, max_level, sizeof(USER_DATA), eosPointer);
}

template <int dim, typename USER_DATA> 
LookUpTableForest<dim,USER_DATA>::LookUpTableForest(double xyz_min[dim], double xyz_max[dim], double constZ, CONST_WHICH_VAR const_which_var, EOS_ENERGY TorH, int max_level, std::map<int, propInfo> name_props, void* eosPointer)
:
m_constZ(constZ),
m_const_which_var(const_which_var),
m_TorH(TorH),
m_num_props(name_props.size()),
m_map_props(name_props)
{
    if(dim!=2)ERROR("This construct function only support dim=2, if you want do 3D table, please get rid of constZ and const_which_var! Note that there is no 1D support!");
    init(xyz_min, xyz_max, max_level, sizeof(USER_DATA), eosPointer);
}

template <int dim, typename USER_DATA> 
LookUpTableForest<dim,USER_DATA>::LookUpTableForest(string filename_forest, void* eosPointer)
{
    m_eosPointer = eosPointer;
    m_num_children = 1<<dim;
    m_num_node_per_quad = m_num_children; //use 4 nodes for 2d and 8 nodes for 3D at this moment
    m_data_size = sizeof(USER_DATA);

    // read from binary file
    if(eosPointer==NULL) //if the eosPointer is NULL, only read header for lutInfo app
    {
        read_forest_from_binary(filename_forest, true);
    }else
    {
        // read forest
        read_forest_from_binary(filename_forest);
        // construct ijk_map2data
        construct_map2dat();
        // read data from binary
        read_props_from_binary(filename_forest);
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
    m_RMSD_RefineCriterion.Rho  = 5;
    m_RMSD_RefineCriterion.H    = 1000;

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
    // release data first
    release_map2data();
    // cout<<"destroy forest"<<endl;
    release_quadrant_data(&m_root);
    release_children(&m_root);
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::release_children(Quadrant<dim,USER_DATA>* quad)
{
    if(quad->isHasChildren)
    {
        bool release_all = false;
        for (int i = 0; i < m_num_children; i++)
        {
            if(quad->children[i]->isHasChildren)
            {
                release_children(quad->children[i]);
            }else
            {
                delete quad->children[i];
                quad->children[i] = NULL;
                quad->isHasChildren = false;
            }
            // release_all = (release_all || quad->children[i]);
        }
        // if(!release_all)
        // {
        //     // cout<<"level: "<<quad->level<<endl;
        //     delete[] quad->children;
        //     quad->children = NULL;
        // }
    }else
    {
        return;
    }
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::release_quadrant_data(Quadrant<dim,USER_DATA>* quad)
{
    // if(quad->pointData)
    // {
    //     delete[] quad->pointData;
    //     quad->pointData = NULL;
    // }

    if(quad->isHasChildren)
    {
        for (int i = 0; i < m_num_children; i++)
        {
            release_quadrant_data(quad->children[i]);
        }
    }else if (quad->user_data)
    {
        delete quad->user_data;
        quad->user_data = NULL;
    }else
    {
        return;
    }
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::init_Root(Quadrant<dim,USER_DATA>& quad)
{
    quad.xyz[0] = m_xyz_min[0];
    quad.xyz[1] = m_xyz_min[1];
    dim == 3 ? quad.xyz[2] = m_xyz_min[2] : 0;
    quad.ijk={0,0,0};
    quad.level  = 0;
    quad.parent = NULL;
    quad.isHasChildren = false;
    // quad.children = NULL;
    if(m_data_size!=0) quad.user_data = new USER_DATA;  //only allocate memory if it is a leaf, this will be released when a quadrent is refined.
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::write_forest(FILE* fpout, Quadrant<dim,USER_DATA>* quad, int order_child, bool isWriteData)
{
    // fwrite(quad, sizeof(Quadrant<dim,USER_DATA>), 1, fpout); //write all quads
    fwrite(&quad->level, sizeof(int), 1, fpout); //write level
    fwrite(&quad->isHasChildren, sizeof(bool), 1, fpout); //write isHasChildren
    // cout<<"write, level: "<<quad->level<<", index: "<<quad->index<<endl;
    if(quad->isHasChildren)
    {
        for (int i = 0; i < m_num_children; i++)
        {
            write_forest(fpout, quad->children[i], i, isWriteData);
        }
    }else
    {
        // write xyz at the lower left corner
        if(order_child==0)
        {
            fwrite(quad->xyz, sizeof(double), dim, fpout); //only write xyz of eldest brother
            fwrite(&quad->ijk, sizeof(Quad_index), 1, fpout); //only write ijk of eldest brother
        }
        // write data
        fwrite(quad->user_data, sizeof(USER_DATA), 1, fpout);
        // fwrite(quad->pointData, sizeof(double), m_num_children*m_num_props, fpout);
    }
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::write_to_binary(string filename, bool isWriteData)
{
    STATUS("Write lookup table forest to binary file ...");
    int dim0 = dim;
    FILE* fpout = NULL;
    fpout = fopen(filename.c_str(), "wb");
    if(fpout == NULL)ERROR("Open file failed: "+filename);
    // write header
    fwrite(&dim0,       sizeof(int),    1,      fpout);
    fwrite(&m_TorH,     sizeof(EOS_ENERGY),    1,      fpout);
    fwrite(&m_const_which_var, sizeof(CONST_WHICH_VAR),    1,      fpout);
    fwrite(m_xyz_min,   sizeof(double), dim,    fpout);
    fwrite(m_xyz_max,   sizeof(double), dim,    fpout);
    fwrite(&m_constZ,    sizeof(double), dim,    fpout);
    fwrite(m_length_scale,  sizeof(double), dim, fpout);
    fwrite(&m_min_level,    sizeof(int), 1, fpout);
    fwrite(&m_max_level,    sizeof(int), 1, fpout);
    fwrite(&m_num_node_per_quad, sizeof(int), 1, fpout);
    // leaves and total quads info 
    // vector<Quadrant<dim,USER_DATA>* > leaves;
    // long int quad_counts = 0;
    // getLeaves(leaves,quad_counts, &m_root);
    // int num_leaves = leaves.size();
    fwrite(&m_num_quads, sizeof(long int), 1, fpout);
    fwrite(&m_num_leaves, sizeof(int), 1, fpout);
    // write props info
    int num_props = m_map_props.size();
    fwrite(&num_props, sizeof(int), 1, fpout);
    for(auto &m : m_map_props)
    {
        int len_str = 0;
        fwrite(&m.first, sizeof(int), 1, fpout);
        // short name
        len_str = m.second.shortName.length();
        fwrite(&len_str, sizeof(int), 1, fpout);
        fwrite(m.second.shortName.c_str(), sizeof(char), len_str, fpout);
        // long name
        len_str = m.second.longName.length();
        fwrite(&len_str, sizeof(int), 1, fpout);
        fwrite(m.second.longName.c_str(), sizeof(char), len_str, fpout);
        // unit
        len_str = m.second.unit.length();
        fwrite(&len_str, sizeof(int), 1, fpout);
        fwrite(m.second.unit.c_str(), sizeof(char), len_str, fpout);
        // const char* name = m.second.c_str();
        // cout<<"length of name: "<<m.second.length()<<" "<<m.second<<" "<<sizeof(m.second.c_str())<<endl;
    }
    fwrite(&m_RMSD_RefineCriterion, sizeof(RMSD_RefineCriterion), 1, fpout);
    // recursion write forest and data
    write_forest(fpout, &m_root, 0, isWriteData);
    // close file
    fclose(fpout);
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
        for(auto &map_ijk2data : m_map_ijk2data)
        {
            fwrite(&map_ijk2data.second[ind_prop], sizeof(double), 1, fpout_prop);
        }
        fclose(fpout_prop);
        ind_prop++;
    }
}

/**
 * @brief Get xyz of a quad according to the xyz of his eldest brother. 
 * Call this function when reading the bin file.
 * 
 * @tparam dim 
 * @tparam USER_DATA 
 * @param eldest_brother 
 * @param order_youngerBrother 
 * @param xyz 
 */
template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::cal_xyz_quad(double* xyz_lower_left, int order_youngerBrother, Quadrant<dim,USER_DATA>* quad)
{
    int length_child = 1<<(MAX_FOREST_LEVEL - quad->level);
    // cout<<((order_youngerBrother - order_youngerBrother/4 * 4)/2)<<", "<<order_youngerBrother<<endl;
    // cout<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<endl;
    quad->xyz[0] = xyz_lower_left[0] + (order_youngerBrother%2)*length_child*m_length_scale[0];
    quad->xyz[1] = xyz_lower_left[1] + ((order_youngerBrother - order_youngerBrother/4 * 4)/2) *length_child*m_length_scale[1];
    if(dim==3)
    quad->xyz[2] = xyz_lower_left[2] + (order_youngerBrother/4)*length_child*m_length_scale[2];
    // cout<<"  "<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<endl;
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::cal_ijk_quad(Quad_index ijk_lower_left, int order_youngerBrother, Quadrant<dim,USER_DATA>* quad)
{
    int length_child = 1<<(MAX_FOREST_LEVEL - quad->level);
    // cout<<((order_youngerBrother - order_youngerBrother/4 * 4)/2)<<", "<<order_youngerBrother<<endl;
    // cout<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<endl;
    quad->ijk.i = ijk_lower_left.i + (order_youngerBrother%2)*length_child;
    quad->ijk.j = ijk_lower_left.j + ((order_youngerBrother - order_youngerBrother/4 * 4)/2) *length_child;
    if(dim==3)
    quad->ijk.k = ijk_lower_left.k + (order_youngerBrother/4)*length_child;
    // cout<<"  "<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<endl;
}

template <int dim, typename USER_DATA> 
double* LookUpTableForest<dim,USER_DATA>::get_lowerleft_xyz(Quadrant<dim,USER_DATA>* quad)
{
    if(quad->isHasChildren)
    {
        return get_lowerleft_xyz(quad->children[0]);
    }else
    {
        return quad->xyz;
    }
}

template <int dim, typename USER_DATA> 
Quad_index LookUpTableForest<dim,USER_DATA>::get_lowerleft_ijk(Quadrant<dim,USER_DATA>* quad)
{
    if(quad->isHasChildren)
    {
        return get_lowerleft_ijk(quad->children[0]);
    }else
    {
        return quad->ijk;
    }
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::ijk2xyz(const Quad_index* ijk, double& x, double& y, double& z)
{
    x = m_root.xyz[0] + ijk->i*m_length_scale[0];
    y = m_root.xyz[1] + ijk->j*m_length_scale[1];
    if(dim==3) z = m_root.xyz[2] + ijk->k*m_length_scale[2];
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::read_forest(FILE* fpin, Quadrant<dim,USER_DATA>* quad, int order_child)
{
    // fread(quad, sizeof(Quadrant<dim,USER_DATA>), 1, fpin); //write all quads
    fread(&quad->level, sizeof(int), 1, fpin); //write level
    fread(&quad->isHasChildren, sizeof(bool), 1, fpin); //write isHasChildren
    // cout<<"read, level: "<<quad->level<<", index: "<<quad->index<<", child: "<<quad->isHasChildren<<endl;
    
    if(quad->isHasChildren)
    {
        // quad->children = new Quadrant<dim,USER_DATA>*[1<dim];
        for (int i = 0; i < m_num_children; i++)
        {
            quad->children[i] = new Quadrant<dim,USER_DATA>;
            quad->children[i]->parent = quad;
            read_forest(fpin, quad->children[i], i); 
        }
    }else
    {
        if(order_child==0) //only read xyz of eldest brother
        {
            fread(quad->xyz, sizeof(double), dim, fpin); 
            fread(&quad->ijk, sizeof(Quad_index), 1, fpin); //only read ijk of eldest brother
        }else
        {
            // get left corner xyz, find it's elderest brother->first child->first child ...
            // 找大哥要xyz，如果大哥有儿子那就找大哥大儿子要xyz；如果大哥的大儿子还有儿子那就找大哥的大孙子要xyz... 直到最后一代一定能找到左下角的xyz
            // get_lowerleft_xyz第一个参数就是大哥的内存地址
            cal_xyz_quad(get_lowerleft_xyz(quad->parent->children[0]), order_child, quad);
            cal_ijk_quad(get_lowerleft_ijk(quad->parent->children[0]), order_child, quad);
        }
        // load data
        quad->user_data = new USER_DATA;
        fread(quad->user_data, sizeof(USER_DATA), 1, fpin);
        // quad->pointData = new double[m_num_props*m_num_node_per_quad];
        // fread(quad->pointData, sizeof(double), m_num_props*m_num_node_per_quad, fpin);
        // return;
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
        for(auto &map_ijk2data : m_map_ijk2data)
        {
            fread(&map_ijk2data.second[ind_prop], sizeof(double), 1, fpin);
        }
        fclose(fpin);
        ind_prop++;
    }
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::read_forest_from_binary(string filename, bool read_only_header)
{
    STATUS("Read lookup table forest from binary file ...");
    FILE* fpin = NULL;
    fpin = fopen(filename.c_str(), "rb");
    if(!fpin)ERROR("Open file failed: "+filename);
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
    // read props info
    fread(&m_num_props, sizeof(int), 1, fpin);
    for (int i = 0; i < m_num_props; i++)
    {
        int ind_prop;
        fread(&ind_prop, sizeof(int), 1, fpin);
        int len_str;
        // short name
        fread(&len_str, sizeof(int), 1, fpin);
        char* str_shortName = new char[len_str];
        fread(str_shortName, sizeof(char), len_str, fpin);
        m_map_props[ind_prop].shortName = string(str_shortName);
        // long name
        fread(&len_str, sizeof(int), 1, fpin);
        char* str_longName = new char[len_str];
        fread(str_longName, sizeof(char), len_str, fpin);
        m_map_props[ind_prop].longName = string(str_longName);
        // unit
        fread(&len_str, sizeof(int), 1, fpin);
        char* str_unit = new char[len_str];
        fread(str_unit, sizeof(char), len_str, fpin);
        m_map_props[ind_prop].unit = string(str_unit);
        // cout<<"len: "<<len_name<<", name: "<<m_map_props[ind_prop].shortName<<endl;
        delete[] str_shortName;
        delete[] str_longName;
        delete[] str_unit;
    }
    fread(&m_RMSD_RefineCriterion, sizeof(RMSD_RefineCriterion), 1, fpin);
    // recursion read forest and data
    if(!read_only_header)read_forest(fpin, &m_root, 0); //child order of the root it self is 0
    // close file
    fclose(fpin);
    STATUS("Reading lookup table forest done");
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::get_ijk_nodes_quadrant(Quadrant<dim,USER_DATA>* quad, int num_nodes_per_quad, Quad_index* ijk)
{
    int length_quad = 1<<(MAX_FOREST_LEVEL - quad->level);
    switch (num_nodes_per_quad)
    {
    case 1<<dim:
        {
            // initialize as lower left corner ijk
            for (int i = 0; i < 4; i++)
            {
                ijk[i] = quad->ijk;
            }
            // ijk[0] = quad->ijk; //lower left corner of a quad
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
void LookUpTableForest<dim,USER_DATA>::refine(Quadrant<dim,USER_DATA>* quad, bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level))
{
    // printf("refine on thread %d, min_level: %d, max_level: %d, level: %d\n", omp_get_thread_num(), m_min_level, m_max_level, quad->level);
    // 
    // NOTE that do not use for loop in the recursion, it will slow down the calculation
    
    if(is_refine(this, quad, m_max_level))
    {
        // if no children, create children
        if(!quad->isHasChildren)
        {
            // make children and release user_data of parent, because the non-leaf quad doesn't need user data
            int length_child = 1<<(MAX_FOREST_LEVEL - quad->level -1); //Note that must -1, because this is the child length
            // z index = 0

            // quad->children = new Quadrant<dim,USER_DATA>*[1<<dim];

            // 1st child: lower left
            quad->children[0] = new Quadrant<dim,USER_DATA>;
            quad->children[0]->xyz[0] = quad->xyz[0];
            quad->children[0]->xyz[1] = quad->xyz[1];
            quad->children[0]->ijk    = quad->ijk;
            quad->children[0]->level  = quad->level+1; //only calculate once
            quad->children[0]->parent = quad;
            quad->children[0]->isHasChildren = false;
            // 2nd child: lower right
            quad->children[1] = new Quadrant<dim,USER_DATA>;
            quad->children[1]->xyz[0] = quad->xyz[0] + length_child*m_length_scale[0];
            quad->children[1]->xyz[1] = quad->xyz[1];
            quad->children[1]->ijk    = quad->ijk;
            quad->children[1]->ijk.i += length_child;
            quad->children[1]->level  = quad->children[0]->level;
            quad->children[1]->parent = quad;
            quad->children[1]->isHasChildren = false;
            // 3th child: upper left
            quad->children[2] = new Quadrant<dim,USER_DATA>;
            quad->children[2]->xyz[0] = quad->xyz[0];
            quad->children[2]->xyz[1] = quad->xyz[1] + length_child*m_length_scale[1];
            quad->children[2]->ijk    = quad->ijk;
            quad->children[2]->ijk.j += length_child;
            quad->children[2]->level  = quad->children[0]->level;
            quad->children[2]->parent = quad;
            quad->children[2]->isHasChildren = false;
            // 4th child: upper right
            quad->children[3] = new Quadrant<dim,USER_DATA>;
            quad->children[3]->xyz[0] = quad->children[1]->xyz[0];
            quad->children[3]->xyz[1] = quad->children[2]->xyz[1];
            quad->children[3]->ijk    = quad->ijk;
            quad->children[3]->ijk.i += length_child;
            quad->children[3]->ijk.j += length_child;
            quad->children[3]->level  = quad->children[0]->level;
            quad->children[3]->parent = quad;
            quad->children[3]->isHasChildren = false;
            if(m_data_size!=0) // only check once
            {
                quad->children[0]->user_data = new USER_DATA;
                quad->children[1]->user_data = new USER_DATA;
                quad->children[2]->user_data = new USER_DATA;
                quad->children[3]->user_data = new USER_DATA;
            }
            // 3D case: z index =1
            if(dim == 3)
            {
                // complete z coordinate of the first four children
                quad->children[0]->xyz[2] = quad->xyz[2];
                quad->children[1]->xyz[2] = quad->xyz[2];
                quad->children[2]->xyz[2] = quad->xyz[2];
                quad->children[3]->xyz[2] = quad->xyz[2];
                // ------- the last four children ------
                // 5th child: lower left
                quad->children[4] = new Quadrant<dim,USER_DATA>;
                quad->children[4]->xyz[0] = quad->xyz[0];
                quad->children[4]->xyz[1] = quad->xyz[1];
                quad->children[4]->xyz[2] = quad->xyz[2] + length_child*m_length_scale[2]; // only calculate once
                quad->children[4]->ijk    = quad->ijk;
                quad->children[4]->ijk.k += length_child;
                quad->children[4]->level  = quad->children[0]->level;
                quad->children[4]->parent = quad;
                quad->children[4]->isHasChildren = false;
                // 6th child: lower right
                quad->children[5] = new Quadrant<dim,USER_DATA>;
                quad->children[5]->xyz[0] = quad->children[1]->xyz[0];
                quad->children[5]->xyz[1] = quad->xyz[1];
                quad->children[5]->xyz[2] = quad->children[4]->xyz[2];
                quad->children[5]->ijk    = quad->ijk;
                quad->children[5]->ijk.i += length_child;
                quad->children[5]->ijk.k += length_child;
                quad->children[5]->level  = quad->children[0]->level;
                quad->children[5]->parent = quad;
                quad->children[5]->isHasChildren = false;
                // 7th child: upper left
                quad->children[6] = new Quadrant<dim,USER_DATA>;
                quad->children[6]->xyz[0] = quad->xyz[0];
                quad->children[6]->xyz[1] = quad->children[2]->xyz[1];
                quad->children[6]->xyz[2] = quad->children[4]->xyz[2];
                quad->children[6]->ijk    = quad->ijk;
                quad->children[6]->ijk.j += length_child;
                quad->children[6]->ijk.k += length_child;
                quad->children[6]->level  = quad->children[0]->level;
                quad->children[6]->parent = quad;
                quad->children[6]->isHasChildren = false;
                // 8th child: upper right
                quad->children[7] = new Quadrant<dim,USER_DATA>;
                quad->children[7]->xyz[0] = quad->children[1]->xyz[0];
                quad->children[7]->xyz[1] = quad->children[2]->xyz[1];
                quad->children[7]->xyz[2] = quad->children[4]->xyz[2];
                quad->children[7]->ijk    = quad->ijk;
                quad->children[7]->ijk.i += length_child;
                quad->children[7]->ijk.j += length_child;
                quad->children[7]->ijk.k += length_child;
                quad->children[7]->level  = quad->children[0]->level;
                quad->children[7]->parent = quad;
                quad->children[7]->isHasChildren = false;
                if(m_data_size!=0) // only check once
                {
                    quad->children[4]->user_data = new USER_DATA; 
                    quad->children[5]->user_data = new USER_DATA; 
                    quad->children[6]->user_data = new USER_DATA; 
                    quad->children[7]->user_data = new USER_DATA; 
                }
            }

            // delete parent data
            delete quad->user_data;
            quad->user_data = NULL;
            quad->isHasChildren = true;
        }

        // for (int i = 0; i < m_num_children; i++)
        // {
        #if USE_OMP == 1
            #pragma omp task shared(is_refine) //firstprivate(quad, m_max_level) //
        #endif
            refine(quad->children[0], is_refine);

        #if USE_OMP == 1
            #pragma omp task shared(is_refine) //firstprivate(quad, m_max_level)
        #endif
            refine(quad->children[1], is_refine);

        #if USE_OMP == 1
            #pragma omp task shared(is_refine) //firstprivate(quad, m_max_level)
        #endif
            refine(quad->children[2], is_refine);

        #if USE_OMP == 1
            #pragma omp task shared(is_refine) //firstprivate(quad, m_max_level)
        #endif
            refine(quad->children[3], is_refine);

            if(dim==3)
            {
            #if USE_OMP == 1
                #pragma omp task shared(is_refine)
            #endif
                refine(quad->children[4], is_refine);

            #if USE_OMP == 1
                #pragma omp task shared(is_refine)
            #endif
                refine(quad->children[5], is_refine);

            #if USE_OMP == 1
                #pragma omp task shared(is_refine)
            #endif
                refine(quad->children[6], is_refine);

            #if USE_OMP == 1
                #pragma omp task shared(is_refine)
            #endif
                refine(quad->children[7], is_refine);
            }
        #if USE_OMP == 1
            #pragma omp taskwait
        #endif
        // }
    }
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::refine(bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level))
{
    refine(&m_root, is_refine);
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::release_map2data()
{
    if(m_map_ijk2data.size()>0)
    {
        // release data and clear map
        for(auto &ijk2data : m_map_ijk2data)
        {
            delete[] ijk2data.second;
        }
        m_map_ijk2data.clear();
    }
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::construct_map2dat()
{
    // 0. safety check: if the m_map_ijk2data has been already created, release it to avoid repleated creation and memory leak
    release_map2data();

    // 1. get all leaves
    vector<Quadrant<dim,USER_DATA>* > leaves;
    getLeaves(leaves, m_num_quads, &m_root);
    m_num_leaves = leaves.size();
    
    // 2. allocate memory of data at all nodes of leaves
    Quad_index *ijk_nodes_quad = new Quad_index[m_num_node_per_quad]; // \todo 如果使用二阶插值，则需要更多节点，需要通过cellType进行判断：比如二维情况九点quad，那么需要限制max_level必须小于MAX_FOREST_LEVEL-2，不过这个好办，在构造函数里面判断一下进行安全检查就行
    for (size_t i = 0; i < leaves.size(); i++)
    {
        // get ijk of all nodes of a leaf quad
        get_ijk_nodes_quadrant(leaves[i], m_num_node_per_quad, ijk_nodes_quad);
        for (int i_node = 0; i_node < m_num_node_per_quad; i_node++)
        {
            if(m_map_ijk2data.count(ijk_nodes_quad[i_node]) == 0)
            {
                // cout<<ijk_nodes_quad[i_node].i<<" "<<ijk_nodes_quad[i_node].j<<" "<<ijk_nodes_quad[i_node].k<<endl;
                m_map_ijk2data[ijk_nodes_quad[i_node]] = new double[m_num_props];
            }
        }
    }

    // check key counts
    // cout<<"leaves: "<<leaves.size()<<endl;
    // cout<<"map size: "<<m_map_ijk2data.size()<<endl;

    // release pointer
    delete[] ijk_nodes_quad;
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::assemble_data(void (*cal_prop)(LookUpTableForest<dim,USER_DATA>* forest, std::map<Quad_index, double*>& map_ijk2data))
{
    STATUS("Assemble properties data to nodes of leaves ...");

    construct_map2dat();
    
    //3. call property calculation function to fill data to array
    cal_prop(this, m_map_ijk2data);

    STATUS("Assemble properties data done.");
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::getLeaves(vector<Quadrant<dim,USER_DATA>* >& leaves, long int& quad_counts, Quadrant<dim,USER_DATA>* quad)
{
    quad_counts++;

    if(quad->isHasChildren)
    {
        for (int i = 0; i < m_num_children; i++)
        {
            getLeaves(leaves, quad_counts, quad->children[i]);
        }
    }else
    {
        // quad->index = leaves.size();
        leaves.push_back(quad);
    }
}

// ASCII version
template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::write_to_vtk_v1(string filename, bool write_data, bool isNormalizeXYZ)
{
    clock_t start = clock();
    STATUS("Write to vtu file starting ...");
    Quadrant<dim,USER_DATA> *targetLeaf = NULL;
    vector<Quadrant<dim,USER_DATA>* > leaves;
    long int quad_counts = 0;
    getLeaves(leaves, quad_counts, &m_root);
    int num_points_per_cell = 1<<dim;
    int num_cells = leaves.size();
    int num_points = num_cells*num_points_per_cell;
    double length_cell[dim]; //physical length
    double length_ref_cell = 0;

    ofstream fout(filename);
    int CELL_TYPE = dim==2 ? 8 : 11;
    // write head
    STATUS("write header");
    fout<<"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<<endl;
    fout<<"  <UnstructuredGrid>"<<endl;
    fout<<"    <Piece NumberOfPoints=\""<<num_points<<"\" NumberOfCells=\""<<num_cells<<"\">"<<endl;
    // write point data 
    STATUS("write point data");
    fout<<"      <PointData>"<<endl;
    // // ---------- . phase index
    // fout<<"        <DataArray type=\"Int32\" Name=\"phaseIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Region;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . Rho
    // fout<<"        <DataArray type=\"Float32\" Name=\"T\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].T;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . Rho
    // fout<<"        <DataArray type=\"Float32\" Name=\"rho\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Rho;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . Rho_l
    // fout<<"        <DataArray type=\"Float32\" Name=\"rho_l\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Rho_l;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . Rho_v
    // fout<<"        <DataArray type=\"Float32\" Name=\"rho_v\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Rho_v;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . Rho_h
    // fout<<"        <DataArray type=\"Float32\" Name=\"rho_h\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Rho_h;}}
    // fout<<"\n        </DataArray>"<<endl;
    // ---------- . H
    // fout<<"        <DataArray type=\"Float32\" Name=\"H\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].H;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . H_l
    // fout<<"        <DataArray type=\"Float32\" Name=\"H_l\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].H_l;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . H_v
    // fout<<"        <DataArray type=\"Float32\" Name=\"H_v\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].H_v;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . H_h
    // fout<<"        <DataArray type=\"Float32\" Name=\"H_h\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].H_h;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . X
    // fout<<"        <DataArray type=\"Float32\" Name=\"X\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].X_wt;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . X_l
    // fout<<"        <DataArray type=\"Float32\" Name=\"X_l\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].X_l;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . X_v
    // fout<<"        <DataArray type=\"Float32\" Name=\"X_v\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].X_v;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . Mu
    // fout<<"        <DataArray type=\"Float32\" Name=\"Mu\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Mu;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . Mu_l
    // fout<<"        <DataArray type=\"Float32\" Name=\"Mu_l\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Mu_l;}}
    // fout<<"\n        </DataArray>"<<endl;
    // // ---------- . Mu_v
    // fout<<"        <DataArray type=\"Float32\" Name=\"Mu_v\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Mu_v;}}
    // fout<<"\n        </DataArray>"<<endl;
    // ----------
    fout<<"      </PointData>"<<endl;
    // write cell data 
    STATUS("write cell data");
    fout<<"      <CellData>"<<endl;
    // ---------- 1. phase index
    fout<<"        <DataArray type=\"Int32\" Name=\"phaseIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ fout<<" "<<leaves[i]->user_data->phaseRegion_cell;}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- 2. need refine
    fout<<"        <DataArray type=\"Int32\" Name=\"needRefine\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ fout<<" "<<leaves[i]->user_data->need_refine;}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- DEBUG. leaf index
    // fout<<"        <DataArray type=\"Int32\" Name=\"quadIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    // for (size_t i = 0; i < leaves.size(); i++){ fout<<" "<<leaves[i]->index;}
    // fout<<"\n        </DataArray>"<<endl;
    // ----------
    
    fout<<"      </CellData>"<<endl;
    // write points
    STATUS("write points (xyz)");
    fout<<"      <Points>"<<endl;
    fout<<"        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0.008838834896540746\" RangeMax=\"1.4053746938907843\">"<<endl;
    double scale=0; //scale = 0.001, keep some space between each quadrant
    double physical_length[dim]; //={m_xyz_max[0] - m_xyz_min[0], m_xyz_max[1] - m_xyz_min[1], m_xyz_max[2] - m_xyz_min[2]};
    for (int i = 0; i < dim; i++){ physical_length[i] = m_xyz_max[i] - m_xyz_min[i];}
    
    for (int i = 0; i < num_cells; i++)
    {
        length_ref_cell = 1<<(MAX_FOREST_LEVEL - leaves[i]->level);
        for (int j = 0; j < dim; j++)
        {
            length_cell[j] = m_length_scale[j] * length_ref_cell;
        }
        // ---------
        fout<<"        ";
        if (isNormalizeXYZ)
        {
            double z = m_constZ;
            if(dim==3)z = (leaves[i]->xyz[2]   - physical_length[2]*scale - m_xyz_min[2])/(m_xyz_max[2] - m_xyz_min[2]);
            fout<<" "<<(leaves[i]->xyz[0]                    + physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1]                  + physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<z; //Lower left
            fout<<" "<<(leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1]                  + physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<z; //lower right
            fout<<" "<<(leaves[i]->xyz[0]                    + physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<z; //upper left
            fout<<" "<<(leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<z; //upper right
            if(dim==3)
            {
                double z = (leaves[i]->xyz[2] + length_cell[2]   - physical_length[2]*scale - m_xyz_min[2])/(m_xyz_max[2] - m_xyz_min[2]);
                fout<<" "<<(leaves[i]->xyz[0]                    + physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1]                  + physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<z; //Lower left at z=zmax
                fout<<" "<<(leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1]                  + physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<z; //lower right at z=zmax
                fout<<" "<<(leaves[i]->xyz[0]                    + physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<z; //upper left at z=zmax
                fout<<" "<<(leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<z; //upper right at z=zmax
            }
        }else
        {
            fout<<" "<<leaves[i]->xyz[0]                    + physical_length[0]*scale<<" "<<leaves[i]->xyz[1]                  + physical_length[1]*scale<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //Lower left
            fout<<" "<<leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale<<" "<<leaves[i]->xyz[1]                  + physical_length[1]*scale<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //lower right
            fout<<" "<<leaves[i]->xyz[0]                    + physical_length[0]*scale<<" "<<leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //upper left
            fout<<" "<<leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale<<" "<<leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //upper right
            if(dim==3)
            {
                double z = leaves[i]->xyz[2] + length_cell[2]   - physical_length[2]*scale;
                fout<<" "<<leaves[i]->xyz[0]                    + physical_length[0]*scale<<" "<<leaves[i]->xyz[1]                  + physical_length[1]*scale<<" "<<z; //Lower left at z=zmax
                fout<<" "<<leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale<<" "<<leaves[i]->xyz[1]                  + physical_length[1]*scale<<" "<<z; //lower right at z=zmax
                fout<<" "<<leaves[i]->xyz[0]                    + physical_length[0]*scale<<" "<<leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale<<" "<<z; //upper left at z=zmax
                fout<<" "<<leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale<<" "<<leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale<<" "<<z; //upper right at z=zmax
            }
        }
        
        fout<<endl;
    }
    fout<<"        </DataArray>"<<endl;
    fout<<"      </Points>"<<endl;
    // write cells
    STATUS("write cells");
    fout<<"      <Cells>"<<endl;
    fout<<"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<num_cells-1<<"\">"<<endl;
    fout<<"        ";
    int point_id_cell = 0;
    for (int i = 0; i < num_cells; i++)
    {
        point_id_cell = i*num_points_per_cell;
        fout<<" "<<point_id_cell;
        fout<<" "<<point_id_cell+1;
        fout<<" "<<point_id_cell+2;
        fout<<" "<<point_id_cell+3;
        if(dim==3)
        {
            fout<<" "<<point_id_cell+4;
            fout<<" "<<point_id_cell+5;
            fout<<" "<<point_id_cell+6;
            fout<<" "<<point_id_cell+7;
        }
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
    STATUS_time("Write to vtu file done: "+filename, clock()-start);
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
    int num_points_per_cell = 1<<dim;
    int num_cells = leaves.size();
    int num_points = m_map_ijk2data.size();
    std::map<Quad_index, int> index_map_ijk2data;
    int ind =0;
    for(auto & m : m_map_ijk2data)
    {
        index_map_ijk2data[m.first] = ind;
        ind++;
    }
    double length_cell[dim]; //physical length
    double length_ref_cell = 0;
    Quad_index *ijk_nodes_quad = new Quad_index[m_num_node_per_quad]; // \todo 如果使用二阶插值，则需要更多节点，需要通过cellType进行判断：比如二维情况九点quad，那么需要限制max_level必须小于MAX_FOREST_LEVEL-2，不过这个好办，在构造函数里面判断一下进行安全检查就行

    ofstream fout(filename);
    int CELL_TYPE = dim==2 ? 8 : 11;
    // write head
    STATUS("write header");
    fout<<"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<<endl;
    fout<<"  <UnstructuredGrid>"<<endl;
    fout<<"    <Piece NumberOfPoints=\""<<num_points<<"\" NumberOfCells=\""<<num_cells<<"\">"<<endl;
    // write point data 
    STATUS("write point data");
    fout<<"      <PointData>"<<endl;
    ind = 0;
    for(auto &m : m_map_props)
    {
        fout<<"        <DataArray type=\"Float32\" Name=\""<<m.second.longName<<"\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
        fout<<"        ";
        for(auto &map_ijk2data : m_map_ijk2data)
        {
            fout<<" "<<map_ijk2data.second[ind];
        }
        fout<<"\n        </DataArray>"<<endl;
        ind++;
    }
    fout<<"      </PointData>"<<endl;
    // write cell data 
    STATUS("write cell data");
    fout<<"      <CellData>"<<endl;
    // ---------- 1. phase index
    fout<<"        <DataArray type=\"Int32\" Name=\"phaseIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ fout<<" "<<leaves[i]->user_data->phaseRegion_cell;}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- 2. need refine
    fout<<"        <DataArray type=\"Int32\" Name=\"needRefine\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ fout<<" "<<leaves[i]->user_data->need_refine;}
    fout<<"\n        </DataArray>"<<endl;
    fout<<"      </CellData>"<<endl;
    // write points
    STATUS("write points (xyz)");
    fout<<"      <Points>"<<endl;
    fout<<"        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0.008838834896540746\" RangeMax=\"1.4053746938907843\">"<<endl;
    double scale=0; //scale = 0.001, keep some space between each quadrant
    double physical_length[dim]; //={m_xyz_max[0] - m_xyz_min[0], m_xyz_max[1] - m_xyz_min[1], m_xyz_max[2] - m_xyz_min[2]};
    for (int i = 0; i < dim; i++){ physical_length[i] = m_xyz_max[i] - m_xyz_min[i];}
    for(auto &ijk2data : m_map_ijk2data)
    {
        fout<<"         "<<ijk2data.first.i<<" "<<ijk2data.first.j<<" "<<ijk2data.first.k<<endl;
    }
    fout<<"        </DataArray>"<<endl;
    fout<<"      </Points>"<<endl;
    // write cells
    STATUS("write cells");
    fout<<"      <Cells>"<<endl;
    fout<<"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<num_cells-1<<"\">"<<endl;
    // fout<<"        ";
    int point_id_cell = 0;
    for (int i = 0; i < num_cells; i++)
    {
        // get ijk of all nodes of a leaf quad
        get_ijk_nodes_quadrant(leaves[i], m_num_node_per_quad, ijk_nodes_quad);
        fout<<"         ";
        for (int i_node = 0; i_node < m_num_node_per_quad; i_node++)
        {
            fout<<index_map_ijk2data[ijk_nodes_quad[i_node]]<<" "; //low performance: distance(m_map_ijk2data.begin(),m_map_ijk2data.find(ijk_nodes_quad[i_node]))
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
    cout<<"All "<<m_num_leaves<<" leaves, "<<m_num_quads-m_num_leaves<<" non-leaf quads. "
    // <<"Estimate size "<<ceil(leaves.size()*(sizeof(FIELD_DATA<dim>) + sizeof(bool) + sizeof(double)*3 + sizeof(int))/1024.0/1024.0)<<" Mb"
    <<endl;
    // props info
    cout<<"Include "<<m_num_props<<" properties on each node."<<endl;
    int ind =0;
    for(auto & m : m_map_props)
    {
        cout<<"  "<<ind<<": "<<m.second.longName<<endl;
        ind++;
    }
    cout<<"================== Summary end ==================="<<endl;
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
void LookUpTableForest<dim,USER_DATA>::searchQuadrant(Quadrant<dim,USER_DATA> *quad_source, Quadrant<dim,USER_DATA> *&quad_target, double x_ref, double y_ref, double z_ref)
{
    // NOTE that do not use for loop in the recursion, it will slow down the calculation
    if(!quad_source->isHasChildren)
    {
        quad_target = quad_source;
        return;
    }else
    {
        int length_child = 1<<(MAX_FOREST_LEVEL - quad_source->children[0]->level); //any of a child to access the child level number, e.g., 0
        int childID_x = x_ref/length_child;
        int childID_y = y_ref/length_child;
        int childID_z = dim == 3 ? z_ref/length_child : 0;
        // some time if the input point at the upper right corner point, the childID_x==2, or childID_y==2, so need to make it to 1
        childID_x = childID_x >=2 ? 1 : childID_x;
        childID_y = childID_y >=2 ? 1 : childID_y;
        childID_z = childID_z >=2 ? 1 : childID_z;
        // cout<<"  childID_x: "<<childID_x<<", childID_y: "<<childID_y<<", level: "<<quad_source->level<<endl;
        searchQuadrant(quad_source->children[childID_z*4 + childID_y*2 + childID_x], quad_target,
            childID_x == 1 ? x_ref - length_child : x_ref, 
            childID_y == 1 ? y_ref - length_child : y_ref, 
            childID_z == 1 ? z_ref - length_child : z_ref
            );
    }
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::searchQuadrant(Quadrant<dim,USER_DATA> *&targetLeaf, double x, double y, double z)
{
    // NOTE that do not use for loop in the recursion, it will slow down the calculation
    double x_ref = (x - m_xyz_min[0])/m_length_scale[0];
    double y_ref = (y - m_xyz_min[1])/m_length_scale[1];
    double z_ref = dim == 3 ? (z - m_xyz_min[2])/m_length_scale[2] : 0;
    searchQuadrant(&m_root, targetLeaf, x_ref, y_ref, z_ref);
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
