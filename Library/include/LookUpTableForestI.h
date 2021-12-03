
template <int dim, typename USER_DATA> 
LookUpTableForest<dim,USER_DATA>::LookUpTableForest(double xyz_min[dim], double xyz_max[dim], int max_level, size_t data_size, void* eosPointer)
:
m_constZ(0) //this is not used in 3D table
{
    init(xyz_min, xyz_max, max_level, data_size, eosPointer);
}

template <int dim, typename USER_DATA> 
LookUpTableForest<dim,USER_DATA>::LookUpTableForest(double xyz_min[dim], double xyz_max[dim], double constZ, int max_level, size_t data_size, void* eosPointer)
:
m_constZ(constZ)
{
    init(xyz_min, xyz_max, max_level, data_size, eosPointer);
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::init(double xyz_min[dim], double xyz_max[dim], int max_level, size_t data_size, void* eosPointer)
{
    m_eosPointer = eosPointer;
    m_num_children = 1<<dim;
    m_data_size = data_size;
    m_min_level = 0;
    m_max_level = max_level;
    RMSD_Rho_min = 5;
    RMSD_H_min = 1000;

    double length_forest = (1<<MAX_FOREST_LEVEL);
    for (size_t i = 0; i < dim; i++)
    {
        m_xyz_max[i] = xyz_max[i];
        m_xyz_min[i] = xyz_min[i];
        length_scale[i] = (m_xyz_max[i] - m_xyz_min[i])/length_forest;
    }
    init_Root(m_root);
}

template <int dim, typename USER_DATA> 
LookUpTableForest<dim,USER_DATA>::~LookUpTableForest()
{
    release_quadrant_data(&m_root);
    release_children(&m_root);
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::release_children(Quadrant<dim,USER_DATA>* quad)
{
    if(quad->isHasChildren)
    {
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
        }
        
    }else
    {
        return;
    }
}

template <int dim, typename USER_DATA> 
void LookUpTableForest<dim,USER_DATA>::release_quadrant_data(Quadrant<dim,USER_DATA>* quad)
{
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
    quad.level  = 0;
    quad.isHasChildren = false;
    if(m_data_size!=0) quad.user_data = new USER_DATA;  //only allocate memory if it is a leaf, this will be released when a quadrent is refined.
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::refine(Quadrant<dim,USER_DATA>* quad, bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level))
{
    // NOTE that do not use for loop in the recursion, it will slow down the calculation
    if(is_refine(this, quad, m_max_level))
    {
        // if no children, create children
        if(!quad->isHasChildren)
        {
            // make children and release user_data of parent, because the non-leaf quad doesn't need user data
            int length_child = 1<<(MAX_FOREST_LEVEL - quad->level -1); //Note that must -1, because this is the child length
            // z index = 0
            // 1st child: lower left
            quad->children[0] = new Quadrant<dim,USER_DATA>;
            quad->children[0]->xyz[0] = quad->xyz[0];
            quad->children[0]->xyz[1] = quad->xyz[1];
            quad->children[0]->level  = quad->level+1; //only calculate once
            quad->children[0]->isHasChildren = false;
            // 2nd child: lower right
            quad->children[1] = new Quadrant<dim,USER_DATA>;
            quad->children[1]->xyz[0] = quad->xyz[0] + length_child*length_scale[0];
            quad->children[1]->xyz[1] = quad->xyz[1];
            quad->children[1]->level  = quad->children[0]->level;
            quad->children[1]->isHasChildren = false;
            // 3th child: upper left
            quad->children[2] = new Quadrant<dim,USER_DATA>;
            quad->children[2]->xyz[0] = quad->xyz[0];
            quad->children[2]->xyz[1] = quad->xyz[1] + length_child*length_scale[1];
            quad->children[2]->level  = quad->children[0]->level;
            quad->children[2]->isHasChildren = false;
            // 4th child: upper right
            quad->children[3] = new Quadrant<dim,USER_DATA>;
            quad->children[3]->xyz[0] = quad->children[1]->xyz[0];
            quad->children[3]->xyz[1] = quad->children[2]->xyz[1];
            quad->children[3]->level  = quad->children[0]->level;
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
                quad->children[4]->xyz[2] = quad->xyz[2] + length_child*length_scale[2]; // only calculate once
                quad->children[4]->level  = quad->children[0]->level;
                quad->children[4]->isHasChildren = false;
                // 6th child: lower right
                quad->children[5] = new Quadrant<dim,USER_DATA>;
                quad->children[5]->xyz[0] = quad->children[1]->xyz[0];
                quad->children[5]->xyz[1] = quad->xyz[1];
                quad->children[5]->xyz[2] = quad->children[4]->xyz[2];
                quad->children[5]->level  = quad->children[0]->level;
                quad->children[5]->isHasChildren = false;
                // 7th child: upper left
                quad->children[6] = new Quadrant<dim,USER_DATA>;
                quad->children[6]->xyz[0] = quad->xyz[0];
                quad->children[6]->xyz[1] = quad->children[2]->xyz[1];
                quad->children[6]->xyz[2] = quad->children[4]->xyz[2];
                quad->children[6]->level  = quad->children[0]->level;
                quad->children[6]->isHasChildren = false;
                // 8th child: upper right
                quad->children[7] = new Quadrant<dim,USER_DATA>;
                quad->children[7]->xyz[0] = quad->children[1]->xyz[0];
                quad->children[7]->xyz[1] = quad->children[2]->xyz[1];
                quad->children[7]->xyz[2] = quad->children[4]->xyz[2];
                quad->children[7]->level  = quad->children[0]->level;
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

        for (int i = 0; i < m_num_children; i++)
        {
            refine(quad->children[i], is_refine);
        }
    }
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::refine(bool (*is_refine)(LookUpTableForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level))
{
    refine(&m_root, is_refine);
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::getLeaves(vector<Quadrant<dim,USER_DATA>* >& leaves, Quadrant<dim,USER_DATA>* quad)
{
    if(quad->isHasChildren)
    {
        for (int i = 0; i < m_num_children; i++)
        {
            getLeaves(leaves, quad->children[i]);
        }
    }else
    {
        quad->index = leaves.size();
        leaves.push_back(quad);
    }
}

// ASCII version
template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::write_to_vtk(string filename, bool write_data, bool isNormalizeXYZ)
{
    clock_t start = clock();
    STATUS("Write to vtu file starting ...");
    Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    vector<Quadrant<dim,USER_DATA>* > leaves;
    getLeaves(leaves, &m_root);
    int num_points_per_cell = 1<<dim;
    int num_cells = leaves.size();
    int num_points = num_cells*num_points_per_cell;
    double length_cell[dim]; //physical length
    double length_ref_cell = 0;

    ofstream fout(filename);
    // write head
    fout<<"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<<endl;
    fout<<"  <UnstructuredGrid>"<<endl;
    fout<<"    <Piece NumberOfPoints=\""<<num_points<<"\" NumberOfCells=\""<<num_cells<<"\">"<<endl;
    // write point data 
    fout<<"      <PointData>"<<endl;
    // ---------- . phase index
    fout<<"        <DataArray type=\"Int32\" Name=\"phaseIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->phaseRegion_point[j];}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . Rho
    fout<<"        <DataArray type=\"Float32\" Name=\"rho\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Rho;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . Rho_l
    fout<<"        <DataArray type=\"Float32\" Name=\"rho_l\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Rho_l;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . Rho_v
    fout<<"        <DataArray type=\"Float32\" Name=\"rho_v\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Rho_v;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . Rho_h
    fout<<"        <DataArray type=\"Float32\" Name=\"rho_h\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Rho_h;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . H
    fout<<"        <DataArray type=\"Float32\" Name=\"H\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].H;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . H_l
    fout<<"        <DataArray type=\"Float32\" Name=\"H_l\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].H_l;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . H_v
    fout<<"        <DataArray type=\"Float32\" Name=\"H_v\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].H_v;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . H_h
    fout<<"        <DataArray type=\"Float32\" Name=\"H_h\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].H_h;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . X
    fout<<"        <DataArray type=\"Float32\" Name=\"X\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].X_wt;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . X_l
    fout<<"        <DataArray type=\"Float32\" Name=\"X_l\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].X_l;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . X_v
    fout<<"        <DataArray type=\"Float32\" Name=\"X_v\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].X_v;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . Mu
    fout<<"        <DataArray type=\"Float32\" Name=\"Mu\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Mu;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . Mu_l
    fout<<"        <DataArray type=\"Float32\" Name=\"Mu_l\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Mu_l;}}
    fout<<"\n        </DataArray>"<<endl;
    // ---------- . Mu_v
    fout<<"        <DataArray type=\"Float32\" Name=\"Mu_v\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ for (int j = 0; j < num_points_per_cell; j++){fout<<" "<<leaves[i]->user_data->prop_point[j].Mu_v;}}
    fout<<"\n        </DataArray>"<<endl;
    // ----------
    fout<<"      </PointData>"<<endl;
    // write cell data 
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
    fout<<"        <DataArray type=\"Int32\" Name=\"quadIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n        ";
    for (size_t i = 0; i < leaves.size(); i++){ fout<<" "<<leaves[i]->index;}
    fout<<"\n        </DataArray>"<<endl;
    // ----------
    
    fout<<"      </CellData>"<<endl;
    // write points
    fout<<"      <Points>"<<endl;
    fout<<"        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0.008838834896540746\" RangeMax=\"1.4053746938907843\">"<<endl;
    double scale=0; //scale = 0.001, keep some space between each quadrant
    double physical_length[dim]; //={m_xyz_max[0] - m_xyz_min[0], m_xyz_max[1] - m_xyz_min[1], m_xyz_max[2] - m_xyz_min[2]};
    for (int i = 0; i < dim; i++){ physical_length[i] = m_xyz_max[i] - m_xyz_min[i];}
    
    for (int i = 0; i < num_cells; i++)
    {
        length_ref_cell = 1<<(MAX_FOREST_LEVEL - leaves[i]->level);
        for (int j = 0; j < 3; j++)
        {
            length_cell[j] = length_scale[j] * length_ref_cell;
        }
        // ---------
        fout<<"        ";
        if (isNormalizeXYZ)
        {
            fout<<" "<<(leaves[i]->xyz[0]                    + physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1]                  + physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //Lower left
            fout<<" "<<(leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1]                  + physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //lower right
            fout<<" "<<(leaves[i]->xyz[0]                    + physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //upper left
            fout<<" "<<(leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //upper right
        }else
        {
            fout<<" "<<leaves[i]->xyz[0]                    + physical_length[0]*scale<<" "<<leaves[i]->xyz[1]                  + physical_length[1]*scale<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //Lower left
            fout<<" "<<leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale<<" "<<leaves[i]->xyz[1]                  + physical_length[1]*scale<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //lower right
            fout<<" "<<leaves[i]->xyz[0]                    + physical_length[0]*scale<<" "<<leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //upper left
            fout<<" "<<leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale<<" "<<leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale<<" "<<(dim == 3 ? leaves[i]->xyz[2] : m_constZ); //upper right
        }
        
        fout<<endl;
    }
    fout<<"        </DataArray>"<<endl;
    fout<<"      </Points>"<<endl;
    // write cells
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
    for (int i = 0; i < num_cells; i++)fout<<" "<<8; //dim dependent
    fout<<endl;
    fout<<"        </DataArray>"<<endl;
    fout<<"      </Cells>"<<endl;

    // write footer
    fout<<"    </Piece>"<<endl;
    fout<<"  </UnstructuredGrid>"<<endl;
    fout<<"</VTKFile>"<<endl;
    fout.close();
    STATUS_time("Write to vtu file done", clock()-start);
}

template <int dim, typename USER_DATA>
void LookUpTableForest<dim,USER_DATA>::get_quadrant_physical_length(int level, double physical_length[dim])
{
    int length_ref = 1<<(MAX_FOREST_LEVEL - level);
    for (int i = 0; i < dim; i++)
    {
        physical_length[i] = length_scale[i] * length_ref;
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
    double x_ref = (x - m_xyz_min[0])/length_scale[0];
    double y_ref = (y - m_xyz_min[1])/length_scale[1];
    double z_ref = dim == 3 ? (z - m_xyz_min[2])/length_scale[2] : 0;
    searchQuadrant(&m_root, targetLeaf, x_ref, y_ref, z_ref);
}

template <int dim, typename USER_DATA>
int LookUpTableForest<dim,USER_DATA>::searchQuadrant(double x, double y, double z)
{
    double x_ref = (x - m_xyz_min[0])/length_scale[0];
    double y_ref = (y - m_xyz_min[1])/length_scale[1];
    double z_ref = dim == 3 ? (z - m_xyz_min[2])/length_scale[2] : 0;
    Quadrant<dim,USER_DATA>* targetLeaf;
    searchQuadrant(&m_root, targetLeaf, x_ref, y_ref, z_ref);
    return targetLeaf->index;
}