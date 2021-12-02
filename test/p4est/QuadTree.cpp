
#include "QuadTree.h"

template <int dim, typename USER_DATA> 
cForest<dim,USER_DATA>::cForest(double xyz_min[3], double xyz_max[3], int max_level, size_t data_size)
:
m_num_children(1<<dim),
m_data_size(data_size),
m_max_level(max_level)
{
    double length_forest = (1<<MAX_FOREST_LEVEL);
    for (size_t i = 0; i < 3; i++)
    {
        m_xyz_max[i] = xyz_max[i];
        m_xyz_min[i] = xyz_min[i];
        length_scale[i] = (m_xyz_max[i] - m_xyz_min[i])/length_forest;
    }
    init_Root(m_root);
}

template <int dim, typename USER_DATA> 
cForest<dim,USER_DATA>::~cForest()
{
    release_quadrant_data(&m_root);
    release_children(&m_root);
}

template <int dim, typename USER_DATA> 
void cForest<dim,USER_DATA>::release_children(Quadrant<dim,USER_DATA>* quad)
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
            }
        }
        
    }else
    {
        return;
    }
}

template <int dim, typename USER_DATA> 
void cForest<dim,USER_DATA>::release_quadrant_data(Quadrant<dim,USER_DATA>* quad)
{
    if(quad->isHasChildren)
    {
        for (size_t i = 0; i < m_num_children; i++)
        {
            release_quadrant_data(quad->children[i]);
        }
    }else if (quad->user_data)
    {
        delete quad->user_data;
    }
}

template <int dim, typename USER_DATA> 
void cForest<dim,USER_DATA>::init_Root(Quadrant<dim,USER_DATA>& quad)
{
    quad.xyz[0] = m_xyz_min[0];
    quad.xyz[1] = m_xyz_min[1];
    quad.xyz[2] = m_xyz_min[2];
    quad.level  = 0;
    quad.isHasChildren = false;
    if(m_data_size!=0) quad.user_data = new USER_DATA;  //only allocate memory if it is a leaf, this will be released when a quadrent is refined.
}

template <int dim, typename USER_DATA>
void cForest<dim,USER_DATA>::refine(Quadrant<dim,USER_DATA>* quad, bool (*is_refine)(cForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level))
{
    // if(quad->level > m_max_level)return;

    if(is_refine(this, quad, m_max_level))
    {
        // if no children, create children
        if(!quad->isHasChildren)
        {
            // make children and release user_data, because the non-leaf quad doesn't need user data
            int length_child = 1<<(MAX_FOREST_LEVEL - quad->level -1); //Note that must -1, because this is the child length
            // first child: lower left
            quad->children[0] = new Quadrant<dim,USER_DATA>;
            quad->children[0]->xyz[0] = quad->xyz[0];
            quad->children[0]->xyz[1] = quad->xyz[1];
            quad->children[0]->xyz[2] = quad->xyz[2];
            quad->children[0]->level = quad->level+1;
            quad->children[0]->isHasChildren = false;
            if(m_data_size!=0) quad->children[0]->user_data = new USER_DATA; 
            // second child: lower right
            quad->children[1] = new Quadrant<dim,USER_DATA>;
            quad->children[1]->xyz[0] = quad->xyz[0] + length_child*length_scale[0];
            quad->children[1]->xyz[1] = quad->xyz[1];
            quad->children[1]->xyz[2] = quad->xyz[2];
            quad->children[1]->level = quad->level+1;
            quad->children[1]->isHasChildren = false;
            if(m_data_size!=0) quad->children[1]->user_data = new USER_DATA; 
            // third child: upper left
            quad->children[2] = new Quadrant<dim,USER_DATA>;
            quad->children[2]->xyz[0] = quad->xyz[0];
            quad->children[2]->xyz[1] = quad->xyz[1] + length_child*length_scale[1];
            quad->children[2]->xyz[2] = quad->xyz[2];
            quad->children[2]->level = quad->level+1;
            quad->children[2]->isHasChildren = false;
            if(m_data_size!=0) quad->children[2]->user_data = new USER_DATA; 
            // second child: upper right
            quad->children[3] = new Quadrant<dim,USER_DATA>;
            quad->children[3]->xyz[0] = quad->xyz[0] + length_child*length_scale[0];
            quad->children[3]->xyz[1] = quad->xyz[1] + length_child*length_scale[1];
            quad->children[3]->xyz[2] = quad->xyz[2];
            quad->children[3]->level = quad->level+1;
            quad->children[3]->isHasChildren = false;
            if(m_data_size!=0) quad->children[3]->user_data = new USER_DATA; 
            //  if dim==3
            //\todo 补充三维的情况

            // delete parent data
            delete quad->user_data;
            quad->user_data = NULL;
            quad->isHasChildren = true;
        }

        for (size_t i = 0; i < m_num_children; i++)
        {
            refine(quad->children[i], is_refine);
        }
    }
}

template <int dim, typename USER_DATA>
void cForest<dim,USER_DATA>::refine(bool (*is_refine)(cForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level))
{
    refine(&m_root, is_refine);
}

template <int dim, typename USER_DATA>
void cForest<dim,USER_DATA>::getLeaves(vector<Quadrant<dim,USER_DATA>* >& leaves, Quadrant<dim,USER_DATA>* quad)
{
    if(quad->isHasChildren)
    {
        for (size_t i = 0; i < m_num_children; i++)
        {
            getLeaves(leaves, quad->children[i]);
        }
    }else
    {
        leaves.push_back(quad);
    }
}

// ASCII version
template <int dim, typename USER_DATA>
void cForest<dim,USER_DATA>::write_to_vtk(string filename, bool write_data, bool isNormalizeXYZ)
{
    vector<Quadrant<dim,USER_DATA>* > leaves;
    getLeaves(leaves, &m_root);
    int num_points_per_cell = 1<<dim;
    int num_cells = leaves.size();
    int num_points = num_cells*num_points_per_cell;
    double length_cell[3]; //physical length
    double length_ref_cell = 0;

    ofstream fout(filename);
    // write head
    fout<<"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<<endl;
    fout<<"  <UnstructuredGrid>"<<endl;
    fout<<"    <Piece NumberOfPoints=\""<<num_points<<"\" NumberOfCells=\""<<num_cells<<"\">"<<endl;
    // write point data 
    fout<<"      <PointData>"<<endl;
    fout<<"      </PointData>"<<endl;
    // write cell data 
    fout<<"      <CellData>"<<endl;
    fout<<"        <DataArray type=\"Int32\" Name=\"phaseIndex\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">"<<endl;
    fout<<"        "<<endl;
    for (size_t i = 0; i < leaves.size(); i++)
    {
        fout<<" "<<leaves[i]->user_data->phaseRegion_cell;
    }
    fout<<endl;
    fout<<"        </DataArray>"<<endl;
    fout<<"      </CellData>"<<endl;
    // write points
    fout<<"      <Points>"<<endl;
    fout<<"        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0.008838834896540746\" RangeMax=\"1.4053746938907843\">"<<endl;
    double scale=0;
    double physical_length[3] ={m_xyz_max[0] - m_xyz_min[0], m_xyz_max[1] - m_xyz_min[1], m_xyz_max[2] - m_xyz_min[2]};
    for (size_t i = 0; i < num_cells; i++)
    {
        length_ref_cell = 1<<(MAX_FOREST_LEVEL - leaves[i]->level);
        for (int j = 0; j < 3; j++)
        {
            length_cell[j] = length_scale[j] * length_ref_cell;
        }
        fout<<"        ";
        if (isNormalizeXYZ)
        {
            fout<<" "<<(leaves[i]->xyz[0]                    + physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1]                  + physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<leaves[i]->xyz[2]; //Lower left
            fout<<" "<<(leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1]                  + physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<leaves[i]->xyz[2]; //lower right
            fout<<" "<<(leaves[i]->xyz[0]                    + physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<leaves[i]->xyz[2]; //upper left
            fout<<" "<<(leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale - m_xyz_min[0])/(m_xyz_max[0] - m_xyz_min[0])<<" "<<(leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale - m_xyz_min[1])/(m_xyz_max[1] - m_xyz_min[1])<<" "<<leaves[i]->xyz[2]; //upper right
        }else
        {
            fout<<" "<<leaves[i]->xyz[0]                    + physical_length[0]*scale<<" "<<leaves[i]->xyz[1]                  + physical_length[1]*scale<<" "<<leaves[i]->xyz[2]; //Lower left
            fout<<" "<<leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale<<" "<<leaves[i]->xyz[1]                  + physical_length[1]*scale<<" "<<leaves[i]->xyz[2]; //lower right
            fout<<" "<<leaves[i]->xyz[0]                    + physical_length[0]*scale<<" "<<leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale<<" "<<leaves[i]->xyz[2]; //upper left
            fout<<" "<<leaves[i]->xyz[0] + length_cell[0]   - physical_length[0]*scale<<" "<<leaves[i]->xyz[1] + length_cell[1] - physical_length[1]*scale<<" "<<leaves[i]->xyz[2]; //upper right
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
}

template <int dim, typename USER_DATA>
void cForest<dim,USER_DATA>::get_quadrant_physical_length(int level, double physical_length[3])
{
    int length_ref = 1<<(MAX_FOREST_LEVEL - level);
    for (int i = 0; i < 3; i++)
    {
        physical_length[i] = length_scale[i] * length_ref;
    }
}

int phaseRegion_sin_func(double x, double y)
{
    double value = sin(x);
    int region = 0;
    if(value==y)
    {
        region = 0;
    }else if(y > value)
    {
        region = 1;
    }else
    {
        region = -1;
    }
    return region;
}

int phaseRegion_H2ONaCl_constantX(double T_C, double p_bar, double X = 3.2)
{
    double xv, xl;
    H2ONaCl::PhaseRegion phaseregion=eos.findPhaseRegion(T_C, p_bar, eos.Wt2Mol(X/100.0),xl,xv);

    return phaseregion;
}

void calProp_H2ONaCl_consX(H2ONaCl::PROP_H2ONaCl& prop, double T_C, double p_bar, double X = 3.2)
{
    prop = eos.prop_pTX(p_bar*1E5, T_C+273.15, X/100.0);
}

template <int dim, typename USER_DATA>
bool refine_fn(cForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level)
{
    if(quad->isHasChildren) return true; //if a quad has children, of course it need refine, but we don't need do anything at here, just return true.

    USER_DATA       *data = (USER_DATA *) quad->user_data;
    bool need_refine = false;
    double physical_length[3];
    forest->get_quadrant_physical_length(quad->level, physical_length);
    const int num_sample_x =2; //最简单的情况就是只取xmin, xmax作为采样点判断这些采样点的函数计算返回值(flat)是否全部相等.但是有时候会有漏掉的情况，所以可以考虑在这里加密采样
    const int num_sample_y =2;
    double dx_qua = physical_length[0] / (num_sample_x - 1.0);
    double dy_qua = physical_length[1] / (num_sample_y - 1.0);
    double x_qua, x_ref, y_qua, y_ref;
    double xyz_tmp[3];
    int regionIndex[num_sample_x*num_sample_y];
    for (int iy = 0; iy < num_sample_y; iy++)
    {
        y_qua = quad->xyz[1] + dy_qua*iy;
        for (int ix = 0; ix < num_sample_x; ix++)
        {
            x_qua = quad->xyz[0] + dx_qua*ix;
            // regionIndex[iy*num_sample_x + ix] = phaseRegion_sin_func(x_qua, y_qua); 
            regionIndex[iy*num_sample_x + ix] = phaseRegion_H2ONaCl_constantX(x_qua, y_qua);
            // cout<<"***** "<<phaseRegion_H2ONaCl_constantX(x_qua, y_qua)<<endl;
        }
    }
    // ========== 1. refinement check for phase index ============
    bool isSame_phaseIndex = true;
    for (int i = 1; i < num_sample_x*num_sample_y; i++)
    {
        isSame_phaseIndex = (isSame_phaseIndex && (regionIndex[0] == regionIndex[i]));
    }
    if(!isSame_phaseIndex) //如果采样点的phase index全相等，则不refine；否则refine
    {
        need_refine = true;
    }else
    {
        need_refine = false;
    }
    // ========================================================
    // calculate properties: four vertices and one midpoint
    calProp_H2ONaCl_consX(data->prop_point[0], quad->xyz[0],                      quad->xyz[1]);                         //xmin,ymin
    calProp_H2ONaCl_consX(data->prop_point[1], quad->xyz[0] + physical_length[0], quad->xyz[1]);                         //xmax,ymin
    calProp_H2ONaCl_consX(data->prop_point[2], quad->xyz[0],                      quad->xyz[1] + physical_length[1]);    //xmin,ymax
    calProp_H2ONaCl_consX(data->prop_point[3], quad->xyz[0] + physical_length[0], quad->xyz[1] + physical_length[1]);    //xmax,ymax
    calProp_H2ONaCl_consX(data->prop_cell, quad->xyz[0] + physical_length[0]/2.0, quad->xyz[1] + physical_length[1]/2.0); //xc,yc
    data->phaseRegion_point[0] = regionIndex[0]; //phase index
    data->phaseRegion_point[1] = regionIndex[num_sample_x-1];
    data->phaseRegion_point[2] = regionIndex[num_sample_x*num_sample_y-num_sample_x];
    data->phaseRegion_point[3] = regionIndex[num_sample_x*num_sample_y-1];
    data->phaseRegion_cell     = phaseRegion_H2ONaCl_constantX(quad->xyz[0] + physical_length[0]/2.0, quad->xyz[1] + physical_length[1]/2.0);
    
    // set some special indicator if cell need refine
    if(need_refine)
    {
        data->phaseRegion_cell = H2ONaCl::MixPhaseRegion;
        data->need_refine = NeedRefine_PB_L2V; //\todo use enum to mark it to different phase boundaries

    }else
    {
        data->need_refine = NeedRefine_NoNeed;
    }

    // return refine indicator
    if(quad->level > forest->m_max_level)return false;
    if(need_refine)return true;

    return false;
}

template <int dim, typename USER_DATA>
bool refine_uniform(cForest<dim,USER_DATA>* forest, Quadrant<dim,USER_DATA>* quad, int max_level)
{
    if(quad->level <= 3)return true;

    return false;
}

int main()
{
    // double xyzmin[3] = {0,-2,0};
    // double xyzmax[3] = {2*3.141592653, 1.5, 100};
    double xyzmin[3] = {1,5, 0};
    double xyzmax[3] = {700, 400, 0};

    int max_level = 12;
    const int dim =2;
    cForest<dim, FIELD_DATA<dim> > forest(xyzmin, xyzmax, max_level, sizeof(FIELD_DATA<dim>));
    // refine 
    forest.refine(refine_uniform);
    forest.refine(refine_fn);
    
    forest.write_to_vtk("quadTree.vtu");

    // std::string dummy;
    // std::cout << "Enter to continue..." << std::endl;
    // std::getline(std::cin, dummy);
    return 0;
}