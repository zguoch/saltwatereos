
#include "QuadTree.h"

template <int dim, typename USER_DATA> 
cForest<dim,USER_DATA>::cForest(double xyz_min[3], double xyz_max[3], int max_level, size_t data_size)
:
m_num_children(1<<dim),
m_data_size(data_size),
m_max_level(max_level)
{
    for (size_t i = 0; i < 3; i++)
    {
        m_xyz_max[i] = xyz_max[i];
        m_xyz_min[i] = xyz_min[i];
    }
    
    init_Root(m_root);
}
template <int dim, typename USER_DATA> 
cForest<dim,USER_DATA>::~cForest()
{
    release_quadrant_data(&m_root);
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
    quad.xyz[0] = 0;
    quad.xyz[1] = 0;
    quad.xyz[2] = 0;
    quad.level  = 0;
    quad.isHasChildren = false;
    if(m_data_size!=0) quad.user_data = new USER_DATA;  //only allocate memory if it is a leaf, this will be released when a quadrent is refined.
}

template <int dim, typename USER_DATA>
void cForest<dim,USER_DATA>::refine_uniform(Quadrant<dim,USER_DATA>* quad, int max_level_refine)
{
    if(quad->level >= max_level_refine)return;

    // if no children, create children
    if(!quad->isHasChildren)
    {
        // make children and release user_data, because the non-leaf quad doesn't need user data
        int length_child = 1<<(m_max_level - quad->level -1);
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
        quad->children[1]->xyz[0] = quad->xyz[0] + length_child;
        quad->children[1]->xyz[1] = quad->xyz[1];
        quad->children[1]->xyz[2] = quad->xyz[2];
        quad->children[1]->level = quad->level+1;
        quad->children[1]->isHasChildren = false;
        if(m_data_size!=0) quad->children[1]->user_data = new USER_DATA; 
        // third child: upper left
        quad->children[2] = new Quadrant<dim,USER_DATA>;
        quad->children[2]->xyz[0] = quad->xyz[0];
        quad->children[2]->xyz[1] = quad->xyz[1] + length_child;
        quad->children[2]->xyz[2] = quad->xyz[2];
        quad->children[2]->level = quad->level+1;
        quad->children[2]->isHasChildren = false;
        if(m_data_size!=0) quad->children[2]->user_data = new USER_DATA; 
        // second child: upper right
        quad->children[3] = new Quadrant<dim,USER_DATA>;
        quad->children[3]->xyz[0] = quad->xyz[0] + length_child;
        quad->children[3]->xyz[1] = quad->xyz[1] + length_child;
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
        refine_uniform(quad->children[i], max_level_refine);
    }
}

template <int dim, typename USER_DATA>
void cForest<dim,USER_DATA>::refine_uniform(int max_level_refine)
{
    refine_uniform(&m_root, max_level_refine);
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

template <int dim, typename USER_DATA>
void cForest<dim,USER_DATA>::write_to_vtk(const Quadrant<dim,USER_DATA>* quad)
{

}

// ASCII version
template <int dim, typename USER_DATA>
void cForest<dim,USER_DATA>::write_to_vtk(string filename)
{
    double length_forest = 1<<m_max_level; //reference system is a square
    double length_scale[3]; //map real xyz to reference system
    double length_cell[3];
    for (int i = 0; i < 3; i++)
    {
        length_scale[i] = (m_xyz_max[i] - m_xyz_min[i])/length_forest;
    }
    vector<Quadrant<dim,USER_DATA>* > leaves;
    getLeaves(leaves, &m_root);
    int num_points_per_cell = 1<<dim;
    int num_cells = leaves.size();
    int num_points = num_cells*num_points_per_cell;

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
    fout<<"      </CellData>"<<endl;
    // write points
    fout<<"      <Points>"<<endl;
    fout<<"        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0.008838834896540746\" RangeMax=\"1.4053746938907843\">"<<endl;
    double length_ref_cell = 0;
    double xyz_cell[3]; //lower left corner
    for (size_t i = 0; i < num_cells; i++)
    {
        length_ref_cell = 1<<(m_max_level - leaves[i]->level);
        for (int j = 0; j < 3; j++)
        {
            length_cell[j] = length_scale[j] * length_ref_cell;
            xyz_cell[j] = leaves[i]->xyz[j]*length_scale[j] + m_xyz_min[j];
        }
        fout<<"         "<<xyz_cell[0]<<" "<<xyz_cell[1]<<" "<<xyz_cell[2]; //Lower left
        fout<<" "<<xyz_cell[0] + length_cell[0]<<" "<<xyz_cell[1]<<" "<<xyz_cell[2]; //lower right
        fout<<" "<<xyz_cell[0]<<" "<<xyz_cell[1] + length_cell[1]<<" "<<xyz_cell[2]; //upper left
        fout<<" "<<xyz_cell[0] + length_cell[0]<<" "<<xyz_cell[1] + length_cell[1]<<" "<<xyz_cell[2]; //upper right
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

int main()
{
    double xyzmin[3] = {0,0,0};
    double xyzmax[3] = {700, 400, 100};
    typedef double type_user_data;
    int max_level = 5;
    cForest<2, type_user_data> forest(xyzmin, xyzmax, max_level, sizeof(type_user_data));

    // refine 
    forest.refine_uniform(4);

    forest.write_to_vtk("quadTree.vtu");

    return 0;
}