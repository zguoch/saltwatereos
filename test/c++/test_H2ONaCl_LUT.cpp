#include "H2ONaCl.H"
#include <iostream>
H2ONaCl::cH2ONaCl eos;

typedef LOOKUPTABLE_FOREST::FIELD_DATA<2> FIELD_DATA_2D;
typedef LOOKUPTABLE_FOREST::Quadrant<2,FIELD_DATA_2D> Quad_2D;

double m_num_children = 1<<2;

void writeForest(FILE* fpout, Quad_2D* quad)
{
    fwrite(quad, sizeof(Quad_2D), 1, fpout); //write all quads
    cout<<"write, level: "<<quad->level<<", index: "<<quad->index<<endl;
    if(quad->isHasChildren)
    {
        for (int i = 0; i < m_num_children; i++)
        {
            writeForest(fpout, quad->children[i]);
        }
    }else{
        fwrite(quad->user_data, sizeof(FIELD_DATA_2D), 1, fpout);
    }
}

void readForest(FILE* fpin, Quad_2D* quad)
{
    fread(quad, sizeof(Quad_2D), 1, fpin); //write all quads
    cout<<"read, level: "<<quad->level<<", index: "<<quad->index<<", child: "<<quad->isHasChildren<<endl;

    if(quad->isHasChildren)
    {
        for (int i = 0; i < m_num_children; i++)
        {
            quad->children[i] = new Quad_2D;
            readForest(fpin, quad->children[i]);
        }
    }else
    {
        quad->user_data = new FIELD_DATA_2D;
        fread(quad->user_data, sizeof(FIELD_DATA_2D), 1, fpin);
        return;
    }
}

void write_binary(H2ONaCl::cH2ONaCl& eos, double TP_min[2], double TP_max[2], int max_level)
{
    const int dim = 2;
    // write to binary file
    cout<<eos.m_lut_PTX->get_root()->xyz[0]<<endl;
    string filename = "lut_PTX.bin";
    FILE* fpout = NULL;
    fpout = fopen(filename.c_str(), "wb");
    if(fpout == NULL)ERROR("Open file failed: "+filename);
    fwrite(&dim, sizeof(dim), 1, fpout);
    fwrite(TP_min, sizeof(double), 2, fpout);
    fwrite(TP_max, sizeof(double), 2, fpout);
    fwrite(&max_level, sizeof(int), 1, fpout);
    // write root
    // fwrite(eos.m_lut_PTX->get_root(), sizeof(LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim>>), 1, fpout);
    writeForest(fpout, eos.m_lut_PTX->get_root());
    fclose(fpout);

    // test read 
    FILE* fpin = NULL;
    fpin = fopen(filename.c_str(), "rb");
    if(!fpin)ERROR("Open file failed: "+filename);
    int dim0;
    fread(&dim0, sizeof(dim0), 1, fpin);
    cout<<"dim: "<<dim0<<endl;
    if(dim0==2)
    {
        const int dim = 2;
        double xy_min[dim], xy_max[dim];
        int max_level;
        Quad_2D root;
        fread(xy_min, sizeof(double), dim, fpin);
        fread(xy_max, sizeof(double), dim, fpin);
        fread(&max_level, sizeof(int), 1, fpin);
        LOOKUPTABLE_FOREST::LookUpTableForest<dim, LOOKUPTABLE_FOREST::FIELD_DATA<dim> > forest_in(xy_min, xy_max, max_level, sizeof(LOOKUPTABLE_FOREST::FIELD_DATA<dim>));

        Quad_2D quad;
        readForest(fpin, &quad);
        // for(int j=0;j<8;j++)
        // {
        //     Quad_2D quad;
        //     fread(&quad, sizeof(Quad_2D), 1, fpin);
        //     // readForest(fpin, forest_in.get_root());
        //     cout<<"level: "<<quad.level
        //         <<", index: "<<quad.index
        //         <<", children: "<<quad.children[0]
        //         <<endl;
        // }
        
    }else if(dim0==3)
    {

    }else
    {
        ERROR("dim in the file is neither 2 nor 3");
    }
    
    // cout<<"isHaveChild: "<<root.isHasChildren<<endl;
    // cout<<"xyz: "<<root.xyz[0]<<", "<<root.xyz[1]<<endl;
    fclose(fpin);

}
int main()
{
    int ind = 0;
    std::string dummy;

    clock_t start = clock();
    H2ONaCl::cH2ONaCl eos;
    const int dim =2;
    double TP_min[2] = {1 + 273.15, 5E5}; //T [K], P[Pa]
    double TP_max[2] = {700 + 273.15, 400E5};
    double X_wt = 0.2; //wt% NaCl [0,1]
    int min_level = 1;
    int max_level = 2;

    eos.createLUT_2D_PTX("constX", TP_min, TP_max, X_wt, min_level, max_level, "lut_PTX.vtu");
    
    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E4;
    for (size_t i = 0; i < n_randSample; i++)
    {
        double T_K = (rand()/(double)RAND_MAX)*(TP_max[0] - TP_min[0]) + TP_min[0];
        double p_Pa = (rand()/(double)RAND_MAX)*(TP_max[1] - TP_min[1]) + TP_min[1];
        // int ind_targetLeaf = forest.searchQuadrant(T_C, p_bar, 0);
        eos.m_lut_PTX->searchQuadrant(targetLeaf, T_K, p_Pa, X_wt);
        // cout<<"T_C: "<<T_K - 273.15<<", p_bar: "<<p_Pa/1E5<<", index: "<<targetLeaf->index<<endl;
        // cout<<eos.getPhaseRegionName(targetLeaf->user_data->phaseRegion_cell)<<endl;
        // compare to directally calculation
        H2ONaCl::PhaseRegion phaseRegion_cal = eos.findPhaseRegion(T_K - 273.15, p_Pa/1E5, eos.Wt2Mol(X_wt));
        // if(targetLeaf->user_data->phaseRegion_cell != phaseRegion_cal)
        // {
        //     cout<<"Need refine point "<<ind++<<": "<<targetLeaf->user_data->need_refine<<endl;
        // }
    }
    // // forest.searchQuadrant(targetLeaf, 200, 300, 0);
    // // cout<<targetLeaf->level<<endl;
    STATUS_time("Searching done", clock() - start);

    // ======= test write to binary =====
    // write_binary(eos, TP_min, TP_max, max_level);
    // call member of forest class
    eos.m_lut_PTX->write_to_binary("lut_TPX.bin");

    // destroy by hand
    // eos.destroyLUT_2D_PTX();

    // std::cout << "Enter to continue..." << std::endl;
    // std::getline(std::cin, dummy);
}