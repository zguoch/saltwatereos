#include "H2ONaCl.H"

int main(int argc, char** argv)
{
    if(argc!=3)ERROR("Usage: lutLookup myLUT.bin xyz.txt");

    // 1. load lut
    H2ONaCl::cH2ONaCl sw;
    sw.loadLUT(argv[1]);
    // 2. read xyz 
    STATUS("Reading xyz data ...");
    vector<double> x, y, z;
    double tmp_x, tmp_y, tmp_z;
    ifstream fin(argv[2]);
    if(!fin)ERROR("Open file failed: "+string(argv[2]));
    while (!fin.eof())
    {
        fin>>tmp_x>>tmp_y>>tmp_z;
        x.push_back(tmp_x);
        y.push_back(tmp_y);
        z.push_back(tmp_z);
    }
    fin.close();
    // 3. lookup
    STATUS("Start looking up ... ");
    H2ONaCl::PROP_H2ONaCl prop_lookup;
    vector<H2ONaCl::PhaseRegion> pb_lookup(x.size(), H2ONaCl::SinglePhase_L);
    int ind = 0;
    time_t start = clock();
    switch (sw.m_dim_lut)
    {
    case 2:
        {
            const int dim = 2;
            H2ONaCl::LookUpTableForest_2D* pLUT = (H2ONaCl::LookUpTableForest_2D*)sw.m_pLUT;
            LOOKUPTABLE_FOREST::Quadrant<dim,H2ONaCl::FIELD_DATA<dim> > *targetLeaf = NULL;
            FILE* fp = NULL;
            fp = fopen("lookup_result.csv", "w");
            if(!fp)ERROR("Open file failed: lookup_result.txt");
            double* props = new double[pLUT->m_map_props.size()];
            fprintf(fp, "Phase region\tNeed_refine");
            for (auto &m : pLUT->m_map_props)
            {
                fprintf(fp, "\t%s", m.second.longName);
            }
            fprintf(fp, "\n");
            
            for (size_t i = 0; i < x.size(); i++)
            {
                targetLeaf = sw.lookup(props, x[i], y[i]);
                fprintf(fp, "%d\t%d", targetLeaf->qData.leaf->user_data->phaseRegion_cell, targetLeaf->qData.leaf->user_data->need_refine);
                for (size_t j = 0; j < pLUT->m_map_props.size(); j++)
                {
                    fprintf(fp, "\t%f", props[j]);
                }
                fprintf(fp, "\n");
                if(targetLeaf->qData.leaf->user_data->phaseRegion_cell == H2ONaCl::MixPhaseRegion)ind++;
            }
            fclose(fp);
            delete[] props;
        }
        break;
    case 3:
        {
            const int dim = 3;
            H2ONaCl::LookUpTableForest_3D* pLUT = (H2ONaCl::LookUpTableForest_3D*)sw.m_pLUT;
            LOOKUPTABLE_FOREST::Quadrant<dim,H2ONaCl::FIELD_DATA<dim> > *targetLeaf = NULL;
            FILE* fp = NULL;
            fp = fopen("lookup_result.csv", "w");
            if(!fp)ERROR("Open file failed: lookup_result.txt");
            double* props = new double[pLUT->m_map_props.size()];
            fprintf(fp, "Phase region\tNeed_refine");
            for (auto &m : pLUT->m_map_props)
            {
                fprintf(fp, "\t%s", m.second.longName);
            }
            fprintf(fp, "\n");

            for (size_t i = 0; i < x.size(); i++)
            {
                targetLeaf = sw.lookup(props, x[i], y[i], z[i]);
                fprintf(fp, "%d\t%d", targetLeaf->qData.leaf->user_data->phaseRegion_cell, targetLeaf->qData.leaf->user_data->need_refine);
                for (size_t j = 0; j < pLUT->m_map_props.size(); j++)
                {
                    fprintf(fp, "\t%f", props[j]);
                }
                fprintf(fp, "\n");
                if(targetLeaf->qData.leaf->user_data->phaseRegion_cell == H2ONaCl::MixPhaseRegion)ind++;
            }
            fclose(fp);
            delete[] props;
        }
        break;
    default:
        break;
    }
    printf("All %d/%ld (%.2f %%) points close to phase boundary.\n", ind, x.size(), ind/(double)x.size()*100);
    STATUS_time("Searching done", clock() - start);
    // // write result to file
    // STATUS("Writting results to file");
    // ofstream fout("lookup_result.csv");
    // if(!fout)ERROR("Open file failed: lookup_result.csv");
    // fout<<"x\ty\tz\tphase region index\tphase region name"<<endl;
    // for (size_t i = 0; i < x.size(); i++)
    // {
    //     fout<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<"\t"<<pb_lookup[i]<<"\t"<<sw.m_phaseRegion_name[pb_lookup[i]]<<endl;
    // }
    
    // fout.close();

    return 0;
}