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
            for (size_t i = 0; i < x.size(); i++)
            {
                prop_lookup = sw.lookup_only(x[i], y[i]);
                pb_lookup[i] = prop_lookup.Region;
                if(prop_lookup.Region == H2ONaCl::MixPhaseRegion)ind++;
            }
        }
        break;
    case 3:
        {
            for (size_t i = 0; i < x.size(); i++)
            {
                prop_lookup = sw.lookup_only(x[i], y[i], z[i]);
                pb_lookup[i] = prop_lookup.Region;
                if(prop_lookup.Region == H2ONaCl::MixPhaseRegion)ind++;
            }
        }
        break;
    default:
        break;
    }
    printf("All %d/%ld (%.2f %%) points close to phase boundary.\n", ind, x.size(), ind/(double)x.size()*100);
    STATUS_time("Searching done", clock() - start);
    // write result to file
    STATUS("Writting results to file");
    ofstream fout("lookup_result.csv");
    if(!fout)ERROR("Open file failed: lookup_result.csv");
    fout<<"x\ty\tz\tphase region index\tphase region name"<<endl;
    for (size_t i = 0; i < x.size(); i++)
    {
        fout<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<"\t"<<pb_lookup[i]<<"\t"<<sw.m_phaseRegion_name[pb_lookup[i]]<<endl;
    }
    
    fout.close();

    return 0;
}