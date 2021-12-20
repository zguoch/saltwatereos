#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;
#include "H2ONaCl.H"

#define MAX(a, b) ((a)>(b) ? (a):(b))
#define MIN(a, b) ((a)<(b) ? (a):(b))

void mGetVector(const mxArray *prhs, double **out, const char *varname, mwSize *out_m)
{
  mwSize m, n;
  double *temp;

  m = mxGetM(prhs);
  n = mxGetN(prhs);

  if(!(m==1 || n==1)){
    char buff[256];
    sprintf(buff, "'%s' must be a vector.\n", varname);
    mexErrMsgTxt(buff);
  }

  if(!mxIsDouble(prhs)){    char buff[256];
    sprintf(buff, "'%s' must be of type 'double'.\n", varname);
    mexErrMsgTxt(buff);
  }

  temp = (double*)mxGetData(prhs);
  *out = temp;
  *out_m = MAX(m,n);
}

void mGetMatrix(const mxArray *prhs, double **out, const char *varname, mwSize *out_m, mwSize *out_n)
{
  mwSize m, n;
  double *temp;

  m = mxGetM(prhs);
  n = mxGetN(prhs);

  if(!mxIsDouble(prhs)){
    char buff[256];
    sprintf(buff, "'%s' must be of type 'double'.\n", varname);
    mexErrMsgTxt(buff);
  }

  temp   = (double*)mxGetData(prhs);
  *out   = temp;
  *out_m = m;
  *out_n = n;
}

void load_binary_3d(string filename)
{
    int ind = 0;
    H2ONaCl::cH2ONaCl eos;
    H2ONaCl::PROP_H2ONaCl prop_cal, prop_lookup;
    clock_t start = clock();
    const int dim = 3;
    eos.loadLUT(filename);
    cout<<"dim of the bin file: "<<eos.m_dim_lut<<endl;
    // eos.save_lut_to_vtk("lut_TPX.vtu");
    H2ONaCl::LookUpTableForest_3D* pLUT = (H2ONaCl::LookUpTableForest_3D*)eos.m_pLUT;

    double TorHmin = pLUT->m_xyz_min[0], TorHmax = pLUT->m_xyz_max[0], Pmin = pLUT->m_xyz_min[1], Pmax = pLUT->m_xyz_max[1];
    double Xmin = pLUT->m_xyz_min[2], Xmax = pLUT->m_xyz_max[2];
    double constX = 0.2; //lookup x=0.2 slice
    double constP = 200E5;
    int index_rho = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_rho));
    int index_h = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_h));
    int index_T = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_T));
    int nP = 100, nTorH = 150, nX = 100;
    double dTorH = (TorHmax - TorHmin)/(nTorH - 1);
    double dP = (Pmax - Pmin)/(nP - 1);
    double dX = (Xmax - Xmin)/(nX - 1);
    double TorH, P, X;
    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,H2ONaCl::FIELD_DATA<dim> > *targetLeaf = NULL;
    double* props = new double[pLUT->m_map_props.size()];
    int ind_rho = 0;
    ofstream fpout_xx("XX.txt");
    ofstream fpout_yy("YY.txt");
    ofstream fpout_zz("ZZ.txt");
    ofstream fpout_zz_lut("ZZ_lut.txt");
    ofstream fpout_phase("phase.txt");
    ofstream fpout_needRefine("refine.txt");
    ofstream fpout_T_lut("T_lut.txt");
    bool isCal = true;
    for (int i = 0; i < nTorH; i++)
    {
        TorH = TorHmin + i*dTorH;
        for (int j = 0; j < nX; j++)
        {
            X = Xmin + j*dX;
            if(isCal)
            {
                switch (pLUT->m_TorH)
                {
                case LOOKUPTABLE_FOREST::EOS_ENERGY_T:
                    prop_cal = eos.prop_pTX(constP, TorH, X);
                    break;
                case LOOKUPTABLE_FOREST::EOS_ENERGY_H:
                    prop_cal = eos.prop_pHX(constP, TorH, X);
                    break;
                default:
                    break;
                }
                fpout_zz<<prop_cal.Rho<<" ";
            }
            
            targetLeaf = eos.lookup(props, TorH, constP, X);
            fpout_xx<<X<<" ";
            fpout_yy<<TorH<<" ";
            fpout_zz_lut<<props[index_rho]<<" ";
            fpout_T_lut<<props[index_T]<<" ";
            fpout_phase<<targetLeaf->user_data->phaseRegion_cell<<" ";
            fpout_needRefine<<targetLeaf->user_data->need_refine<<" ";
        }
        fpout_xx<<endl;
        fpout_yy<<endl;
        fpout_zz<<endl;
        fpout_zz_lut<<endl;
        fpout_needRefine<<endl;
        fpout_phase<<endl;
        fpout_T_lut<<endl;
    }
    fpout_xx.close();
    fpout_yy.close();
    fpout_zz.close();
    fpout_zz_lut.close();
    fpout_phase.close();
    fpout_needRefine.close();
    fpout_T_lut.close();
    delete[] props;

    STATUS_time("Searching done", clock() - start);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check arguments
    if(nrhs!=4){
    mexErrMsgTxt("Usage: lookupLUT_3D(filename, x, y, z)\nNote that the x, y, z order MUST BE T/H, P, X.");
    }
    /* allocate input variables */
    double *xMat, *yMat, *zMat;
    mwSize  m, n;
    /* allocate output variables */
    double *Rho_out, *T_out;
    LOOKUPTABLE_FOREST::NeedRefine *needRefine_out;
    H2ONaCl::PhaseRegion *PhaseRegion;

    string filename = mxArrayToString(prhs[0]);
    mGetMatrix(prhs[1], &xMat, "px", &m, &n);
    mGetMatrix(prhs[2], &yMat, "py", &m, &n);
    mGetMatrix(prhs[3], &zMat, "pz", &m, &n);

    plhs[0] = mxCreateNumericMatrix(m, n, mxINT32_CLASS, mxREAL);
    needRefine_out = (LOOKUPTABLE_FOREST::NeedRefine*)mxGetData(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
    Rho_out = (double*)mxGetData(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);
    T_out = (double*)mxGetData(plhs[2]);
    // plhs[2] = mxCreateNumericMatrix(m, n, mxINT32_CLASS, mxREAL);
    // PhaseRegion = (H2ONaCl::PhaseRegion*)mxGetData(plhs[0]);
    // test
    // load_binary_3d(filename);

    // maybe make a safety check, input matrix dimension consistency
    clock_t start = clock();
    STATUS("Passing data end, start loading data ...");
    H2ONaCl::cH2ONaCl sw;
    sw.loadLUT(filename); //\Todo pass this sw pointer as input arg from outside, don't always load the LUT file.

    H2ONaCl::LookUpTableForest_3D* pLUT = (H2ONaCl::LookUpTableForest_3D*)sw.m_pLUT;
    double* props = new double[pLUT->m_map_props.size()];
    // ====== get order of props in the data array =====
    int index_Rho = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_rho));
    int index_T = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_T));
    // =================================================
    STATUS_time("Load LUT end. ", clock()-start);
    start = clock();
    STATUS("Start looking up ...");
    const int dim = 3;
    LOOKUPTABLE_FOREST::Quadrant<dim,H2ONaCl::FIELD_DATA<dim> > *targetLeaf = NULL;
    int index = 0;
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
      {
        index = i*n +j;
        targetLeaf = sw.lookup(props, xMat[index], yMat[index], zMat[index]); 
        needRefine_out[index] = targetLeaf->user_data->need_refine;
        T_out[index] = props[index_T];
        Rho_out[index] = props[index_Rho];
        // PhaseRegion[index] = targetLeaf->user_data->phaseRegion_cell;
      }
    }
    delete[] props;
    STATUS_time("Lookup end. ", clock()-start);
}
