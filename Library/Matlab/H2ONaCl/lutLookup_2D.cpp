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
// this function just used to test and debug, please ignore it.
void load_binary_2d(string filename)
{
    int ind = 0;
    H2ONaCl::cH2ONaCl eos;
    H2ONaCl::PROP_H2ONaCl prop_cal;
    clock_t start = clock();
    const int dim = 2;
    int dim_file = LOOKUPTABLE_FOREST::get_dim_from_binary(filename);
    if(dim!=dim_file)ERROR("The input LUT file is not 2D.");
    eos.loadLUT(filename);
    H2ONaCl::LookUpTableForest_2D* pLUT = (H2ONaCl::LookUpTableForest_2D*)eos.m_pLUT;
    // eos.save_lut_to_vtk("lut_load_save.vtu");
    double Xmin = pLUT->m_xyz_min[0], Xmax = pLUT->m_xyz_max[0], Ymin = pLUT->m_xyz_min[1], Ymax = pLUT->m_xyz_max[1];
    double constZ = pLUT->m_constZ;
    int index_rho = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_rho));
    int index_h = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_h));
    int index_T = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_T));

    STATUS("Start search ... ");
    start = clock();
    LOOKUPTABLE_FOREST::Quadrant<dim,H2ONaCl::FIELD_DATA<dim> > *targetLeaf = NULL;
    int n_randSample = 1E3;
    double* props = new double[pLUT->m_map_props.size()];
    int ind_rho = 0;
    int nx = 100, ny = 150;
    double dx = (Xmax - Xmin)/(nx - 1);
    double dy = (Ymax - Ymin)/(ny - 1);
    double x, y;
    ofstream fpout_xx("XX.txt");
    ofstream fpout_yy("YY.txt");
    ofstream fpout_zz("ZZ.txt");
    ofstream fpout_zz_lut("ZZ_lut.txt");
    bool isCal = true;
    for (int i = 0; i < ny; i++)
    {
        y = Ymin + i*dy;
        for (int j = 0; j < nx; j++)
        {
            x = Xmin + j*dx;
            if(isCal)
            {
                switch (pLUT->m_TorH)
                {
                case LOOKUPTABLE_FOREST::EOS_ENERGY_T:
                    {
                        switch (pLUT->m_const_which_var)
                        {
                        case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
                            prop_cal = eos.prop_pTX(y, x, constZ);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
                            prop_cal = eos.prop_pTX(constZ, y, x);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
                            prop_cal = eos.prop_pTX(y, constZ, x);
                            break;
                        default:
                            break;
                        }
                    }
                    break;
                case LOOKUPTABLE_FOREST::EOS_ENERGY_H:
                    {
                        switch (pLUT->m_const_which_var)
                        {
                        case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
                            prop_cal = eos.prop_pHX(y, x, constZ);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
                            prop_cal = eos.prop_pHX(constZ, y, x);
                            break;
                        case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
                            prop_cal = eos.prop_pHX(y, constZ, x);
                            break;
                        default:
                            break;
                        }
                    }
                    break;
                default:
                    break;
                }
                fpout_zz<<prop_cal.Rho<<" ";
            }
            
            targetLeaf = eos.lookup(props, x, y, false);
            fpout_xx<<x<<" ";
            fpout_yy<<y<<" ";
            fpout_zz_lut<<props[index_rho]<<" ";
        }
        fpout_xx<<endl;
        fpout_yy<<endl;
        fpout_zz<<endl;
        fpout_zz_lut<<endl;
    }
    fpout_xx.close();
    fpout_yy.close();
    fpout_zz.close();
    fpout_zz_lut.close();
    delete[] props;

    STATUS_time("Searching done", clock() - start);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check arguments
    if(nrhs!=3){
    mexErrMsgTxt("Usage: lookupLUT_2D(filename, x, y)\nNote that the x, y order MUST BE consistent with the lookup table file.");
    }
    /* allocate input variables */
    double *xMat, *yMat;
    mwSize  m, n;
    /* allocate output variables */
    double *Rho_out, *T_out;
    LOOKUPTABLE_FOREST::NeedRefine *needRefine_out;

    string filename = mxArrayToString(prhs[0]);
    mGetMatrix(prhs[1], &xMat, "px", &m, &n);
    mGetMatrix(prhs[2], &yMat, "py", &m, &n);

    plhs[0] = mxCreateNumericMatrix(m, n, mxINT32_CLASS, mxREAL);
    needRefine_out = (LOOKUPTABLE_FOREST::NeedRefine*)mxGetData(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
    Rho_out = (double*)mxGetData(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);
    T_out = (double*)mxGetData(plhs[2]);
    // maybe make a safety check, input matrix dimension consistency
    clock_t start = clock();
    STATUS("Passing data end, start loading data ...");
    H2ONaCl::cH2ONaCl sw;
    // load_binary_2d(filename);

    sw.loadLUT(filename); //\Todo pass this sw pointer as input arg from outside, don't always load the LUT file.
    STATUS_time("Load LUT end. ", clock()-start);
    start = clock();
    STATUS("Start looking up ...");
    const int dim = 2;
    H2ONaCl::LookUpTableForest_2D* pLUT = (H2ONaCl::LookUpTableForest_2D*)sw.m_pLUT;
    LOOKUPTABLE_FOREST::Quadrant<dim,H2ONaCl::FIELD_DATA<dim> > *targetLeaf = NULL;
    double* props = new double[pLUT->m_map_props.size()];
    // ====== get order of props in the data array =====
    int index_Rho = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_rho));
    int index_T = distance(pLUT->m_map_props.begin(),pLUT->m_map_props.find(Update_prop_T));
    // =================================================
    // ofstream fpout_x("XX.txt");
    // ofstream fpout_y("YY.txt");
    // ofstream fpout_z("ZZ.txt");
    // ofstream fpout_prop0("prop0.txt");
    // ofstream fpout_prop1("prop1.txt");
    // ofstream fpout_prop2("prop2.txt");
    int index = 0;
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
      {
        index = i*n +j;
        targetLeaf = sw.lookup(props, xMat[index], yMat[index]); 
        needRefine_out[index] = targetLeaf->user_data->need_refine;
        T_out[index] = props[index_T];
        Rho_out[index] = props[index_Rho]; 
        // fpout_x<<xMat[index]<<" ";
        // fpout_y<<yMat[index]<<" ";
        // fpout_z<<yMat[index]<<" ";
        // fpout_prop0<<props[0]<<" ";
        // fpout_prop1<<props[1]<<" ";
        // fpout_prop2<<props[2]<<" ";
      }
    //   fpout_x<<endl;
    //   fpout_y<<endl;
    //   fpout_z<<endl;
    //   fpout_prop0<<endl;
    //   fpout_prop1<<endl;
    //   fpout_prop2<<endl;
    }
    // fpout_x.close();
    // fpout_y.close();
    // fpout_z.close();
    // fpout_prop0.close();
    // fpout_prop1.close();
    // fpout_prop2.close();


    delete[] props;
    STATUS_time("Lookup end. ", clock()-start);
}
