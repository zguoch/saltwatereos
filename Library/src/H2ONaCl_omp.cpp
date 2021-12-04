// because swig can not pass the OpenMP, so make all omp-related code in this cpp file. and do not compile this file in the swig CMakeLists.txt

#include "H2ONaCl.H"
namespace H2ONaCl
{
    void cH2ONaCl::set_num_threads(int num_threads)
    {
        omp_set_num_threads(num_threads);
        m_num_threads = num_threads;
    };
    int cH2ONaCl::get_num_threads()
    {
        return omp_get_num_threads();
    };
}