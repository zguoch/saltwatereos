#include "H2ONaCl.H"
#include <iostream>
H2ONaCl::cH2ONaCl eos;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    std::cout<<eos.m_NaCl.T_Melting(200)<<std::endl;
}