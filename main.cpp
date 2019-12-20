#include<iostream>

using namespace std;

#include"H2ONaCl.h"

int main()
{
    cout.precision(8);//control cout precision of float
    cH2ONaCl eos(1e7,100,0.03);
    eos.Calculate();
    return 0;
}