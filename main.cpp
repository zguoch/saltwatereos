#include<iostream>

using namespace std;

#include"H2ONaCl.h"

int main()
{
    cH2ONaCl eos(1e7,100,0.03);
    eos.Calculate();
    return 0;
}