// /* File : example.i */
// %module example
// %{
// /* Put headers and other declarations here */
// extern double My_variable;
// extern int    fact(int);
// extern int    my_mod(int n, int m);
// %}

// extern double My_variable;
// extern int    fact(int);
// extern int    my_mod(int n, int m);

/* File: example.i */
%module example

%{
#define SWIG_FILE_WITH_INIT
#include "example.h"
%}

int fact(int n);
#define aaa 100
// Very simple C++ example for linked list
class Foo {
public:
      Foo(){};
     ~Foo(){};
      int  bar(int x);
      int  x;
};