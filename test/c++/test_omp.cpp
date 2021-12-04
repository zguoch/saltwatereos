#include <stdio.h>
#include <omp.h>
#include <unistd.h>

int fib(int n)
{
    // sleep(1);
  int i, j;
  if (n<2)
    {
        return n;
    }
  else
    {
       #pragma omp task shared(i) firstprivate(n)
       i=fib(n-1);
    printf("1 on thread %d \n", omp_get_thread_num());

       #pragma omp task shared(j) firstprivate(n)
       j=fib(n-2);
    printf("2 on thread %d \n", omp_get_thread_num());

       #pragma omp taskwait
    printf("Current int %d is on thread %d \n", i + j, omp_get_thread_num());
       return i+j;
    }
}

int main()
{
  int n = 10;

  #pragma omp parallel shared(n)
  {
    #pragma omp single
    {
        printf("%d\n", omp_get_num_threads());
        printf ("fib(%d) = %d\n", n, fib(n));
    }
  }
}