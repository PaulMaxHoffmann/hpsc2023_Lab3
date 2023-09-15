#include<omp.h>     // Required for OpenMP
#include<mpi.h>     // Required for MPI
#include<stdio.h>
#include<chrono>
#include<ctime>

#include "main.h"
// ----------------------------------------------------------------
//
// Notes
//
// (1) How to determine what hardware you currently have in Ubuntu:
//
//     cat /proc/cpuinfo | grep -e processor -e cores
//
// (2) How to compile with  OpenMP only:
//
//     g++ -fopenmp main.cpp -o main
//
// (3) How to compile with MPI and OpenMP:
//
//     mpicxx -fopenmp main.cpp -o main
//
// ----------------------------------------------------------------


// ----------------------------------------------------------------
// 
//   MPI Lesson Plan:
//
//   (1) Modifying the Makefile
//   (2) Running a parallel code with mpirun
//   (3) Initializing MPI
//   (4) MPI commands...add one at a time
//   (5) OpenMP command on loop
// 
// 
// ----------------------------------------------------------------

void printArray(string arrayName, int n , int *array , int myPE)
{
  for ( int i = 0 ; i < n ; ++i ) cout << "myPE: " << myPE << " " << arrayName << "[" << i << "] = " << array[i] << endl;
}

void printVal(string valName, int val , int myPE)
{
  cout << "myPE: " << myPE << " " << valName << " = " << val << endl;
}


int main(  int argc, char *argv[] )
{

  // (2)
  
  printf("Hello World!\n");

  // (3)

  MPI_Init     (&argc         , &argv       );

  int numPE = 4;                              // Must match mpirun command
  int myPE;
  MPI_Comm_size(MPI_COMM_WORLD,  &numPE);     // Set the number of processors
  MPI_Comm_rank(MPI_COMM_WORLD, &myPE  );     // Get my PE number


  cout << "myPE: " << myPE << " Hello, world!" << endl;
  
  // (4)   Broadcast

  int n = 10;
  int a[n];

  if ( myPE == 0 ) for ( int i = 0 ; i < 10; ++i ) a[i] = 333;
  
  MPI_Bcast( a , n , MPI_INT, 0, MPI_COMM_WORLD);

  printArray("BCast a",n,a,myPE);


  // (4)   Scatter

  int m = 40;
  int b[m];

  if ( myPE == 0 ) for ( int i = 0 ; i < 40; ++i ) b[i] = i;
  
  MPI_Scatter( b , n , MPI_INT,
	       a , n , MPI_INT, 0 , MPI_COMM_WORLD);

  printArray("Scatter a",n,a,myPE);  // Results good on all PE

  // (4)  Gather


  for ( int i = 0 ; i < n ; ++i ) a[i] = myPE;
  
  MPI_Gather( a , n , MPI_INT,
	      b , n , MPI_INT, 0 , MPI_COMM_WORLD);

  printArray("Gather b",m,b,myPE);  // Results only valid on PE 0

  //  (4) Reduce

  int myVal = myPE;

  int myReducedVal;

  MPI_Reduce(&myVal, &myReducedVal, 1, MPI_INT, MPI_SUM, 0 , MPI_COMM_WORLD);

  printVal("SummedValue",myReducedVal,myPE);  // Result only valid on PE 0

  // Also can do MPI_MAX, MPI_MIN, MPI_PROD...

  //  (4) All-Reduce

  MPI_Allreduce(&myVal, &myReducedVal, 1 , MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
  printVal("SummedValue",myReducedVal,myPE);  // Result valid on all PE

  //  (5) All-Gather

  int myContribution[n];
  for ( int i = 0 ; i < n ; ++i ) myContribution[i] = myPE;
  int everyonesContribution[m];

  printArray("myContribution",n,myContribution,myPE);  // Result valid on all PE
  
  MPI_Allgather( myContribution        ,  n , MPI_INT,
                 everyonesContribution ,  n , MPI_INT, MPI_COMM_WORLD );
  
  printArray("everyonesContribution",m,everyonesContribution,myPE);  // Result valid on all PE

  

  // ~LabBegin~ "LabStrip~
  
  // Part 2: OpenMP ---- Threading
  
#pragma omp parallel
  {
    cout << "myPE: " << myPE << " Hello, world from thread " << omp_get_thread_num() << endl;;
  }


  cout << "myPE: " << myPE << " Num threads = " << omp_get_num_threads() << endl;
  cout << "myPE: " << myPE << " Max threads = " << omp_get_max_threads() << endl;

  // ~LabEnd~
  
  

  MPI_Barrier(MPI_COMM_WORLD);  
  MPI_Finalize();

  
  return 0;



  
}



