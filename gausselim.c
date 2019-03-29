/* * *
 * Hriddhi Kulkarni
 * * * *
 * hriddhi.kulkarni@ttu.edu
 * * * */

/* Header files  */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <limits.h>
#include <pthread.h>
#include <unistd.h>

#define _ENTRANT
#define minimum(a, b)            ((a) < (b)) ? (a) :  (b)

#define NMAX 3000 /*  We define the Maximum value of N (matrix dimension size). */
#define PMAX 10 /*  We define the Maximum matrix size for printing. */

int M; /* Let the Size of the Matrix be M. */
int numth; /* Let the Number of threads to use be numth */
volatile float A[NMAX][NMAX], B[NMAX], X[NMAX]; /*  The values for A * X = B are given and we have to solve for X.  */
int norm, cRow, flag; 		/* Here cRow is the current row, flag is the counter variable*/
void threadc(void); 		/* This function spawns off <NumThreads>. */
void *gausselim(void *); 	/* This function runs concurrently in <NumThreads> threads. */
int rchunkSize(int *, int ); 	/* The rows chunk size is determined by this function. */
void barriersync(int * ); 		/* We get the barrier with synchronization using this function*/
void threadWait(void); /* This function awaits the termination of all threads.        */
unsigned int seedRet(void); /* Returns a seed for <srand()> based on the time.           */
void param(int , char **); /* The program parameters are set from the command-line arguments.     */
void init(void); /* We initialize the matrices A, B and X and set X to 0.0s.                 */

/* These print the input matrix and the solution matrix  */
void inputm(void);
void outputm(void);

/*  We use this to store the identity of each thread. */
pthread_t threads[_POSIX_THREAD_THREADS_MAX];

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t fLock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t fNext = PTHREAD_COND_INITIALIZER;

/* Main function starts here*/
void main(int ArgC, char **ArgV)
{
/* We get the elapsed times using <gettimeofday()>. */
 struct timeval eStart, eStop;
 struct timezone Dummytz;
 clock_t eStartt, eStoptt;
                   
 unsigned long uStart, uStop;  /* The elapsed times using <times()>. */

 struct tms cstart, cstop; /* This gives the CPU start and stop times for the threads. */

 int row, column;
 float CLK_TCK;

 param(ArgC, ArgV);
 init();
 inputm();

 cRow = norm+1;
 flag = numth-1;

 printf("Starting the clock ....\n");
 gettimeofday(&eStart, &Dummytz);
 eStartt = times(&cstart);

 threadc();

 threadWait();

/* Back substitution starts here */
 for (row = M-1; row >= 0; row--)
 {
  X[row] = B[row];
  for (column = M-1; column > row; column--)
  X[row] -= A[row][column]*X[column];
  X[row] /= A[row][row];
 }

 gettimeofday(&eStop, &Dummytz);
 eStoptt = times(&cstop);
 printf("The clock is stopped.\n");

 uStart = (unsigned long)eStart.tv_sec*1000000+eStart.tv_usec;
 uStop = (unsigned long)eStop.tv_sec*1000000+eStop.tv_usec;

 outputm();

 printf("The elapsed time = %g ms.\n",
(float)(uStop-uStart)/(float)1000);

 printf("The elapsed time according to <times()> = %g ms.\n",
(eStoptt-eStartt)/(float)CLK_TCK*1000);

 printf("The CPU times are accurate to the nearest %g ms.\n",
1.0/(float)CLK_TCK*1000.0);

 printf("The total CPU time required for the parent = %g ms.\n",
(float)((cstop.tms_utime+cstop.tms_stime)-
(cstart.tms_utime+cstart.tms_stime))/(float)CLK_TCK*1000);

 printf("The system CPU time required for the parent = %g ms.\n\n",
(float)(cstop.tms_stime-cstart.tms_stime)/(float)CLK_TCK*1000);
 }

unsigned int seedRet(void)
{
 struct timeval t;
 struct timezone Dummytz;

 gettimeofday(&t, &Dummytz);
 return (unsigned int)(t.tv_usec);
 }


void param(int ArgC, char **ArgV)
{

 int subparam = 0; /* We set <subparam> to 1 only if we use the submission parameters. */
 int seedVal = 0; /* The random seed value. */

char L_cuserid; /* Here <L_cuserid> is a macro defined in stdio.h */
char uid[62]; /* Here <uid> contains the user name.*/

srand(seedRet()); /* Here randomising is done. */

/* We now start reading the command-line arguments. */
if (ArgC != 3)
   {
    if (ArgC == 2 && !strcmp(ArgV[1], "subparam"))
       {
        subparam = 1;
        M = 4;
        numth = 2;
     //printf("Submission run for \"%s\".\n", cuserid(uid));
        srand(uid[0]);
        }
    else
       {
        if (ArgC == 4)
           {
            seedVal = atoi(ArgV[3]);
            srand(seedVal);
            printf("The random seed value = %d.\n", seedVal);
            }
        else
           {
            printf("The usage is: %s <matrix_dimension> <num_threads> \
[random seed]\n", ArgV[0]);
            printf("       %s subparam\n", ArgV[0]);
            exit(0);
			}
        }
    }

/*  Here we start interpreting the command-line arguments. */
if (!subparam)
   {
    M = atoi(ArgV[1]);
    if (M < 1 || M > NMAX)
       {
        printf("Matrix dimension M = %d is out of range.\n", M);
        exit(0);
        }

    numth = atoi(ArgV[2]);
    if (numth < 1)
       {
        printf("Warning: The invalid number of threads = %d.  Using 1.\n",
numth);
        numth = 1;
        }

        if (numth > _POSIX_THREAD_THREADS_MAX)
       {
        printf("Warning: %d threads requested; but only %d are available.\n",
numth, _POSIX_THREAD_THREADS_MAX);
        numth = _POSIX_THREAD_THREADS_MAX;
        }
    }

/* Here we are printing the parameters. */
printf("The dimension of matrix M = %d.\n", M);
printf("The number of threads = %d.\n", numth);
}


void init(void)
{
 int row, column;

 printf("\nInitializing ...\n");
 for (column = 0; column < M; column++)
 {
/* We add 1 is to ensure that there are non-zero entries in the coefficient matrix A. */
	for (row = 0; row < M; row++)
  	A[row][column] = (float)rand()/RAND_MAX+1;
  	B[column] = (float)rand()/RAND_MAX+1;
  	X[column] = 0.0;
  }
 }


void inputm(void)
{
 int row, column;

 if (M < PMAX)
    {
     printf("\nA =\n\t");
     for (row = 0; row < M; row++)
     {
      for (column = 0; column < M; column++)
      printf("%6.3f  %s", A[row][column], (column < M-1) ? ", " : ";\n\t");
      }
     printf("\nB = [");
     for (column = 0; column < M; column++)
     printf("%6.3f  %s", B[column], (column < M-1) ? "; " : "]\n");
     }
 }

 void outputm(void) 
 {
 	int row;
	if (M < PMAX)
    {
     	printf("\nX = [");
     	for (row = 0; row < M; row++)
     	printf("%6.3f  %s", X[row], (row < M-1) ? "; " : "]\n");
    }
 }


void *gausselim(void *dummy)
{
 int FRow = 0, row, column; /*  <FRow> denotes the first row of the chunk assigned to a thread.  */

 int NormRow = 0; /*   This is the Normalisation row. */

 float multi; /* <multi> is the multiplier*/

 int Sizeofchunk; /*Size of the chunk*/

/* Gaussian elimination begins here. */

 while (NormRow < M-1)
 {

  while (Sizeofchunk = rchunkSize(&FRow, NormRow)) /* The row chunk should be assigned to this thread. */
  {

/* We perform the eliminations across these rows concurrently.   */
   for (row = FRow; row < (minimum(M, FRow+Sizeofchunk)); row++)
   {
    multi = A[row][NormRow]/A[NormRow][NormRow];
    for (column = NormRow; column < M; column++)
    A[row][column] -= A[NormRow][column]*multi;
    B[row] -= B[NormRow]*multi;
    }
   }

  barriersync(&NormRow); /* We wait until all threads are done with this stage and then proceed to handle the next normalisation row.  */
  }
 }

/* Barrier with synchronization begins here */
 void barriersync(int *NormRow)
{

pthread_mutex_lock(&fLock); /* We implement synchronisation using condition variables.    */

 if (flag == 0)
    {
     norm++; /*  Here only the last thread for each value of <norm> reaches. */
     flag = numth-1;
     cRow = norm+1;
     pthread_cond_broadcast(&fNext);
     }
 else
    {
     flag--;
     pthread_cond_wait(&fNext, &fLock);
     }

 *NormRow = norm; /* <*NormRow> is each thread's view of the global variable <norm>.  */

 pthread_mutex_unlock(&fLock);
 }

 int rchunkSize(int *FRow, int NormRow)
{
 int Sizeofchunk;
 pthread_mutex_lock(&mutex);

 *FRow = cRow;

 Sizeofchunk = (*FRow < M) ? (M-NormRow-1)/(2*numth)+1 : 0; /*   We determine the chunk size here for guided-self scheduling.   */
 cRow += Sizeofchunk;

 pthread_mutex_unlock(&mutex);

 return Sizeofchunk;
 }


void threadc(void)
{
 int i;
 for (i = 0; i < numth; i++) {
 pthread_create(&threads[i], NULL, gausselim, NULL);
 }
}

void threadWait(void)
{
 int i;
 for (i = 0; i < numth; i++) {
 pthread_join(threads[i], NULL);
 }
}