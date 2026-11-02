#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/*
 * A LAPW1 MPI rank owns a complete dense eigenproblem.  Allowing every rank
 * to start a multithreaded BLAS oversubscribes both CPUs and memory badly
 * (and can make macOS kill a rank).  Apply conservative defaults for an MPI
 * run.  Users who intentionally want nested MPI+BLAS parallelism can set
 * EDMFT_BLAS_THREADS to the desired positive thread count.
 *
 * This is called after MPI_Init and before LAPW1 makes its first BLAS call.
 */
void limit_blas_threads_(int *nprocs)
{
      const char *requested;
      const char *threads;

      if (nprocs == NULL || *nprocs <= 1)
            return;

      requested = getenv("EDMFT_BLAS_THREADS");
      threads = (requested != NULL && requested[0] != '\0') ? requested : "1";

      setenv("OMP_NUM_THREADS", threads, 1);
      setenv("OPENBLAS_NUM_THREADS", threads, 1);
      setenv("MKL_NUM_THREADS", threads, 1);
      setenv("BLIS_NUM_THREADS", threads, 1);
      setenv("VECLIB_MAXIMUM_THREADS", threads, 1);
}

void cputim_(double *dsec)
{
      struct tms buffer;

      times(&buffer);
      *dsec = (double)buffer.tms_utime/sysconf(_SC_CLK_TCK);
/*      fprintf(stderr,"inside cputim: %lf\n",dsec);*/
}


void walltim_(double *dsec)
{       struct timeval tv;
        gettimeofday(&tv, NULL);
        *dsec = (double) tv.tv_sec + tv.tv_usec/1000000.0;
}

void cputim(double *dsec)
{ cputim_(dsec);
}

void walltim(double *dsec)
{ walltim_(dsec);
}

