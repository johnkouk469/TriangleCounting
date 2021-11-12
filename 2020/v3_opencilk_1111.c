/*
** v3.c -- version 3 of Vertexwise triangle counting
**
** OpenCilk implementation
**
** Representng graph in CSC format
** Loading sparse adjacency matrix using the Matrix Market format
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cilk/cilk.h>
#include <pthread.h>
#include "mmio.h"

void coo2csc(
    int       * const row,       /*!< CSC row start indices */
    int       * const col,       /*!< CSC column indices */
    int const * const row_coo,   /*!< COO row indices */
    int const * const col_coo,   /*!< COO column indices */
    int const         nnz,       /*!< Number of nonzero elements */
    int const         n,         /*!< Number of rows/columns */
    int const         isOneBased /*!< Whether COO is 0- or 1-based */
) {

    // ----- cannot assume that input is already 0!
    for (int l = 0; l < n+1; l++) col[l] = 0;

    // ----- find the correct column sizes
    for (int l = 0; l < nnz; l++)
        col[col_coo[l] - isOneBased]++;

    // ----- cumulative sum
    for (int i = 0, cumsum = 0; i < n; i++){
        int temp = col[i];
        col[i] = cumsum;
        cumsum += temp;
    }
    col[n] = nnz;

    // ----- copy the row indices to the correct place
    for (int l=0; l < nnz; l++){
        int col_l;
        col_l = col_coo[l] - isOneBased;

        int dst = col[col_l];
        row[dst] = row_coo[l] - isOneBased;

        col[col_l]++;
    }

    // ----- revert the column pointers
    for (int i=0, last=0; i<n; i++){
        int temp = col[i];
        col[i] = last;
        last = temp;
    }
}


int main(int argc, char *argv[]){

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nnz;
    int *coo_row, *coo_col;
    double *val;

    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0){
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) || mm_is_dense(matcode) || mm_is_hermitian(matcode)){
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nnz)) !=0)
        exit(1);


    /* reseve memory for matrices */

    coo_row = (int *) malloc(nnz * sizeof(int));
    coo_col = (int *) malloc(nnz * sizeof(int));
    val = (double *) malloc(nnz * sizeof(double));

    /* Replace missing val column with 1s and change the fscanf to match patter matrices*/

    if (!mm_is_pattern(matcode)){

        for (int i=0; i<nnz; i++)
            fscanf(f, "%d %d %lg\n", &coo_row[i], &coo_col[i], &val[i]);

    }
    else{

        for (int i=0; i<nnz; i++){
            fscanf(f, "%d %d\n", &coo_row[i], &coo_col[i]);
            val[i]=1;
        }
    }

    if (f !=stdin) fclose(f);

    /* Write out the matrix */
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nnz);
    //for (int i=0; i<nnz; i++)
    //    fprintf(stdout, "%d %d %20.19g\n", coo_row[i], coo_col[i], val[i]);

    int *csc_row = (int *)malloc(nnz*sizeof(int));
    int *csc_col = (int *)malloc((N+1)*sizeof(int));

    coo2csc(csc_row, csc_col,
            coo_row, coo_col,
            nnz, N, 1);

    printf("csc_col: ");
    for(int i=0; i<N+1; i++)
        printf("%d ", csc_col[i]);
    printf("\ncsc_row: ");
    for(int i=0; i<nnz; i++)
        printf("%d ", csc_row[i]);
    printf("\n");

    free(coo_row);
    free(coo_col);

    int triangles_found = 0;
    int c3[N];
    for(int i=0; i<N; i++)
        c3[i]=0;

    struct timespec start;
    struct timespec stop;
    struct timespec duration;

    clock_gettime(CLOCK_MONOTONIC, &start);

    pthread_mutex_t m;
    pthread_mutex_init(&m,NULL);

    //find the triangles
    cilk_for(int i=0; i<N-2; i++){
        cilk_for(int j=csc_col[i]; j<csc_col[i+1]; j++){
            cilk_for(int k=csc_col[csc_row[j]]; k<csc_col[csc_row[j]+1]; k++){
                if(j<csc_col[i+1]){
                    cilk_for(int l=j+1; l<csc_col[i+1]; l++){
                        if(csc_row[k] == csc_row[l]){
                            pthread_mutex_lock(&m);
                            triangles_found++;
                            c3[i]++;
                            c3[csc_row[j]]++;
                            c3[csc_row[l]]++;
                            printf("triangle found between nodes: (%d, %d, %d)\n", i, csc_row[j], csc_row[l]);
                            pthread_mutex_unlock(&m);
                        }
                    }
                }
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &stop);

    duration.tv_sec = stop.tv_sec - start.tv_sec;
    duration.tv_nsec = stop.tv_nsec - start.tv_nsec;

    if((duration.tv_nsec < 0) && (duration.tv_sec > 0)){
        duration.tv_sec--;
        duration.tv_nsec+=1000000000;
    }else if((duration.tv_nsec < 0) && (duration.tv_sec < 0)){
        printf("\n%s can't be running for negative time", argv[0]);
        exit(1);
    }

    //print results
    printf("id c3\n");
    for(int i=0; i<N; i++)
        printf("%d  %d\n", i, c3[i]);
    printf("\nTotal number of triangles: %d\n", triangles_found);
    printf("The counting of the triangles took %ld seconds and %ld nanoseconds.\n", duration.tv_sec, duration.tv_nsec);

    free(val);
    free(csc_col);
    free(csc_row);

    return 0;
}
