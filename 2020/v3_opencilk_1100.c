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
#include <stdint.h>
#include <time.h>
#include <cilk/cilk.h>
#include <pthread.h>
#include "mmio.h"


int main(int argc, char *argv[]){

    uint32_t ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t M, N, nnz;
    uint32_t *coo_row, *coo_col;
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

    coo_row = (uint32_t *) malloc(nnz * sizeof(uint32_t));
    coo_col = (uint32_t *) malloc(nnz * sizeof(uint32_t));
    val = (double *) malloc(nnz * sizeof(double));

    /* Replace missing val column with 1s and change the fscanf to match patter matrices*/

    if (!mm_is_pattern(matcode)){

        for (uint32_t i=0; i<nnz; i++)
            fscanf(f, "%d %d %lg\n", &coo_row[i], &coo_col[i], &val[i]);

    }
    else{

        for (uint32_t i=0; i<nnz; i++){
            fscanf(f, "%d %d\n", &coo_row[i], &coo_col[i]);
            val[i]=1;
        }
    }

    if (f !=stdin) fclose(f);

    /* Write out the matrix */
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nnz);
    //for (uint32_t i=0; i<nnz; i++)
    //    fprintf(stdout, "%d %d %20.19g\n", coo_row[i], coo_col[i], val[i]);

    uint32_t *csc_row = (uint32_t *)malloc(nnz*sizeof(uint32_t));
    uint32_t *csc_col = (uint32_t *)malloc((N+1)*sizeof(uint32_t));

        // ----- cannot assume that input is already 0!
    for (uint32_t l = 0; l < N+1; l++) csc_col[l] = 0;

    // ----- find the correct column sizes
    for (uint32_t l = 0; l < nnz; l++)
        csc_col[coo_col[l] - 1]++;

    // ----- cumulative sum
    for (uint32_t i = 0, cumsum = 0; i < N; i++){
        uint32_t temp = csc_col[i];
        csc_col[i] = cumsum;
        cumsum += temp;
    }
    csc_col[N] = nnz;

    // ----- copy the row indices to the correct place
    for (uint32_t l=0; l < nnz; l++){
        uint32_t col_l;
        col_l = coo_col[l] - 1;

        uint32_t dst = csc_col[col_l];
        csc_row[dst] = coo_row[l] - 1;

        csc_col[col_l]++;
    }

    // ----- revert the column pointers
    for (uint32_t i=0, last=0; i<N; i++){
        uint32_t temp = csc_col[i];
        csc_col[i] = last;
        last = temp;
    }

    free(coo_row);
    free(coo_col);

    uint32_t triangles_found = 0;
    uint32_t c3[N];
    for(uint32_t i=0; i<N; i++)
        c3[i]=0;

    struct timespec start;
    struct timespec stop;
    struct timespec duration;

    printf("Vertexwise triangle counting for %s ...", argv[1]);

    clock_gettime(CLOCK_MONOTONIC, &start);

    pthread_mutex_t m;
    //pthread_mutex_init(&m,NULL);

    //find the triangles
    pthread_mutex_init(&m,NULL);
    cilk_for(uint32_t i=0; i<N-2; i++){
        cilk_for(uint32_t j=csc_col[i]; j<csc_col[i+1]; j++){
            //pthread_mutex_lock(&m);
            for(uint32_t k=csc_col[csc_row[j]]; k<csc_col[csc_row[j]+1]; k++){
                if(j<csc_col[i+1]){
                    for(uint32_t l=j+1; l<csc_col[i+1]; l++){
                        pthread_mutex_lock(&m);
                        if(csc_row[k] == csc_row[l]){
                            //pthread_mutex_lock(&m);
                            triangles_found++;
                            c3[i]++;
                            c3[csc_row[j]]++;
                            c3[csc_row[l]]++;
                            //printf("triangle found between nodes: (%d, %d, %d)\n", i, csc_row[j], csc_row[l]);
                            //pthread_mutex_unlock(&m);
                        }
                        pthread_mutex_unlock(&m);
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
    //printf("id c3\n");
    //for(uint32_t i=0; i<N; i++)
    //    printf("%d  %d\n", i, c3[i]);
    printf("\nTotal number of triangles: %d\n", triangles_found);
    printf("The counting of the triangles took %ld seconds and %ld nanoseconds.\n", duration.tv_sec, duration.tv_nsec);

    free(val);
    free(csc_col);
    free(csc_row);

	return 0;
}
