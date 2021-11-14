/*
** Loading sparse adjacency matrix using the Matrix Market format
** and converting the martix from COO to CSC
**
** Representing graph in CSC format
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "mmio.h"

void quicksort(uint32_t element_list[], uint32_t low, uint32_t high){
	uint32_t pivot, value1, value2, temp;
	if (low < high){
		pivot = low;
		value1 = low;
		value2 = high;
		while (value1 < value2){
			while (element_list[value1] <= element_list[pivot] && value1 <= high){
				value1++;
			}
			while (element_list[value2] > element_list[pivot] && value2 >= low){
				value2--;
			}
			if (value1 < value2){
				temp = element_list[value1];
				element_list[value1] = element_list[value2];
				element_list[value2] = temp;
			}
		}
		temp = element_list[value2];
		element_list[value2] = element_list[pivot];
		element_list[pivot] = temp;
		quicksort(element_list, low, value2 - 1);
		quicksort(element_list, value2 + 1, high);
	}
}

void coo2csc(
    uint32_t       * const row,       /*!< CSC row start indices */
    uint32_t       * const col,       /*!< CSC column indices */
    uint32_t const * const row_coo,   /*!< COO row indices */
    uint32_t const * const col_coo,   /*!< COO column indices */
    uint32_t const         nnz,       /*!< Number of nonzero elements */
    uint32_t const         n,         /*!< Number of rows/columns */
    uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
) {

    // ----- cannot assume that input is already 0!
    for (uint32_t l = 0; l < n+1; l++) col[l] = 0;

    // ----- find the correct column sizes
    for (uint32_t l = 0; l < nnz; l++)
        col[col_coo[l] - isOneBased]++;

    // ----- cumulative sum
    for (uint32_t i = 0, cumsum = 0; i < n; i++){
        uint32_t temp = col[i];
        col[i] = cumsum;
        cumsum += temp;
    }
    col[n] = nnz;

    // ----- copy the row indices to the correct place
    for (uint32_t l=0; l < nnz; l++){
        uint32_t col_l;
        col_l = col_coo[l] - isOneBased;

        uint32_t dst = col[col_l];
        row[dst] = row_coo[l] - isOneBased;

        col[col_l]++;
    }

    // ----- revert the column pointers
    for (uint32_t i=0, last=0; i<n; i++){
        uint32_t temp = col[i];
        col[i] = last;
        last = temp;
    }
}


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

    uint32_t *cooFull_row = (uint32_t *) malloc((2*nnz)*sizeof(uint32_t));
    uint32_t *cooFull_col = (uint32_t *) malloc((2*nnz)*sizeof(uint32_t));
    uint32_t *valFull = (uint32_t *) malloc((2*nnz)*sizeof(uint32_t));

    for(uint32_t i=0; i<nnz; i++){
        cooFull_row[i] = coo_row[i];
        cooFull_row[nnz+i] = coo_col[i];
        cooFull_col[i] = coo_col[i];
        cooFull_col[nnz+i] = coo_row[i];
        valFull[i] = 0;
        valFull[nnz+i] = 0;
    }

    /* Write out the matrix */
    /*
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nnz);
    for (int i=0; i<nnz; i++)
        fprintf(stdout, "%d %d %g   %d %d\n",
                coo_row[i], coo_col[i], val[i], cooFull_row[nnz+i], cooFull_col[nnz+i]);
    */

    free(coo_row);
    free(coo_col);
    free(val);

    uint32_t *csc_row = (uint32_t *)malloc((2*nnz)*sizeof(uint32_t));
    uint32_t *csc_col = (uint32_t *)malloc((N+1)*sizeof(uint32_t));

    coo2csc(csc_row, csc_col,
            cooFull_row, cooFull_col,
            2*nnz, N, 1);

    free(cooFull_row);
    free(cooFull_col);

    for(uint32_t i=0; i<N+1; i++){
        quicksort(csc_row, csc_col[i], csc_col[i+1]-1);
    }

    /*
    printf("csc_col: ");
    for(int i=0; i<N+1; i++)
        printf("%d ", csc_col[i]);
    printf("\ncsc_row: ");
    for(int i=0; i<2*nnz; i++)
        printf("%d ", csc_row[i]);
    printf("\n");
    */

    uint32_t c3[N];

    for(uint32_t i=0; i<N; i++)
        c3[i]=0;


    struct timespec start;
    struct timespec stop;
    struct timespec duration;

    clock_gettime(CLOCK_MONOTONIC, &start);

	for(uint32_t i=0; i<N+1; i++){ // i the rows of matrix 1
        uint32_t sum = 0;
        if(csc_col[i+1] > csc_col[i]){
            uint32_t ioc1[csc_col[i+1]-csc_col[i]];
            for(uint32_t k=csc_col[i]; k<csc_col[i+1]; k++)
                ioc1[k-csc_col[i]] = csc_row[k]; // list of the column indices for row i of the non zero values of matrix 1

            for(uint32_t j=0; j<N+1; j++){ // j the columns of matrix 2
                if(csc_col[j+1] > csc_col[j]){

                    uint32_t ior2[csc_col[j+1]-csc_col[j]];
                    for(uint32_t w=csc_col[j]; w<csc_col[j+1]; w++)
                        ior2[w-csc_col[j]] = csc_row[w]; // list of the row indices for column j of the non zero values of matrix 2

                    uint32_t common = 0;
                    uint32_t flag = 0;
                    for(uint32_t l=0; l<csc_col[i+1]-csc_col[i]; l++){
                        uint32_t counter = 0;
                        while((counter + flag) < (csc_col[j+1]-csc_col[j])){
                            if(ioc1[l] < ior2[counter+flag]){
                                counter++;
                                break;
                            }else if(ioc1[l] == ior2[counter+flag]){
                                common++;
                                break;
                            }else
                                flag++;
                        }
                    }
                    // Found A^2(i,j)=common
                    if(common != 0){
                        for(uint32_t k=csc_col[i]; k<csc_col[i+1]; k++){
                            if(csc_row[k] == j){
                                sum += common;
                                break;
                            }
                        }
                    }
                }
            }
        }
        c3[i] = sum;

    }

    clock_gettime(CLOCK_MONOTONIC, &stop);

    free(csc_row);
    free(csc_col);

    printf("\nC3:\n");
    for(uint32_t i=0; i<N; i++){
        if(c3[i]%2 != 0)
            c3[i]++;
        c3[i] = c3[i]/2;
        printf("%d %d\n", i, c3[i]);
    }

    duration.tv_sec = stop.tv_sec - start.tv_sec;
    duration.tv_nsec = stop.tv_nsec - start.tv_nsec;

    if((duration.tv_nsec < 0) && (duration.tv_sec > 0)){
        duration.tv_sec--;
        duration.tv_nsec+=1000000000;
    }else if((duration.tv_nsec < 0) && (duration.tv_sec < 0)){
        printf("\n%s can't be running for negative time", argv[0]);
        exit(1);
    }

    printf("The prossess took %ld seconds and %ld nanoseconds", duration.tv_sec, duration.tv_nsec);


	return 0;
}
