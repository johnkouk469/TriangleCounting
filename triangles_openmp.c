/*
** v4.c -- version 4 of Vertexwise triangle counting
**
** Loading sparse adjacency matrix using the Matrix Market format
** and converting the martix from COO to CSC
**
** Representing graph in CSC format
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "mmio.h"

void quicksort(int element_list[], int low, int high){
	int pivot, value1, value2, temp;
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

    int *cooFull_row = (int *) malloc((2*nnz)*sizeof(int));
    int *cooFull_col = (int *) malloc((2*nnz)*sizeof(int));
    int *valFull = (int *) malloc((2*nnz)*sizeof(int));

    for(int i=0; i<nnz; i++){
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

    int *csc_row = (int *)malloc((2*nnz)*sizeof(int));
    int *csc_col = (int *)malloc((N+1)*sizeof(int));

    coo2csc(csc_row, csc_col,
            cooFull_row, cooFull_col,
            2*nnz, N, 1);

    free(cooFull_row);
    free(cooFull_col);

    for(int i=0; i<N+1; i++){
        quicksort(csc_row, csc_col[i], csc_col[i+1]-1);
    }
    
    // printf("csc_col: ");
    // for(int i=0; i<N+1; i++)
    //     printf("%d ", csc_col[i]);
    // printf("\ncsc_row: ");
    // for(int i=0; i<2*nnz; i++)
    //     printf("%d ", csc_row[i]);
    // printf("\n");
    // printf("nnz: %d\n", nnz);
    

    int c3[N];

    for(int i=0; i<N; i++)
        c3[i]=0;


    struct timespec start;
    struct timespec stop;
    struct timespec duration;

    clock_gettime(CLOCK_MONOTONIC, &start);
    
    #pragma omp parallel
    #pragma omp for
	for(int j=0; j<N; j++){
        
        int nzrangeOfColA = csc_col[j+1]-csc_col[j];
        int colA[nzrangeOfColA];
        // printf("\nj=%d",j);
        // printf("\nnzrangeOfColA: %d \ncolA: ", nzrangeOfColA);
        for(int y=csc_col[j]; y<csc_col[j+1]; y++){
            colA[y-csc_col[j]] = csc_row[y];
            // printf("%d ",csc_row[y]);
        }
        // printf("\n");

        for(int n=csc_col[j]; n<csc_col[j+1]; n++){
            
            int i = csc_row[n];
            /*
            ** Iterate all the non zero values of matrix A
            ** A(i,j) !=  0 
            */
            // printf("\n(i,j)=(%d,%d)\n", i, j);

            int nnzrangeOfRowA = csc_col[i+1]-csc_col[i];       
            int rowA[nnzrangeOfRowA];
            // printf("\nnnzrangeOfRowA: %d \nrowA: ", nnzrangeOfRowA);        
            for(int x=csc_col[i]; x<csc_col[i+1]; x++){
                rowA[x-csc_col[i]] = csc_row[x];
                // printf("%d ",csc_row[x]);
            }
            // printf("\n");
            
        
            // for(int z=0; z<nnzrangeOfRowA; z++){
            //     for(int w=0; w<nzrangeOfColA; w++){
            //         if(rowA[z]<colA[w]){
            //             break;
            //         }else if(rowA[z]=colA[w]){
            //             c3[j]++;
            //             printf("common");
            //             break;
            //         }
            //     }
            // }
        
            int common = 0;
            int flag = 0;
                for(int l=0; l<nnzrangeOfRowA; l++){
                    int counter = 0;
                    while((counter + flag) < nzrangeOfColA){
                    if(rowA[l] < colA[counter+flag]){
                        counter++;
                        break;
                    }else if(rowA[l] == colA[counter+flag]){
                        common++;
                        break;
                    }else
                        flag++;
                }
            }
            c3[j] += common;
            
        }        
    }

    clock_gettime(CLOCK_MONOTONIC, &stop);

    free(csc_row);
    free(csc_col);

    // printf("\nC3:\n");
    // for(int i=0; i<N; i++){
    //     if(c3[i]%2 != 0)
    //         c3[i]++;
    //     c3[i] = c3[i]/2;
    //     printf("%d %d\n", i, c3[i]);
    // }

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
