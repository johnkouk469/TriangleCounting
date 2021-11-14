/*
** v2.c -- version 2 of Vertexwise triangle counting
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_NODE 8000

int main(int argc, char* argv[]){

    int n;
    srand(time(NULL));

    //input management
    if(argc != 2){
        printf("Usage: %s n\n where n is no. of nodes\n", argv[0]);
        exit(1);
    }

    n = atoi(argv[1]);

    if((n < 1) || (n > MAX_NODE)){
        printf("The no of nodes should be between 1 and %d.\n", MAX_NODE);
        exit(1);
    }

    int count[n];
    for(int i=0; i<n; i++){
        count[i]=0;
    }

    //create the graph
    int A[n][n];
    for(int i=0; i<n; i++){
        A[i][i]=0;
        for(int j=0; j<i; j++){
            A[i][j] = rand() % 2;
            A[j][i] = A[i][j];
        }
    }

    //print graph
    /*
    for(int i=0; i<n; i++){
        printf("\n");
        for(int j=0; j<n; j++){
            printf("%d ",A[i][j]);
        }
    }
    printf("\n");
    */

    struct timespec start;
    struct timespec stop;
    struct timespec duration;

    clock_gettime(CLOCK_MONOTONIC, &start);

    //find the triangles
    for(int i=0; i<n-2; i++){
        for(int j=i+1; j<n-1; j++){
            for(int k=j+1; k<n; k++){
                if(A[i][j] == 1 && A[j][k] == 1 && A[k][i] == 1){
                    count[i]++;
                    count[j]++;
                    count[k]++;
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
    }else if((duration.tv_nsec < 0) && (duration.tv_sec <0)){
        printf("%s can't be running for negative time", argv[0]);
        exit(1);
    }

    //print results
    printf("%s was running for %ld seconds and %ld nanoseconds.\n", argv[0], duration.tv_sec, duration.tv_nsec);
    /*
    for(int i=0; i<n; i++){
        printf("%d ",count[i]);
    }
    printf("\n");
    */
    exit(0);
}
