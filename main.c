#include <stdio.h>  
#include <stdlib.h>
#include <sys/time.h>

#define N 3
#define MAX 50

double clearcache[30000000];

float** createMatrix(int value) {
    float** ret = malloc(sizeof(float*) * N);
    int i, j;
    for(i = 0; i < N; i++) {
        ret[i] = malloc(sizeof(float) * N);
        for(j = 0; j < N; j++) {
            if(value == -1) // -1 for random numbers
                ret[i][j] = (float) (rand()) / ((float) (RAND_MAX/MAX));
            else 
                ret[i][j] = value;    
        }      
    }
    return ret;
}

void freeMatrix(float** matrix) {
    int i, j;
    for(i = 0; i < N; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void printMatrix(float** matrix) {
    int i, j;
    for(i = 0; i < N; i++) {
        for(j = 0; j < N; j++)
            printf("%f ",matrix[i][j]); 
        printf("\n");
    }
    printf("\n");
}

void clearCache (void) {
	int i;
	for (i = 0; i < 30000000; i++)
		clearcache[i] = i;
}

void trans(float** matrix) {
    int i, j;
    float** aux = malloc(sizeof(float*) * N);
    for(i = 0; i < N; i++){
        aux[i] = malloc(sizeof(float) * N);
        for(j = 0; j < N; j++) 
            aux[i][j] = matrix[j][i];
    }
    for(i = 0; i < N; i++) {
        for(j = 0; j < N; j++) 
            matrix[i][j] = aux[i][j];
        free(aux[i]);    
    }   
    free(aux);        
}   

void mul(float **a, float **b, float **c, int size) {
    float sum;
    int i, j, k;
    for (i = 0; i < size; i++) 
        for (j = 0; j < size; j++) 
            for (k = 0; k < size; k++)
                //  c[i][j] += a[i][k] * b[k][j]; (trans b = b)
                c[i][j] += a[i][k] * b[k][j];   
    printMatrix(c);
    return;
}

void mulIKJ(float **a, float **b, float **c, int size) {
    int i, j, k;
    for (i = 0; i < size; i++) {
        for (k = 0; k < size; k++) 
            for (j = 0; j < size; j++) 
                c[i][j] += a[i][k] * b[k][j];
    }
    printMatrix(c);
    return;
}

void mulJKI(float **a, float **b, float **c, int size) {
    float sum;
    int i, j, k;

    trans(a);
    // trans(b) = b

    for (j = 0; j < size; j++) {
        for (k = 0; k < size; k++) 
            for (i = 0; i < size; i++) 
                //c[i][j] += a[i][k] * b[k][j];
                c[i][j] += a[k][i] * b[j][k];
    }
    printMatrix(c);
    return;
}

int main(int argc, char const *argv[]) {
    int i, j;
    
    float **fst_matrix = createMatrix(-1); ;
    float **snd_matrix = createMatrix(1);
    float **result_matrix = createMatrix(0);
    
    clearCache();
    
    mul(fst_matrix,snd_matrix,result_matrix,N);
    
    result_matrix = createMatrix(0);
    mulIKJ(fst_matrix,snd_matrix,result_matrix,N);
    
    result_matrix = createMatrix(0);
    mulJKI(fst_matrix,snd_matrix,result_matrix,N);
    
    freeMatrix(fst_matrix);
    freeMatrix(snd_matrix);
    freeMatrix(result_matrix);
    
    return 0;
}
