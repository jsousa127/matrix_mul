#include <stdio.h>  
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#define MAX 50

float **fst_matrix;
float **snd_matrix;
float **result_matrix;
double clearcache[30000000];

float** createMatrix(int value, int size) {
    float** ret = malloc(sizeof(float*) * size);
    int i, j;
    for(i = 0; i < size; i++) {
        ret[i] = malloc(sizeof(float) * size);
        for(j = 0; j < size; j++) {
            if(value == -1) // -1 for random numbers
                ret[i][j] = (float) (rand()) / ((float) (RAND_MAX/MAX));
            else 
                ret[i][j] = value;    
        }      
    }
    return ret;
}

void freeMatrix(float** matrix, int size) {
    int i, j;
    for(i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void printMatrix(float** matrix, int size) {
    int i, j;
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++)
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

void trans(float** matrix, int size) {
    int i, j;
    float** aux = malloc(sizeof(float*) * size);
    for(i = 0; i < size; i++){
        aux[i] = malloc(sizeof(float) * size);
        for(j = 0; j < size; j++) 
            aux[i][j] = matrix[j][i];
    }
    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) 
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
                c[i][j] += a[i][k] * b[k][j]; 
    return;
}

void mul_t(float **a, float **b, float **c, int size) {
    float sum;
    int i, j, k;
    trans(b,size);
    for (i = 0; i < size; i++) 
        for (j = 0; j < size; j++) 
            for (k = 0; k < size; k++)
                c[i][j] += a[i][k] * b[j][k];   
    return;
}

void mulIKJ(float **a, float **b, float **c, int size) {
    int i, j, k;
    for (i = 0; i < size; i++) {
        for (k = 0; k < size; k++) 
            for (j = 0; j < size; j++) 
                c[i][j] += a[i][k] * b[k][j];
    }
    return;
}

void mulJKI(float **a, float **b, float **c, int size) {
    float sum;
    int i, j, k;

    for (j = 0; j < size; j++) {
        for (k = 0; k < size; k++) 
            for (i = 0; i < size; i++) 
                c[i][j] += a[i][k] * b[k][j];
    }
    return;
}

void mulJKI_t(float **a, float **b, float **c, int size) {
    float sum;
    int i, j, k;

    trans(a, size);
    trans(b, size);

    for (j = 0; j < size; j++) {
        for (k = 0; k < size; k++) 
            for (i = 0; i < size; i++) 
                c[j][i] += a[k][i] * b[j][k];
    }
    trans(c,size);
    return;
}

double dtime() {
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime,(struct timezone*)0);
    tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    return( tseconds );
}

static int cmpdouble (const void * a, const void * b)
{
  if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;  
}

int kBest(double* array,int size, int k, double tol) {
    int i;
    double sigma, top;
    if(size < 3) return 0;
    qsort(array, size, sizeof(double), cmpdouble);
    sigma = (double)(array[0] * 0.05);
    top = (double)(array[0] + sigma);
    for(i=1; i<k; i++) {
        if(array[i] > top) return 0;
    }
    return 1;
}

void runFunc(void (*f)(float**, float**, float **, int), int size, double* t) {
    int i,run, sizet=0;
    double t0, t1, ret;
    for(run = 0; run < 8; run++){    
        clearCache();
        t0 = dtime();
        (*f)(fst_matrix,snd_matrix,result_matrix,size);
        t1 = (dtime() - t0);
        t[sizet] = (double) (t1 * 1000.0);
        sizet++;
        free(result_matrix);
        result_matrix = createMatrix(0,size);
        if(kBest(t,sizet,3,5.0)) {
            printf("%f;",t[0]);
            return;
        }
    }
    printf("%f*;",t[0]);
    return;
}   





int main(int argc, char const *argv[]) {
    int i, j;
    double *t = malloc(sizeof(double)*8);
    double ijk, ikj, jki;
    int size = atoi(argv[1]);
    srand(getpid());
    
    fst_matrix = createMatrix(-1, size); 
    snd_matrix = createMatrix(1, size);
    result_matrix = createMatrix(0,size);
    
    runFunc(mul,size,t);
    runFunc(mul_t,size,t);
    runFunc(mulIKJ,size,t);
    runFunc(mulJKI,size,t);
    runFunc(mulJKI_t,size,t);

    freeMatrix(fst_matrix,size);
    freeMatrix(snd_matrix,size);
    freeMatrix(result_matrix,size);
    free(t);
    return 0;
}
