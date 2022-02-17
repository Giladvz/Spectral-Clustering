#ifndef spkmeans_H_
#define spkmeans_H_

/* Used in hw1 and 2. keeps the cluster's data */
typedef struct Clusters {
    int d;
    int size;
    double* pointsum;
} clus;

/* Used for index i,j values*/
typedef struct IntTuple {
    int i;
    int j;
} itup;

/* Used for index s,c values*/
typedef struct DoubleTuple {
    double s;
    double c;
} dtup;

/* Used to sort the eigen vectors*/
typedef struct eigen {
    double value;
    double* vector;
}eigen;

/* Used to return datapoint matrix and k to python */
typedef struct T {
    int k;
    double ** matrix;
}T;

void printMatrix(double** matrix, int numofx, int isJacobi);
void kmeans(double ** matrix,int k, int numofx, double * centroids,int fromC);
T fitC(double ** x, int d, int k, int numofx,char * goal, int fromC);


#endif
