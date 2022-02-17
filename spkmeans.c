#define _CRT_SECURE_NO_WARNINGS
#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include "spkmeans.h"

#define size_of_type 8
#define max(a,b) (((a) > (b)) ? (a) : (b))

/* Returns distance between 2 points with d dimensions*/
double distance(double* x, double* y,int d)
{
    double v = 0.0;
    int k;
    for (k = 0; k<d;k++) {
        v += (double)((x[k]-y[k])*(x[k]-y[k]));
    }
    return v;
}

/* Creates WAM matrix, saved at WAMSpace*/
void getWAM(double ** x, int d,int numofx,double ** WAMSpace){
    int i =0;
    int j = 0;
    for (; i<numofx;i++){
        WAMSpace[i] = malloc(numofx*size_of_type);
        assert(WAMSpace[i]!=NULL && "An Error Has Occured");
        WAMSpace[i][i] = 0;
        for (j=0; j<i; j++){
            WAMSpace[i][j] = exp((-(sqrt(distance(x[i],x[j],d))))/2);
            WAMSpace[j][i] = WAMSpace[i][j];
           
        }
    }
}

/* Creates DDG values - returns only array of diag and not a full matrix
only in DDG will we need a matrix, so in order to save space and time
we chose to implement this way */
void getDDG(int numofx,double * DDGSpace,double ** WAMSpace){
    int i = 0;
    int j;
    double sum;
    sum = 0.0;
    for (; i<numofx;i++){
        for (j=0; j<numofx;j++){
            sum += WAMSpace[i][j];
        }
        DDGSpace[i] = sum;
        sum = 0.0;
    }
}

/* Takes DDG array and creates D^(-.5) */
double * getDMH(double * DDGSpace,int numofx) {
    int i =0;
    double * DMHSpace = malloc(numofx*size_of_type);
    assert(DMHSpace != NULL && "An Error Has Occured");
    for (; i <numofx;i++) {
        DMHSpace[i] = (1/sqrt(DDGSpace[i]));
    }
    return DMHSpace;
}

/*Takes DDG^(-.5) array and Wam matrix and creates Symmetric normalized
Laplacian matrix */
void getLnorm(int numofx,double * DDGSpace,
    double ** WAMSpace, double ** LSpace){
    int i=0;
    int j=0;
    double * DMHSpace = getDMH(DDGSpace,numofx);
    for (;i<numofx;i++){
        LSpace[i] = malloc(numofx*size_of_type);
        assert(LSpace[i] != NULL && "An Error Has Occured");
        for (j=0;j<numofx;j++){
            /*I - D^(-.5)*WAM*D^(-.5) */
            if (i==j){
                LSpace[i][j] = 1 - DMHSpace[i]*WAMSpace[i][j]*DMHSpace[j];
            }
            else {
                 LSpace[i][j] = - DMHSpace[i]*WAMSpace[i][j]*DMHSpace[j];
            }
            
        }
    }
    free(DMHSpace);
}

/* Get indexes of max value not on diagonal in matrix*/
itup getmax(double ** A,int numofx){
    int i = 0;
    int j;
    double max;
    itup ret;
    ret.i=0;
    ret.j=1;
    max = fabs(A[0][1]);
    for (; i <numofx;i++) {
        for (j=i+1;j<numofx;j++) {
            if (fabs(A[i][j]) > max) {
                max = fabs(A[i][j]);
                ret.i = i;
                ret.j = j;
                
            }
        }
    }
    return ret;
}

/* Gets values of s and c */
dtup getSC(double ** A, itup index){
    double degree; 
    double sign;
    dtup ret;
    double t;
    degree = (double)((A[index.j][index.j] - A[index.i][index.i]) /
         (2.0*A[index.i][index.j]));
    if (degree < 0) {
        sign = -1.0;
    }
    else {
        sign = 1.0;
    }
    t = (double)(sign /(fabs(degree) + sqrt((degree*degree) + 1.0)));
    ret.c = 1.0 / (sqrt((t*t) + 1.0));
    ret.s = t*ret.c;
    return ret;
}

/* Creates P matrix from i,j,s,c, saved at PSpace*/
void setP(double ** PSpace,int numofx, itup max, dtup SC){
    int i;
    for (i =0; i<numofx;i++){
        if (i == max.i || i == max.j){
            PSpace[i][i] = SC.c;
        }
        else {
            PSpace[i][i] = 1;
        }
    }
    PSpace[max.i][max.j] = SC.s;
    PSpace[max.j][max.i] = -SC.s;
}
/* Multiplies matrices */
double ** setV(double ** VSpace,double ** PSpace,int numofx){
    int i;
    int j;
    int r;
    double ** ret = calloc(numofx,size_of_type);
    assert(ret != NULL && "An Error Has Occured");
    for (i =0; i<numofx;i++){
        ret[i] = calloc(numofx,size_of_type);
        assert(ret[i]!=NULL && "An Error Has Occured");
        for (j=0; j<numofx;j++){
            for (r=0; r<numofx;r++){
                ret[i][j] += VSpace[i][r] * PSpace[r][j];
            }
        }
    }
    for (i=0;i<numofx;i++) {
        free(VSpace[i]);
    }
    free(VSpace);
    return ret;
}

/* Gets offset of symmetric matrix */
double getOf(double ** matrix, int numofx) {
    int i,j;
    double sum;
    sum= 0.0;
    for (i = 0; i<numofx;i++){
        for (j=i+1; j<numofx;j++){
            sum += (matrix[i][j]*matrix[i][j]);
        }
    }
    return sum*2;
}

/* Changes matrix on the spot (set newA = (P^t)AP) without the need of mult */
void setNewA(double ** A,int numofx,itup index,dtup SC){
    int r = 0;
    int i = index.i;
    int j = index.j;
    double tempri;
    double temprj;
    double tempii = A[i][i];
    double tempij = A[i][j];
    double tempjj = A[j][j];
    for (; r<numofx;r++){
        if (r != i && r != j){
            tempri = A[r][i];
            temprj = A[r][j];
            A[r][i] = SC.c*tempri - SC.s*temprj;
            /* Symmetic*/
            A[i][r] = A[r][i];
            A[r][j] = SC.c*temprj + SC.s*tempri;
            /* Symmetic*/
            A[j][r] = A[r][j];
        }
    }
    A[i][i] = (SC.c*SC.c*tempii) + (SC.s*SC.s*tempjj) - (2*SC.s*SC.c*tempij);
    A[j][j] = (SC.s*SC.s*tempii) + (SC.c*SC.c*tempjj) + (2*SC.s*SC.c*tempij);
    A[i][j] = 0;
    /* Symmetic*/
    A[j][i] = 0;
}

/* Takes A symmetric matrix and returns eigenvalues (A) 
   and eigenvectors (VSpace)*/
double ** getJacobi(double ** A,int numofx, double ** VSpace){
    int i;
    double offsetA;
    double offsetB;
    double ** PSpace;
    int count = 1;
    itup max;
    dtup SC;
    /* sets space for v matrix which is mostly 0 in the first run*/
    for (i = 0; i<numofx;i++){
        VSpace[i] = calloc(numofx,size_of_type);
        assert(VSpace[i] != NULL && "An Error Has Occured");
    }
    max = getmax(A,numofx);
    SC = getSC(A,max);
    /* First time, v=p */
    setP(VSpace ,numofx,max,SC);
    offsetA = getOf(A,numofx);
    setNewA(A,numofx,max,SC);
    offsetB = getOf(A,numofx);
    while ((offsetA-offsetB > pow(10,-15)) && count < 100) {
        count++;
        PSpace = calloc(numofx,size_of_type);
        assert(PSpace != NULL && "An Error Has Occured");
        for (i =0; i<numofx;i++){
            PSpace[i] = calloc(numofx,size_of_type);
            assert(PSpace[i] != NULL && "An Error Has Occured");
        }
        max = getmax(A,numofx);
        SC = getSC(A,max);
        setP(PSpace,numofx,max,SC);
        VSpace = setV(VSpace,PSpace,numofx);
        offsetA = getOf(A,numofx);
        setNewA(A,numofx,max,SC);
        offsetB = getOf(A,numofx);
        for (i = 0; i<numofx; i++){
            free(PSpace[i]);
        }
        free(PSpace);
    }
    return VSpace;
}

/* Create array of eigen structs, will be used to sort by eigenvalue*/
eigen* createEigen(double** A,double** VSpace,int numofx){
    int i,j;
    eigen* ret;
    ret = calloc(numofx,sizeof(eigen));
    assert(ret != NULL && "An Error Has Occured");
    for (i = 0; i<numofx;i++) {
        ret[i].vector = calloc(numofx,size_of_type);
        assert(ret[i].vector != NULL && "An Error Has Occured");
        for(j=0;j<numofx;j++){
            ret[i].vector[j] = VSpace[j][i];
        }
        ret[i].value = A[i][i]; 
    }
    return ret; 
}

/* Referenced by geeksforgeeks, part of mergeSort*/
void merge(eigen* eigens, int indexL, int indexR, int middle){
    int i = 0,j =0,k=indexL,len1,len2;
    eigen *list1, *list2;
    len1=(middle - indexL + 1);
    len2=(indexR - middle);
    list1 = (eigen*)calloc(len1,sizeof(eigen));
    assert(list1 != NULL && "An Error Has Occured");
    list2 = (eigen*)calloc(len2,sizeof(eigen));
    assert(list2 != NULL && "An Error Has Occured");

    
    while (i<max(len1,len2)) {
        if (i<len1) {
            list1[i] = eigens[indexL + i];
        }
        if (i<len2) {
            list2[i] = eigens[middle + i +1];
        }
        i++;
    }
    i = 0;
    while (i<len1 && j<len2) {
        if (list1[i].value <= list2[j].value) {
            eigens[k] = list1[i];
            i++;
        }
        else {
            eigens[k] = list2[j];
            j++;
        }
        k++;
    }
    while (i < len1) {
        eigens[k] = list1[i];
        i++;
        k++;
    }
    while (j < len2) {
        eigens[k] = list2[j];
        j++;
        k++;
    }
    free(list1);
    free(list2);
}

/*Decided to use MergeSort since it's a stable sort with 
big O notation of nlogn*/
void mergeSort(eigen* eigens, int indexL, int indexR) {
    int i;
    if (indexL<indexR) {
        i = (indexL+indexR)/2;
        mergeSort(eigens,indexL,i);
        mergeSort(eigens,i+1,indexR);
        merge(eigens,indexL,indexR,i);
    }
}
/* Prints matrices to screen,if isjacobi is set to 1,than transposes the
matrix */
void printMatrix(double** matrix, int numofx, int isJacobi){
    int i,j;
    if (isJacobi == 0) {
        for (i = 0; i<numofx; i++){
            for (j=0; j<numofx-1; j++) {
                if ((matrix[i][j] <= 0) && (matrix[i][j] > -0.00005)){
                    printf("0.0000,");
                }
                else {
                    printf("%.4f,", matrix[i][j]);
                }
            }
            if ((matrix[i][numofx-1] <= 0)&&(matrix[i][numofx-1] > -0.00005)){
                    printf("0.0000");
            }
            else {
                printf("%.4f", matrix[i][numofx-1]);
            }
            if (i!=numofx-1){
                printf("\n");
            }
        }
    }
    else {
        for (i = 0; i<numofx; i++){
            for (j=0; j<numofx-1; j++) {
                if ((matrix[j][i] <= 0) && (matrix[j][i] > -0.00005)){
                    printf("0.0000,");
                }
                else {
                    printf("%.4f,", matrix[j][i]);
                }
            }
            if ((matrix[numofx-1][i] <= 0)&&(matrix[numofx-1][i] > -0.00005)){
                    printf("0.0000");
            }
            else {
                printf("%.4f", matrix[numofx-1][i]);
            }
            if (i!=numofx-1){
                printf("\n");
            }
        }
    }
}

/* Used only when goal == ddg, creates a matrix from the array
   and prints to screen*/
void createDDGMatrix(double * array,int numofx){
    int i,j;
    double ** ret = calloc(numofx,size_of_type);
    assert(ret != NULL && "An Error Has Occured");
    for (i =0; i<numofx;i++){
        ret[i] = calloc(numofx,size_of_type);
        assert(ret[i]!= NULL && "An Error Has Occured");
        for (j=0; j<numofx;j++){
            if (i == j) {
                ret[i][j] = array[i];
            }
            else {
                ret[i][j] = 0;
            }
        }
    }
    printMatrix(ret,numofx,0);
    for (i = 0; i<numofx; i++){
        free(ret[i]);
    }
    free(ret);
}

/* Prints the first line of goal=jacobi, which includes the eigenvalues*/
void printEigenValues(double ** x,int numofx){
    int i;
    for (i = 0; i<numofx-1;i++){
        if ((x[i][i] <= 0) && (x[i][i] > -0.00005)){
            printf("0.0000");
        }
        else {
            printf("%.4f,", x[i][i]);
        }
    }
    if ((x[numofx-1][numofx-1] <= 0) && (x[numofx-1][numofx-1] > -0.00005)){
        printf("0.0000");
    }
    else {
        printf("%.4f", x[numofx-1][numofx-1]);
    }
    printf("\n");
}

/* If k == 0, gets K using the Eigengap Heuristic */
int getK(eigen * eigens,int numofx){
    int i;
    double max = 0.0;
    int ret;
    ret=0;
    for (i=0;i<floor((numofx) / 2) - 1;i++) {
        if (fabs(eigens[i].value - eigens[i+1].value) > max ) {
            max = fabs(eigens[i].value - eigens[i+1].value);
            ret = i;
        }
    }
    return ret+1;
}

/* If goal == spk, sorts the eigen structs by order of the values
Creates matrix U and T, and calls Kmeans or sends the info back to python
so kmeans++ algorithem would start*/
T spk(double ** VSpace, double ** LSpace, int numofx, int k,int fromC) {
    int i,j;
    double sum;
    T ret;
    eigen * eigens;
    double ** U;
    sum = 0.0;
    VSpace = getJacobi(LSpace,numofx,VSpace);
    eigens = createEigen(LSpace,VSpace,numofx);
    U = (double**)calloc(numofx,size_of_type);
    assert(U != NULL && "An Error Has Occured");

    mergeSort(eigens,0,numofx-1);

    if (k == 0){
        k = getK(eigens,numofx);
    }
    for (i=0;i<numofx;i++){
        U[i] = (double*)calloc(k,size_of_type);
        assert(U[i] != NULL && "An Error Has Occured");
        for (j=0;j<k;j++) {
            U[i][j] = eigens[j].vector[i]; 
            sum += eigens[j].vector[i]*eigens[j].vector[i];
        }
        /* Creates T from U*/
        for (j=0;j<k;j++) {
            U[i][j] = (U[i][j]) / (sqrt(sum));
        }
        sum=0.0;
    }
    if (fromC) {
        kmeans(U,k,numofx,NULL,1);
    }
    ret.k = k;
    ret.matrix = U;

    for (i=0; i<numofx;i++){
        free(VSpace[i]);
        free(eigens[i].vector);
        free(LSpace[i]);
    }
    free(VSpace);
    free(eigens);
    free(LSpace);
    
    return ret;
}

/* Acoording to goal, creates spaces needed for the matrices
and makes the matrices needed. Frees spaces not needed anymore and
calls spk if goal==spk. Implemented this way since except of jacobi,
every goal relies on the the one that comes before:
WAM-->DDG-->LNORM-->JACOBI-->SPK */
T fitC(double ** x, int d, int k, int numofx,char * goal, int fromC){
    int i;
    T ret;

    /* Jacobi has different input and doesn't rely on other matrices*/
    if (strcmp(goal,"jacobi") == 0) {
        double ** VSpace = calloc(numofx,size_of_type);
        assert(VSpace != NULL && "An Error Has Occured");
        VSpace = getJacobi(x,numofx,VSpace);
        printEigenValues(x,numofx);
        printMatrix(VSpace,numofx,1);
        for(i=0; i<numofx;i++) {
            free(VSpace[i]);
        }
        free(VSpace);
    }
    else {
        double ** WAMSpace = malloc(numofx*size_of_type);
        assert(WAMSpace != NULL && "An Error Has Occured");
        getWAM(x,d,numofx,WAMSpace);
        if (strcmp(goal,"wam") == 0){
            printMatrix(WAMSpace,numofx,0);
            for(i=0; i<numofx;i++) {
                free(WAMSpace[i]);
            }
            free(WAMSpace);
        }
        else {
            double * DDGSpace = malloc(numofx*size_of_type);
            assert(DDGSpace!= NULL && "An Error Has Occured");
            getDDG(numofx,DDGSpace,WAMSpace);       
            if ((strcmp(goal,"ddg") == 0)) {
                createDDGMatrix(DDGSpace,numofx);
                for(i=0; i<numofx;i++) {
                    free(WAMSpace[i]);
                }
                free(WAMSpace);
                free(DDGSpace);
            }
            else {
                double ** LSpace = malloc(numofx*size_of_type);
                assert(LSpace != NULL && "An Error Has Occured");
                getLnorm(numofx,DDGSpace,WAMSpace,LSpace);
                free(DDGSpace);
                if (strcmp(goal,"lnorm") == 0) {
                    printMatrix(LSpace,numofx,0);
                    for(i=0; i<numofx;i++) {
                        free(WAMSpace[i]);
                        free(LSpace[i]);
                    }
                    free(WAMSpace);
                    free(LSpace);
                }
                else {
                    /* goal is spk*/
                    double ** VSpace = calloc(numofx,size_of_type);
                    assert(VSpace != NULL && "An Error Has Occured");
                    ret = spk(VSpace,LSpace,numofx,k,fromC); 
                    for(i=0; i<numofx;i++) {
                        free(WAMSpace[i]);
                        free(x[i]);
                    }
                    free(WAMSpace);
                    free(x);

                    return ret;
                }
            }
        }
    }
    /* Has to return something*/
    ret.matrix = x;
    ret.k = k;
    return ret;
}


/* FUNCTIONS used by KMEANS*/
/* Adds the datapoints into the clusters */
void addallfromx(clus *cluslist, double **x, int i, int d,int n,int k){
    int clustoadd;
    double v[100];
    double u[100];
    int j;
    int t;
    double dis;
    double newdis;
    for(;i<n;i++){
        clustoadd = 0;
        for(j=0;j<d;j++){
            v[j] = x[i][j];    
            u[j] = cluslist[0].pointsum[j]; 
        }
        dis = distance(v,u,d);
        newdis = 0.0;
        for(j=0; j<k;j++)
        {
            for(t=0;t<d;t++){
                u[t] = cluslist[j].pointsum[t];
            }
            newdis = distance(u,v,d);
            if(newdis <= dis){
                clustoadd = j;
                dis = newdis;
            }
        }
        for(j=0;j<d;j++)
        {
             cluslist[clustoadd].pointsum[j+d] += v[j];
        }
        cluslist[clustoadd].size += 1;
    }
    for (i =0;i<k;i++){
        for (j=0;j<d;j++){
            if (cluslist[i].size != 0){
                cluslist[i].pointsum[j+(2*d)] = ((cluslist[i].pointsum[j+d]) 
                / (cluslist[i].size));
            }
        }
    }
}

/* Resets the clusters before a new cycle of adding datapoints*/
void setMid(clus* cluster,int d,int k)
{
    int i;
    int j;
    for (i =0;i<k;i++){
        for (j=0;j<d;j++){
            cluster[i].pointsum[j] = cluster[i].pointsum[j+(2*d)];
            cluster[i].size=0;
            cluster[i].pointsum[j+d]=0.0;
        }
    }
}

/* Checks if the clusters havn't changed after cycle, if so, the kmeans
algorithem finished */
int fincheck(clus* cluslist,int d,int k) {
    int i;
    int j;
    for (j = 0; j < k; j++) {
        for (i =0;i < d; i++) {
            if (cluslist[j].pointsum[i] != cluslist[j].pointsum[i+2*d]) {
                return 0;
            }
        }
    }
    return 1;
}

/* Creates the clusters from the first k points and adds remaining datapoints.
After that runs in cycles untill algorithem finishes and prints to screen*/
void kmeans(double ** x,int k,int numofx,double * centroids, int fromC)
{
    int r,l,s;
    clus* clusters;
    int maxiter = 300;
    clusters = (clus*)calloc(k,sizeof(clus));
    assert(clusters != NULL && "An Error Has Occured");
    
    for (r=0;r < k;r++){
        /* pointsum:0->d point, d->2d sum, 2d->3d mid*/
        (clusters[r]).pointsum = (double*)calloc(3 * (k) , size_of_type);
        assert(clusters[r].pointsum != NULL && "An Error Has Occured");
        (clusters[r]).d = k;
        if (fromC == 1) {
            (clusters[r]).size = 1;
            for (l = 0; l < k;l++) {
                (clusters[r]).pointsum[(k+l)] = x[r][l];
                (clusters[r]).pointsum[l] = x[r][l];
            }
            
        }
        else {
            (clusters[r]).size = 0;
            for (l = 0; l < k;l++) {
                (clusters[r]).pointsum[l] = centroids[(r*k+l)];
            }
        }
    }

    if (fromC == 1) {
        addallfromx(clusters, x,k, k,numofx,k);
    }
    else {
        addallfromx(clusters, x,0, k,numofx,k);
    }
    maxiter--;

    while (maxiter > 0 && (fincheck(clusters,k,k)==0)){
        setMid(clusters,k,k);
        addallfromx(clusters,x,0, k,numofx,k);
        maxiter--;
    }
    for (s = 0; s<k; s++){
        for (r=0; r<k-1; r++) {
            printf("%.4f,", clusters[s].pointsum[r]);
        }
        printf("%.4f", clusters[s].pointsum[k-1]);
        if (s != (k-1)){
                printf("\n");
            }
    }

    for (s=0; s<k;s++){
        free(clusters[s].pointsum);
    }
    free(clusters);
}

/* Main function, gets inputs from user and creates a matrix of the datapoints,
*/
int main(int argc,char * argv[]){
    /* 9 chars for number, max 10 numbers, max 9 commas, enter char*/
    char line[100];
    int d;
    int buff_size;
    int buff_used;
    int i;
    int j;
    T ret;
    
    char *token;
    double **x;
    double * observation;
    FILE * input;

    d=argc;
    d=0;
    i=0;
    j=0;
    buff_used=0;
    buff_size=25;

    input = fopen(argv[3],"r");
    assert(input!=NULL && "An Error Has Occured");
    fgets(line, 100 , input);
    token = strtok(line, ",");
    while (token != NULL)
    {
        token = strtok(NULL, ",\n");
        d++;
    }
    /* x is an array of pointers, each pointer points to an observation*/
    x = (double**)calloc(buff_size,size_of_type);
    assert(x!=NULL && "An Error Has Occured");
    rewind(input);

    while(fgets(line, 1000 , input) != NULL) {
        
        token = strtok(line, ",");
        if (buff_size == buff_used) {
            buff_size*=2;
            x = realloc(x,buff_size*size_of_type);
            assert (x != NULL && "An Error Has Occured");
        }
        observation = (double*)calloc(d,size_of_type);
        assert(observation!=NULL && "An Error Has Occured");
        while (token != NULL)
        {
            observation[j] = atof(token);
            token = strtok(NULL, ",\n");
            j++;
        }
        j = 0;
        x[i] = observation;
        i++;
        buff_used++;
    }
    ret = fitC(x,d,atoi(argv[1]),buff_used,argv[2],1);
    for (i=0;i<buff_used;i++){
        free(ret.matrix[i]);
    }
    free(ret.matrix);
    fclose(input);
    return 0;
}


















