// have all the functs used to call from python
#include <stdio.h>
#include <stdlib.h>
#include <Python.h>
#include "spkmeans.h"
#include <string.h>
#include <assert.h>

#define size_of_type 8

double ** toMatrix(PyObject * pylist, Py_ssize_t size, Py_ssize_t len) {
    int i;
    int j;
    PyObject * item;
    double ** clist = malloc(size*size_of_type);
    assert(clist != NULL && "An Error Has Occured");
    for (i=0;i<size;i++) {
        clist[i] = malloc(size_of_type*len);
        assert(clist[i] != NULL && "An Error Has Occured");

        for (j=0;j<len;j++) {
            item = PyList_GetItem(PyList_GetItem(pylist,i),j);
            clist[i][j] = PyFloat_AsDouble(item);
        }
    }
    return clist;
}

double * toList(PyObject * pylist, Py_ssize_t size, Py_ssize_t len) {
    int i;
    int j;
    PyObject * item;
    double * clist = malloc(size*size_of_type*len);
    assert(clist != NULL);
    for (i=0;i<size;i++) {
        for (j=0;j<len;j++) {
            item = PyList_GetItem(PyList_GetItem(pylist,i),j);
            clist[i*len+j] = PyFloat_AsDouble(item);
        }
    }
    return clist;
}

static PyObject * fit(PyObject *self, PyObject *args){
    PyObject *x_p, * pylist,*pyd,*pynumpy;
    Py_ssize_t d_p,k_p,numofx_p;
    char * goal;
    double ** x_c;
    T ret_c;
    int d_c,k_c,numofx_c,s,r;
    if (!(PyArg_ParseTuple(args,"Osii",&x_p,&goal,&k_p,&d_p))) {
        printf("An Error Has Occured");
        return NULL;
    }
    /* set python integers, will need to convert to c integer*/
    numofx_p = PyList_Size(x_p);
    /* Sets c integers*/
    numofx_c = (int)numofx_p;
    k_c = (int)k_p;
    d_c = (int)d_p;
    /* Sets c lists from python arg*/
    x_c = toMatrix(x_p,numofx_c,d_c);
    ret_c = fitC(x_c,d_c,k_c,numofx_c,goal,0);
    pylist = PyList_New(numofx_c);
    if (!pylist) {
        printf("An Error Has Occured");
        return NULL;
    }
    d_c = ret_c.k;
    /* Create 2d array to return to python*/
    for (s = 0; s<numofx_c; s++){
        pyd = PyList_New(d_c);
        if (!pyd){
            printf("An Error Has Occured");
            return NULL;
        }
        for (r=0; r<d_c; r++) {
            pynumpy = PyFloat_FromDouble(ret_c.matrix[s][r]);
            if (!pynumpy) {
                printf("An Error Has Occured");
                return NULL;
            }
            PyList_SET_ITEM(pyd,r,pynumpy);
        } 
        PyList_SET_ITEM(pylist,s,pyd);
    }
    for (s=0;s<numofx_c;s++){
        free(ret_c.matrix[s]);
    }
    free(ret_c.matrix);
    return pylist;
}

static PyObject * Kfit(PyObject *self, PyObject *args){
    PyObject *x_p,*centroids_p;
    Py_ssize_t d_p,k_p,numofx_p;
    double **x_c;
    double *centroids_c;
    int d_c,k_c,numofx_c;

    if (!(PyArg_ParseTuple(args,"OO",&x_p,&centroids_p))) {
        printf("An Error Has Occured");
        return NULL;
    }
    else {    
        /* set python integers, will need to convert to c integer*/
        numofx_p = PyList_Size(x_p);
        k_p = PyList_Size(centroids_p);
        d_p = PyList_Size(PyList_GetItem(x_p,0));
        /* Sets c integers*/
        numofx_c = (int)numofx_p;
        k_c = (int)k_p;
        d_c = (int)d_p;
        /* Sets c lists from python arg*/
        x_c = toMatrix(x_p,numofx_c,d_c);
        centroids_c = toList(centroids_p,k_c,d_c);

        kmeans(x_c,k_c,numofx_c,centroids_c,0);
        free(centroids_c);
        for (d_c = 0; d_c < k_c; d_c++){
            free(x_c[d_c]);
        }
        free(x_c);
        return x_p;
    }
    return x_p;
}

static PyMethodDef mykmeansMethods[] = {
    {"fit",
    (PyCFunction) fit,
    METH_VARARGS
    },
    {"Kfit",
    (PyCFunction) Kfit,
    METH_VARARGS
    },
    {NULL,NULL,0,NULL}
};

static struct PyModuleDef moduleDef = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule",
    NULL,
    -1,
    mykmeansMethods
};
PyMODINIT_FUNC
PyInit_spkmeansmodule(void) {
    PyObject * m;
    m = PyModule_Create(&moduleDef);
    if (!m) {
        return NULL;
    }
    return m;
}











