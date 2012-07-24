#include <Python.h>

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include "flamingo-4.1/src/filtertree/src/wrappers/wrappers.h"
#include "flamingo-4.1/src/stringmap/src/editdistance.h"

using namespace std;


typedef struct {
    PyObject_HEAD
    GramGenFixedLen * gramGen;
    StringContainerVector * strContainer;
    WrapperSimpleEd * index;
} flamingo_WrapperSimpleEd;


static int WrapperSimpleEd_init(flamingo_WrapperSimpleEd *self, PyObject *args, PyObject *kwds)
{
    int numLines;
    char * line;
    PyObject * listObj; /* the list of strings */
    PyObject * strObj;  /* one string in the list */
    
    /* the O! parses for a Python object (listObj) checked
       to be of type PyList_Type */
    if (! PyArg_ParseTuple( args, "O!", &PyList_Type, &listObj))
        return NULL;

    numLines = PyList_Size(listObj);

    if (numLines < 0)
        return NULL; /* Not a list */

    vector<string> barcodes;

    /* iterate over items of the list, grabbing strings, and parsing
       for numbers */
    int i;
    for (i=0; i<numLines; i++){
	    strObj = PyList_GetItem(listObj, i);
	    line = PyString_AsString( strObj );
	    barcodes.push_back(line);
    }

    self->gramGen = new GramGenFixedLen(2);
    self->strContainer = new StringContainerVector(true);
    self->strContainer->initStatsCollector(self->gramGen);
    self->strContainer->fillContainer(barcodes.begin(), barcodes.end());
    self->index = new WrapperSimpleEd(self->strContainer,
                                            self->gramGen, true);
    self->index->buildIndex();
    return 0;
}

static PyObject *
WrapperSimpleEd_search(flamingo_WrapperSimpleEd* self, PyObject *args) {
    const char * query;
    float editDistance;
    if (!PyArg_ParseTuple(args, "sf", &query, &editDistance))
        return NULL;

    vector<unsigned> resultStringIDs;
    string tmp;
    self->index->search(query, editDistance, resultStringIDs);
    PyObject * ret = PyList_New(resultStringIDs.size());
    for(unsigned i = 0; i < resultStringIDs.size(); i++) {
        self->strContainer->retrieveString(tmp, resultStringIDs[i]);
        PyList_SetItem(ret, i, Py_BuildValue("s", tmp.c_str()));
    }
    return ret;
}


static PyMethodDef WrapperSimpleEdMethods[] = {
    {"search", (PyCFunction)WrapperSimpleEd_search, METH_VARARGS,
     "Return the name, combining the first and last name"
    },
    {NULL}  /* Sentinel */
};


unsigned int levenshtein_distance(const string &s1, const string & s2) {
    const size_t len1 = s1.size(), len2 = s2.size();
    vector<unsigned int> col(len2+1), prevCol(len2+1);

    for (unsigned int i = 0; i < prevCol.size(); i++)
            prevCol[i] = i;
    for (unsigned int i = 0; i < len1; i++) {
            col[0] = i+1;
            for (unsigned int j = 0; j < len2; j++)
                    col[j+1] = min( min( 1 + col[j], 1 + prevCol[1 + j]),
                                                            prevCol[j] + (s1[i]==s2[j] ? 0 : 1) );
            col.swap(prevCol);
    }
    return prevCol[len2];
}


static PyObject *
flamingo_distance(PyObject *self, PyObject *args) {
    const char * c1;
    const char * c2;
    if (!PyArg_ParseTuple(args, "ss", &c1, &c2))
        return NULL;

    string s1(c1);
    string s2(c2);

    return PyInt_FromSsize_t(levenshtein_distance(s1, s2));
}


static PyMethodDef FlamingoMethods[] = {
        {"distance", flamingo_distance, METH_VARARGS,
         "Caculate the Levenshtein edit distance."},
        {NULL, NULL, 0, NULL}               /* Sentinel */
};

static PyTypeObject flamingo_WrapperSimpleEdType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "flamingo.WrapperSimpleEd",             /*tp_name*/
    sizeof(flamingo_WrapperSimpleEd),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "flamingo simple wrapper",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    WrapperSimpleEdMethods,             /* tp_methods */ 
    0,             /* tp_members TODO */
    0,           /* tp_getset TODO */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)WrapperSimpleEd_init,      /* tp_init */
    0,                         /* tp_alloc */
    0,                 /* tp_new TODO */
};


PyMODINIT_FUNC
initflamingo(void)
{
    PyObject* m;
    m = Py_InitModule("flamingo", FlamingoMethods);
    flamingo_WrapperSimpleEdType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&flamingo_WrapperSimpleEdType) < 0)
        return;
    PyModule_AddObject(m, "WrapperSimpleEd", (PyObject *) &flamingo_WrapperSimpleEdType);
}

int
main(int argc, char *argv[])
{
        /* Pass argv[0] to the Python interpreter */
        Py_SetProgramName(argv[0]);

        /* Initialize the Python interpreter.    Required. */
        Py_Initialize();

        /* Add a static module */
        initflamingo();
        return 0;
}
