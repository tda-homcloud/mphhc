#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "mphhc/core.hpp"
#include <vector>
#include <new>
#include <iostream>

// --- mphhc.get_version() ---

static PyObject* mphhc_get_version(PyObject* self, PyObject* args) {
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }
    std::string version = mphhc::get_version();
    return PyUnicode_FromString(version.c_str());
}

// --- mphhc.boundary_matrix class ---

typedef struct {
    PyObject_HEAD
    mphhc::boundary_matrix* bm;
} MatrixObject;

static void Matrix_dealloc(MatrixObject* self) {
    if (self->bm) {
        delete self->bm;
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* Matrix_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    MatrixObject* self;
    self = (MatrixObject*)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->bm = NULL;
    }
    return (PyObject*)self;
}

static int Matrix_init(MatrixObject* self, PyObject* args, PyObject* kwds) {
    int maxdim;
    if (!PyArg_ParseTuple(args, "i", &maxdim)) {
        return -1;
    }
    try {
        self->bm = new mphhc::boundary_matrix(maxdim);
    } catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create boundary_matrix instance");
        return -1;
    }
    return 0;
}

static bool ensure_bm_initialized(MatrixObject* self) {
    if (self->bm == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "boundary_matrix not initialized");
        return false;
    }
    return true;
}

static PyObject* Matrix_max_dim(MatrixObject* self, PyObject* args) {
    if (!ensure_bm_initialized(self)) return NULL;
    return PyLong_FromLong(self->bm->max_dim());
}

static PyObject* Matrix_num_simplices(MatrixObject* self, PyObject* args) {
    if (!ensure_bm_initialized(self)) return NULL;
    return PyLong_FromLong(self->bm->num_simplices());
}

static PyObject* Matrix_is_reduced(MatrixObject* self, PyObject* args) {
    if (!ensure_bm_initialized(self)) return NULL;
    if (self->bm->is_reduced()) {
        Py_RETURN_TRUE;
    } else {
        Py_RETURN_FALSE;
    }
}

static PyObject* Matrix_set_dim_col(MatrixObject* self, PyObject* args) {
    if (!ensure_bm_initialized(self)) return NULL;

    mphhc::index idx;
    int dim;
    PyObject* col_obj;
    if (!PyArg_ParseTuple(args, "iiO", &idx, &dim, &col_obj)) {
        return NULL;
    }

    if (idx != self->bm->num_simplices()) {
        PyErr_SetString(PyExc_TypeError, "index must be an incremental integer");
        return NULL;
    }

    if (!PySequence_Check(col_obj)) {
        PyErr_SetString(PyExc_TypeError, "col must be a sequence");
        return NULL;
    }

    mphhc::column col;
    Py_ssize_t len = PySequence_Size(col_obj);
    col.reserve(len);

    for (Py_ssize_t i = 0; i < len; i++) {
        PyObject* item = PySequence_GetItem(col_obj, i);
        if (!item) {
            return NULL;
        }
        if (!PyLong_Check(item)) {
            Py_DECREF(item);
            PyErr_SetString(PyExc_TypeError, "column elements must be integers");
            return NULL;
        }
        col.push_back((mphhc::index)PyLong_AsLong(item));
        Py_DECREF(item);
    }

    self->bm->set_dim_col(idx, dim, col);
    return PyLong_FromLong(idx);
}

static PyObject* Matrix_reduce_standard(MatrixObject* self, PyObject* args) {
    if (!ensure_bm_initialized(self)) return NULL;
    self->bm->reduce_standard();
    Py_RETURN_NONE;
}

static PyObject* Matrix_reduce_twist(MatrixObject* self, PyObject* args) {
    if (!ensure_bm_initialized(self)) return NULL;
    self->bm->reduce_twist();
    Py_RETURN_NONE;
}

static PyObject* Matrix_birth_death_pairs(MatrixObject* self, PyObject* args) {
    if (!ensure_bm_initialized(self)) return NULL;
    
    if (!self->bm->is_reduced()) {
        PyErr_SetString(PyExc_RuntimeError, "Matrix must be reduced before calling birth_death_pairs. Call reduce_standard() or reduce_twist() first.");
        return NULL;
    }
    
    std::vector<mphhc::birth_death_pair> pairs = self->bm->birth_death_pairs();
    
    PyObject* list = PyList_New(pairs.size());

    if (!list) return NULL;

    for (size_t i = 0; i < pairs.size(); ++i) {
        PyObject* tuple;
        if (pairs[i].death < 0)
          tuple = Py_BuildValue("(iiO)", pairs[i].dim, pairs[i].birth, Py_None);
        else
          tuple = Py_BuildValue("(iii)", pairs[i].dim, pairs[i].birth, pairs[i].death);

        if (!tuple) {
            Py_DECREF(list);
            return NULL;
        }
        PyList_SetItem(list, i, tuple); // Steals reference to tuple
    }

    return list;
}

static PyMethodDef Matrix_methods[] = {
    {"max_dim", (PyCFunction)Matrix_max_dim, METH_NOARGS, "Return max dimension"},
    {"num_simplices", (PyCFunction)Matrix_num_simplices, METH_NOARGS, "Return number of simplices"},
    {"is_reduced", (PyCFunction)Matrix_is_reduced, METH_NOARGS, "Return whether the matrix is reduced"},
    {"set_dim_col", (PyCFunction)Matrix_set_dim_col, METH_VARARGS, "Add a column for a dimension"},
    {"reduce_standard", (PyCFunction)Matrix_reduce_standard, METH_NOARGS, "Perform standard reduction"},
    {"reduce_twist", (PyCFunction)Matrix_reduce_twist, METH_NOARGS, "Perform twist reduction"},
    {"birth_death_pairs", (PyCFunction)Matrix_birth_death_pairs, METH_NOARGS, "Get birth-death pairs"},
    {NULL}  /* Sentinel */
};

static PyTypeObject MatrixType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "mphhc.Matrix",   /* tp_name */
    sizeof(MatrixObject), /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Matrix_dealloc, /* tp_dealloc */
    0,                         /* tp_vectorcall_offset */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_as_async */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
    "Boundary Matrix objects", /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Matrix_methods,    /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Matrix_init, /* tp_init */
    0,                         /* tp_alloc */
    Matrix_new,        /* tp_new */
};

static PyMethodDef mphhc_methods[] = {
    {"get_version", mphhc_get_version, METH_VARARGS, "Get version string"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mphhc_module = {
    PyModuleDef_HEAD_INIT,
    "mphhc",   /* name of module */
    NULL,      /* module documentation, may be NULL */
    -1,        /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    mphhc_methods
};

PyMODINIT_FUNC PyInit_mphhc(void) {
    PyObject* m;

    if (PyType_Ready(&MatrixType) < 0)
        return NULL;

    m = PyModule_Create(&mphhc_module);
    if (m == NULL)
        return NULL;

    Py_INCREF(&MatrixType);
    if (PyModule_AddObject(m, "Matrix", (PyObject*)&MatrixType) < 0) {
        Py_DECREF(&MatrixType);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
