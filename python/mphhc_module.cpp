#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <iostream>
#include <new>
#include <string>
#include <vector>

#include "mphhc/core.hpp"

// --- mphhc.GetVersion() ---

static PyObject *mphhc_get_version(PyObject *self, PyObject *args) {
  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }
  std::string version = mphhc::GetVersion();
  return PyUnicode_FromString(version.c_str());
}

// --- mphhc.BoundaryMatrix class ---

typedef struct {
  PyObject_HEAD mphhc::BoundaryMatrix *bm;
} MatrixObject;

static void Matrix_dealloc(MatrixObject *self) {
  if (self->bm) {
    delete self->bm;
  }
  Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *Matrix_new(PyTypeObject *type, PyObject *args,
                            PyObject *kwds) {
  MatrixObject *self;
  self = (MatrixObject *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->bm = NULL;
  }
  return (PyObject *)self;
}

static int Matrix_init(MatrixObject *self, PyObject *args, PyObject *kwds) {
  int maxdim;
  int save_basis = 0;  // default false
  static char *kwlist[] = {(char *)"maxdim", (char *)"save_basis", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|p", kwlist, &maxdim,
                                   &save_basis)) {
    return -1;
  }
  try {
    self->bm = new mphhc::BoundaryMatrix(maxdim, save_basis);
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Failed to create boundary_matrix instance");
    return -1;
  }
  return 0;
}

static bool ensure_bm_initialized(MatrixObject *self) {
  if (self->bm == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "boundary_matrix not initialized");
    return false;
  }
  return true;
}

static PyObject *Matrix_max_dim(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;
  return PyLong_FromLong(self->bm->MaxDim());
}

static PyObject *Matrix_num_simplices(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;
  return PyLong_FromLong(self->bm->NumSimplices());
}

static PyObject *Matrix_is_reduced(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;
  if (self->bm->IsReduced()) {
    Py_RETURN_TRUE;
  } else {
    Py_RETURN_FALSE;
  }
}

static PyObject *Matrix_is_save_basis(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;
  if (self->bm->IsSaveBasis()) {
    Py_RETURN_TRUE;
  } else {
    Py_RETURN_FALSE;
  }
}

static PyObject *Matrix_set_dim_col(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;

  int idx;
  int dim;
  PyObject *col_obj;
  if (!PyArg_ParseTuple(args, "iiO", &idx, &dim, &col_obj)) {
    return NULL;
  }

  if (idx != self->bm->NumSimplices()) {
    PyErr_SetString(PyExc_TypeError, "index must be an incremental integer");
    return NULL;
  }

  if (!PySequence_Check(col_obj)) {
    PyErr_SetString(PyExc_TypeError, "col must be a sequence");
    return NULL;
  }

  mphhc::Column col;
  Py_ssize_t len = PySequence_Size(col_obj);
  col.reserve(len);

  for (Py_ssize_t i = 0; i < len; i++) {
    PyObject *item = PySequence_GetItem(col_obj, i);
    if (!item) {
      return NULL;
    }
    if (!PyLong_Check(item)) {
      Py_DECREF(item);
      PyErr_SetString(PyExc_TypeError, "column elements must be integers");
      return NULL;
    }
    col.push_back((mphhc::Index)PyLong_AsLong(item));
    Py_DECREF(item);
  }

  self->bm->SetDimCol(idx, dim, col);
  return PyLong_FromLong(idx);
}

static PyObject *Matrix_reduce_standard(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;
  self->bm->ReduceStandard();
  Py_RETURN_NONE;
}

static PyObject *Matrix_reduce_twist(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;
  self->bm->ReduceTwist();
  Py_RETURN_NONE;
}

static PyObject *Matrix_reduce(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;

  char *algorithm;

  if (!PyArg_ParseTuple(args, "z", &algorithm)) return NULL;

  if (algorithm == NULL || std::string(algorithm) == "mphhc-twist") {
    self->bm->ReduceTwist();
  } else if (std::string(algorithm) == "mphhc-standard") {
    self->bm->ReduceStandard();
  } else {
    return PyErr_Format(PyExc_NotImplementedError, "Unknown algoithm: %s",
                        algorithm);
  }

  Py_RETURN_NONE;
}

static PyObject *Matrix_birth_death_pairs(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;

  if (!self->bm->IsReduced()) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Matrix must be reduced before calling birth_death_pairs. "
                    "Call reduce_standard() or reduce_twist() first.");
    return NULL;
  }

  std::vector<mphhc::BirthDeathPair> pairs = self->bm->BirthDeathPairs();

  PyObject *list = PyList_New(pairs.size());

  if (!list) return NULL;

  for (size_t i = 0; i < pairs.size(); ++i) {
    PyObject *tuple;
    if (pairs[i].death < 0)
      tuple = Py_BuildValue("(iiO)", pairs[i].dim, pairs[i].birth, Py_None);
    else
      tuple =
          Py_BuildValue("(iii)", pairs[i].dim, pairs[i].birth, pairs[i].death);

    if (!tuple) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SetItem(list, i, tuple);  // Steals reference to tuple
  }

  return list;
}

static PyObject *Matrix_basis(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;

  if (!self->bm->IsReduced()) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Matrix must be reduced before calling basis.");
    return NULL;
  }

  if (!self->bm->IsSaveBasis()) {
    PyErr_SetString(PyExc_RuntimeError,
                    "save_basis must be True to retrieve basis.");
    return NULL;
  }

  std::vector<mphhc::Column> basis = self->bm->Basis();

  PyObject *list = PyList_New(basis.size());
  if (!list) return NULL;

  for (size_t i = 0; i < basis.size(); ++i) {
    PyObject *col_list = PyList_New(basis[i].size());
    if (!col_list) {
      Py_DECREF(list);
      return NULL;
    }
    for (size_t j = 0; j < basis[i].size(); ++j) {
      PyObject *val = PyLong_FromLong(basis[i][j]);
      if (!val) {
        Py_DECREF(col_list);
        Py_DECREF(list);
        return NULL;
      }
      PyList_SetItem(col_list, j, val);
    }
    PyList_SetItem(list, i, col_list);
  }
  return list;
}

static PyObject *Matrix_lows(MatrixObject *self, PyObject *args) {
  if (!ensure_bm_initialized(self)) return NULL;

  std::vector<mphhc::Index> lows = self->bm->Lows();

  PyObject *list = PyList_New(lows.size());
  if (!list) return NULL;

  for (size_t i = 0; i < lows.size(); ++i) {
    PyObject *val = PyLong_FromLong(lows[i]);
    if (!val) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SetItem(list, i, val);
  }
  return list;
}

static PyMethodDef Matrix_methods[] = {
    {"max_dim", (PyCFunction)Matrix_max_dim, METH_NOARGS,
     "Return max dimension"},
    {"num_simplices", (PyCFunction)Matrix_num_simplices, METH_NOARGS,
     "Return number of simplices"},
    {"is_reduced", (PyCFunction)Matrix_is_reduced, METH_NOARGS,
     "Return whether the matrix is reduced"},
    {"is_save_basis", (PyCFunction)Matrix_is_save_basis, METH_NOARGS,
     "Return whether the matrix saves basis"},
    {"set_dim_col", (PyCFunction)Matrix_set_dim_col, METH_VARARGS,
     "Add a column for a dimension"},
    {"reduce_standard", (PyCFunction)Matrix_reduce_standard, METH_NOARGS,
     "Perform standard reduction"},
    {"reduce_twist", (PyCFunction)Matrix_reduce_twist, METH_NOARGS,
     "Perform twist reduction"},
    {"reduce", (PyCFunction)Matrix_reduce, METH_VARARGS, "Perform reduction"},
    {"birth_death_pairs", (PyCFunction)Matrix_birth_death_pairs, METH_NOARGS,
     "Get birth-death pairs"},
    {"basis", (PyCFunction)Matrix_basis, METH_NOARGS, "Get basis"},
    {"lows", (PyCFunction)Matrix_lows, METH_NOARGS, "Get lows"},
    {NULL} /* Sentinel */
};

static PyTypeObject MatrixType = {
    PyVarObject_HEAD_INIT(NULL, 0) "mphhc.Matrix", /* tp_name */
    sizeof(MatrixObject),                          /* tp_basicsize */
    0,                                             /* tp_itemsize */
    (destructor)Matrix_dealloc,                    /* tp_dealloc */
    0,                                             /* tp_vectorcall_offset */
    0,                                             /* tp_getattr */
    0,                                             /* tp_setattr */
    0,                                             /* tp_as_async */
    0,                                             /* tp_repr */
    0,                                             /* tp_as_number */
    0,                                             /* tp_as_sequence */
    0,                                             /* tp_as_mapping */
    0,                                             /* tp_hash  */
    0,                                             /* tp_call */
    0,                                             /* tp_str */
    0,                                             /* tp_getattro */
    0,                                             /* tp_setattro */
    0,                                             /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,      /* tp_flags */
    "Boundary Matrix objects",                     /* tp_doc */
    0,                                             /* tp_traverse */
    0,                                             /* tp_clear */
    0,                                             /* tp_richcompare */
    0,                                             /* tp_weaklistoffset */
    0,                                             /* tp_iter */
    0,                                             /* tp_iternext */
    Matrix_methods,                                /* tp_methods */
    0,                                             /* tp_members */
    0,                                             /* tp_getset */
    0,                                             /* tp_base */
    0,                                             /* tp_dict */
    0,                                             /* tp_descr_get */
    0,                                             /* tp_descr_set */
    0,                                             /* tp_dictoffset */
    (initproc)Matrix_init,                         /* tp_init */
    0,                                             /* tp_alloc */
    Matrix_new,                                    /* tp_new */
};

// --- mphhc.FlatBoundaryMatrix class ---

typedef struct {
  PyObject_HEAD mphhc::FlatBoundaryMatrix *bm;
} FlatMatrixObject;

static void FlatMatrix_dealloc(FlatMatrixObject *self) {
  if (self->bm) {
    delete self->bm;
  }
  Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *FlatMatrix_new(PyTypeObject *type, PyObject *args,
                                PyObject *kwds) {
  FlatMatrixObject *self;
  self = (FlatMatrixObject *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->bm = NULL;
  }
  return (PyObject *)self;
}

static int FlatMatrix_init(FlatMatrixObject *self, PyObject *args,
                           PyObject *kwds) {
  int maxdim;
  int save_basis = 0;  // default false
  static char *kwlist[] = {(char *)"maxdim", (char *)"save_basis", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|p", kwlist, &maxdim,
                                   &save_basis)) {
    return -1;
  }
  try {
    self->bm = new mphhc::FlatBoundaryMatrix(maxdim, save_basis);
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Failed to create flat_boundary_matrix instance");
    return -1;
  }
  return 0;
}

static bool ensure_flat_bm_initialized(FlatMatrixObject *self) {
  if (self->bm == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "flat_boundary_matrix not initialized");
    return false;
  }
  return true;
}

static PyObject *FlatMatrix_num_simplices(FlatMatrixObject *self,
                                          PyObject *args) {
  if (!ensure_flat_bm_initialized(self)) return NULL;
  return PyLong_FromLong(self->bm->NumSimplices());
}

static PyObject *FlatMatrix_is_reduced(FlatMatrixObject *self, PyObject *args) {
  if (!ensure_flat_bm_initialized(self)) return NULL;
  if (self->bm->IsReduced()) {
    Py_RETURN_TRUE;
  } else {
    Py_RETURN_FALSE;
  }
}

static PyObject *FlatMatrix_is_save_basis(FlatMatrixObject *self,
                                          PyObject *args) {
  if (!ensure_flat_bm_initialized(self)) return NULL;
  if (self->bm->IsSaveBasis()) {
    Py_RETURN_TRUE;
  } else {
    Py_RETURN_FALSE;
  }
}

static PyObject *FlatMatrix_set_dim_col(FlatMatrixObject *self,
                                        PyObject *args) {
  if (!ensure_flat_bm_initialized(self)) return NULL;

  int idx;
  int dim;
  PyObject *col_obj;
  if (!PyArg_ParseTuple(args, "iiO", &idx, &dim, &col_obj)) {
    return NULL;
  }

  if (idx != self->bm->NumSimplices()) {
    PyErr_SetString(PyExc_TypeError, "index must be an incremental integer");
    return NULL;
  }

  if (!PySequence_Check(col_obj)) {
    PyErr_SetString(PyExc_TypeError, "col must be a sequence");
    return NULL;
  }

  mphhc::Column col;
  Py_ssize_t len = PySequence_Size(col_obj);
  col.reserve(len);

  for (Py_ssize_t i = 0; i < len; i++) {
    PyObject *item = PySequence_GetItem(col_obj, i);
    if (!item) {
      return NULL;
    }
    if (!PyLong_Check(item)) {
      Py_DECREF(item);
      PyErr_SetString(PyExc_TypeError, "column elements must be integers");
      return NULL;
    }
    col.push_back((mphhc::Index)PyLong_AsLong(item));
    Py_DECREF(item);
  }

  self->bm->SetDimCol(idx, dim, col);
  return PyLong_FromLong(idx);
}

static PyObject *FlatMatrix_reduce(FlatMatrixObject *self, PyObject *args) {
  if (!ensure_flat_bm_initialized(self)) return NULL;
  self->bm->Reduce();
  Py_RETURN_NONE;
}

static PyObject *FlatMatrix_birth_death_pairs(FlatMatrixObject *self,
                                              PyObject *args) {
  if (!ensure_flat_bm_initialized(self)) return NULL;

  if (!self->bm->IsReduced()) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Matrix must be reduced before calling birth_death_pairs. "
                    "Call reduce() first.");
    return NULL;
  }

  std::vector<mphhc::BirthDeathPair> pairs = self->bm->BirthDeathPairs();

  PyObject *list = PyList_New(pairs.size());

  if (!list) return NULL;

  for (size_t i = 0; i < pairs.size(); ++i) {
    PyObject *tuple;
    if (pairs[i].death < 0)
      tuple = Py_BuildValue("(iiO)", pairs[i].dim, pairs[i].birth, Py_None);
    else
      tuple =
          Py_BuildValue("(iii)", pairs[i].dim, pairs[i].birth, pairs[i].death);

    if (!tuple) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SetItem(list, i, tuple);  // Steals reference to tuple
  }

  return list;
}

static PyObject *FlatMatrix_basis(FlatMatrixObject *self, PyObject *args) {
  if (!ensure_flat_bm_initialized(self)) return NULL;

  if (!self->bm->IsReduced()) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Matrix must be reduced before calling basis.");
    return NULL;
  }

  if (!self->bm->IsSaveBasis()) {
    PyErr_SetString(PyExc_RuntimeError,
                    "save_basis must be True to retrieve basis.");
    return NULL;
  }

  std::vector<mphhc::Column> basis = self->bm->Basis();

  PyObject *list = PyList_New(basis.size());
  if (!list) return NULL;

  for (size_t i = 0; i < basis.size(); ++i) {
    PyObject *col_list = PyList_New(basis[i].size());
    if (!col_list) {
      Py_DECREF(list);
      return NULL;
    }
    for (size_t j = 0; j < basis[i].size(); ++j) {
      PyObject *val = PyLong_FromLong(basis[i][j]);
      if (!val) {
        Py_DECREF(col_list);
        Py_DECREF(list);
        return NULL;
      }
      PyList_SetItem(col_list, j, val);
    }
    PyList_SetItem(list, i, col_list);
  }
  return list;
}

static PyObject *FlatMatrix_lows(FlatMatrixObject *self, PyObject *args) {
  if (!ensure_flat_bm_initialized(self)) return NULL;

  std::vector<mphhc::Index> lows = self->bm->Lows();

  PyObject *list = PyList_New(lows.size());
  if (!list) return NULL;

  for (size_t i = 0; i < lows.size(); ++i) {
    PyObject *val = PyLong_FromLong(lows[i]);
    if (!val) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SetItem(list, i, val);
  }
  return list;
}

static PyMethodDef FlatMatrix_methods[] = {
    {"num_simplices", (PyCFunction)FlatMatrix_num_simplices, METH_NOARGS,
     "Return number of simplices"},
    {"is_reduced", (PyCFunction)FlatMatrix_is_reduced, METH_NOARGS,
     "Return whether the matrix is reduced"},
    {"is_save_basis", (PyCFunction)FlatMatrix_is_save_basis, METH_NOARGS,
     "Return whether the matrix saves basis"},
    {"set_dim_col", (PyCFunction)FlatMatrix_set_dim_col, METH_VARARGS,
     "Add a column for a dimension"},
    {"reduce", (PyCFunction)FlatMatrix_reduce, METH_NOARGS,
     "Perform reduction"},
    {"birth_death_pairs", (PyCFunction)FlatMatrix_birth_death_pairs,
     METH_NOARGS, "Get birth-death pairs"},
    {"basis", (PyCFunction)FlatMatrix_basis, METH_NOARGS, "Get basis"},
    {"lows", (PyCFunction)FlatMatrix_lows, METH_NOARGS, "Get lows"},
    {NULL} /* Sentinel */
};

static PyTypeObject FlatMatrixType = {
    PyVarObject_HEAD_INIT(NULL, 0) "mphhc.FlatMatrix", /* tp_name */
    sizeof(FlatMatrixObject),                          /* tp_basicsize */
    0,                                                 /* tp_itemsize */
    (destructor)FlatMatrix_dealloc,                    /* tp_dealloc */
    0,                                        /* tp_vectorcall_offset */
    0,                                        /* tp_getattr */
    0,                                        /* tp_setattr */
    0,                                        /* tp_as_async */
    0,                                        /* tp_repr */
    0,                                        /* tp_as_number */
    0,                                        /* tp_as_sequence */
    0,                                        /* tp_as_mapping */
    0,                                        /* tp_hash  */
    0,                                        /* tp_call */
    0,                                        /* tp_str */
    0,                                        /* tp_getattro */
    0,                                        /* tp_setattro */
    0,                                        /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
    "Flat Boundary Matrix objects",           /* tp_doc */
    0,                                        /* tp_traverse */
    0,                                        /* tp_clear */
    0,                                        /* tp_richcompare */
    0,                                        /* tp_weaklistoffset */
    0,                                        /* tp_iter */
    0,                                        /* tp_iternext */
    FlatMatrix_methods,                       /* tp_methods */
    0,                                        /* tp_members */
    0,                                        /* tp_getset */
    0,                                        /* tp_base */
    0,                                        /* tp_dict */
    0,                                        /* tp_descr_get */
    0,                                        /* tp_descr_set */
    0,                                        /* tp_dictoffset */
    (initproc)FlatMatrix_init,                /* tp_init */
    0,                                        /* tp_alloc */
    FlatMatrix_new,                           /* tp_new */
};

static PyMethodDef mphhc_methods[] = {
    {"get_version", mphhc_get_version, METH_VARARGS, "Get version string"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef mphhc_module = {
    PyModuleDef_HEAD_INIT, "mphhc", /* name of module */
    NULL,                           /* module documentation, may be NULL */
    -1, /* size of per-interpreter state of the module, or -1 if the module
           keeps state in global variables. */
    mphhc_methods};

PyMODINIT_FUNC PyInit_mphhc(void) {
  PyObject *m;

  if (PyType_Ready(&MatrixType) < 0) return NULL;

  m = PyModule_Create(&mphhc_module);
  if (m == NULL) return NULL;

  Py_INCREF(&MatrixType);
  if (PyModule_AddObject(m, "Matrix", (PyObject *)&MatrixType) < 0) {
    Py_DECREF(&MatrixType);
    Py_DECREF(m);
    return NULL;
  }

  if (PyType_Ready(&FlatMatrixType) < 0) return NULL;
  Py_INCREF(&FlatMatrixType);
  if (PyModule_AddObject(m, "FlatMatrix", (PyObject *)&FlatMatrixType) < 0) {
    Py_DECREF(&FlatMatrixType);
    Py_DECREF(m);
    return NULL;
  }

  return m;
}
