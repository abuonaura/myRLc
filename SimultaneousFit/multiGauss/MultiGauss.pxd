cdef extern from "MultiGauss.cpp":
    pass
    
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string



cdef extern from "MultiGauss.h" namespace "multifit":
    cdef cppclass MultiGauss:
        MultiGauss() except +
        map[string, vector[double] ] defineValues();
        map[string, vector[double] ]  function();
