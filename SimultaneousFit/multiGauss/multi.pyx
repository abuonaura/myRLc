#distutils: language = c++

from MultiGauss cimport MultiGauss
from ROOT import *

cdef class PyMultiGauss:
    cdef MultiGauss* c_multig
    
    def __cinit__(self):
        self.c_multig = new MultiGauss()
    def __dealloc__(self):
        del self.c_multig
        
    def get_parameters(self):
        return self.c_multig.function()