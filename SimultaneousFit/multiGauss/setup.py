from distutils.core import setup
import os, subprocess
from os import environ
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
setup(
    cmdclass={'multiGaus': build_ext},
    ext_modules=[Extension(
        "multi",                        # name of the extension
        sources=["multi.pyx"],  # additional source file(s)
        language="c++",          # generate C++ code
        #include_dirs=["/home/iaro/Programs/ROOT/root-6.14.06/builddir/include", ],
        #library_dirs=["/home/iaro/Programs/ROOT/root-6.14.06/builddir/lib", ],
        include_dirs=[environ['ROOTSYS']+"/include", ],
        library_dirs=[environ['ROOTSYS']+"/lib", ],
        libraries=["RooFit","RooFitCore"],
        # set complier options here:
        #extra_compile_args=["-I/home/iaro/Programs/ROOT/root-6.14.06/builddir/lib","-std=c++11"],
        extra_compile_args=["-I/"+environ['ROOTSYS']+"/lib","-std=c++11"],
    )])

