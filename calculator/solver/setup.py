from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
 
ext_modules = Extension("lu_fast",["lu_fast.pyx"],\
    include_dirs=[numpy.get_include()],
    extra_compile_args=['-llapack -O3'],
    #extra_link_args=['-fopenmp']
    )
 
setup(
    name= 'lu',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext_modules]
)
