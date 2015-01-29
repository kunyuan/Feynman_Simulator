from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
 
ext_modules = Extension("lu",["solver.pyx"],\
    include_dirs=[numpy.get_include()],
    extra_compile_args=['-llapack'],
    #extra_link_args=['-fopenmp']
    )
 
setup(
    name= 'lu',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext_modules]
)
