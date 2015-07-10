from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy,os

#use icc and intel mkl lib
#os.environ["CC"]="icc"
#ext_modules = Extension("lu_fast",["lu_fast.pyx"],\
    #include_dirs=[numpy.get_include()],
    #extra_compile_args=['-O3'],
    ##extra_link_args=['-fopenmp'],
    ##libraries=['mkl_rt']
    #)
 
#use default compiler
ext_modules = Extension("lu_fast",["lu_fast.pyx"],\
    include_dirs=[numpy.get_include()],
    extra_compile_args=['-O3'],
    extra_link_args=['-llapack']
    )
 
setup(
    name= 'lu',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext_modules]
)
