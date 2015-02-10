from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy,os
os.environ["CC"]="icc"

ext_modules = Extension("lu_fast",["lu_fast.pyx"],\
    include_dirs=[numpy.get_include()],
    extra_compile_args=['-O3'],
    extra_link_args=['-fopenmp'],
    #extra_link_args=['-llapack']
    #libraries=['mkl_intel_lp64','mkl_intel_thread','mkl_core']
    libraries=['mkl_rt']
    )
 
setup(
    name= 'lu',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext_modules]
)
