from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "glue",
    ext_modules = cythonize('*.pyx'),
    language="c++"
)
