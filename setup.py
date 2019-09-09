from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

particleinterpolator_extension = Extension(
    name="particleinterpolator",
    sources=["particleinterpolator.pyx"],
    libraries=["interpolate"],
    library_dirs=["lib"],
    include_dirs=["lib"]
)
setup(
    name="particleinterpolator",
    ext_modules=cythonize([particleinterpolator_extension])
)
