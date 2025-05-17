from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        name="calc_dpred",
        sources=["calc_dpred.pyx"],
        include_dirs=[np.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        extra_compile_args=["-O3"],  
        language="c",  
    )
]

setup(
    name="calc_dpred",
    ext_modules=cythonize(
        extensions,
        language_level=3,  
        compiler_directives={
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True
        }
    ),
)
