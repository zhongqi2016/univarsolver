from distutils.core import setup, Extension
from Cython.Build import cythonize

ext1 = Extension(
    name="psl_bounds_cy",
    sources=["./FZCP/psl_bounds_cy.pyx"],
    language="c++"
)
setup(name="psl_bounds_cy", ext_modules=cythonize(ext1, language_level=3))

ext2 = Extension(
    name="psqe_bounds_cy",
    sources=["./FZCP/psqe_bounds_cy.pyx"],
    language="c++"
)
setup(name="psqe_bounds_cy", ext_modules=cythonize(ext2, language_level=3))
