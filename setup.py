from numpy.distutils.core import Extension
from numpy.distutils.core import setup


extension_module = Extension(name='specfab.specfabpy',
                             sources=['specfab/src/specfabpy.f90'],
                             library_dirs=['specfab/src'],
                             libraries=['specfab'],
                             extra_f90_compile_args=['-ffree-line-length-none', '-mcmodel=medium'])

if __name__ == '__main__':
    setup(name='specfab',
          version='1.0.1',
          description='This is a python version of Specfab',
          ext_modules=[extension_module],
          packages=['specfab'])
