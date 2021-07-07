#!/bin/sh


# Install non-python dependencies
yum install -y lapack-devel fftw-devel
export LD_LIBRARY_PATH=/github/workspace/specfab/src:$LD_LIBRARY_PATH 
export LD_LIBRARY_PATH=/opt/rh/devtoolset-9/root/lib/gcc/x86_64-redhat-linux/9:$LD_LIBRARY_PATH 
export FC=/opt/rh/devtoolset-9/root/usr/bin/gfortran
export PLAT=manylinux2014_x86_64

# Build specfab
make lib

# Do a normal build
for PYBIN in /opt/python/cp38-cp*/bin; do
    "${PYBIN}/pip" install numpy==1.19.0 cython;
    "${PYBIN}/pip" wheel --no-deps -w /github/workspace/wheelhouse/ .;
done

# Make the wheels into manylinux
ls wheelhouse/*.whl
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w /github/workspace/wheelhouse/;
done

python -m twine upload /github/workspace/wheelhouse/*manylinux*.whl
