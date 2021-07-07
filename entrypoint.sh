#!/bin/sh


# Install non-python dependencies
yum install -y lapack-devel fftw-devel
export LD_LIBRARY_PATH=`$PWD`/specfab/src:$LD_LIBRARY_PATH 
export LD_LIBRARY_PATH=/opt/rh/devtoolset-9/root/lib/gcc/x86_64-redhat-linux/9:$LD_LIBRARY_PATH 
export FC=/opt/rh/devtoolset-9/root/usr/bin/gfortran

# Build specfab
make lib

# Do a normal build
for PYBIN in /opt/python/cp38-cp*/bin; do
    "${PYBIN}/pip" install numpy==1.19.0 cython;
    "${PYBIN}/pip" wheel --no-deps -w wheelhouse/ .;
done

# Make the wheels into manylinux
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w `$PWD`/wheelhouse/;
done
