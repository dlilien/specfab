#! /bin/sh
#
# entrypoint.sh
# Copyright (C) 2021 dlilien <dlilien@hozideh>
#
# Distributed under terms of the MIT license.
#


#Prereqs for python
yum install -y lapack-devel fftw-devel
export LD_LIBRARY_PATH=`$PWD`/specfab/src:$LD_LIBRARY_PATH 
export LD_LIBRARY_PATH=/opt/rh/devtoolset-9/root/lib/gcc/x86_64-redhat-linux/9:$LD_LIBRARY_PATH 
make lib

for PYBIN in /opt/python/cp37-cp*/bin; do
    "${PYBIN}/pip" install numpy==1.19.0 cython;
done

/opt/python/cp37-cp37m/bin/python /entrypoint.py
