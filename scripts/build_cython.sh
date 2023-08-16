#!/bin/bash

rm -vf *python*so
sage -python <<EOF
from sage.misc.cython import cython
cython("ntl_helpers.pyx", create_local_so_file=True)
EOF