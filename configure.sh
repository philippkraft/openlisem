#!/bin/bash
QWT_DIR=$(realpath ../qwt-?.?.?-ma)
echo "$QWT_DIR"
cmake -S . -B build/ \
	-DQWT_BUILD_DIR=$QWT_DIR \
	-DGDAL_INCLUDE_DIRS=$CONDA_PREFIX/include \
	-DGDAL_LIBRARIES=$CONDA_PREFIX/lib/libgdal.so
