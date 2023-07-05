#!/bin/bash
DIRNAME=$(dirname $0)
echo $DIRNAME
pushd $DIRNAME/build
make
popd
