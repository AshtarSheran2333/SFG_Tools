#!/bin/bash

compiler_path=$(which ifx)

if [ -z "$compiler_path" ]; then
	echo "ifx (intel fortran compiler) not found in path"
	exit 1
else
	cmake -B build -DCMAKE_INSTALL_PREFIX=./ -DCMAKE_Fortran_COMPILER="$compiler_path"
	cmake --build build
	cmake --install build
fi
