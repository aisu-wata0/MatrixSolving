#!/bin/bash

#32 33 64 65 128 129 256 257 512 1000 2000
er() { echo "\$ $@" ; "$@" ; }

for size in 512; do
	for version in "Master" "OPTM" "Tri" "PTri"; do
		for ((i=0; i<2; i++)); do
			er ./MatrixSolving${version}.exe -r$size -i3 -oout${version}${i}.txt
		done
	done
done
