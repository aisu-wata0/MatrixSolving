#!/bin/bash

#32 33 64 65 128 129 256 257 512 1000 2000
er() { echo; echo "\$ $@" ; "$@" ; }

testdir=tests
mkdir -p ${testdir}

size=512

mkdir -p ${testdir}/tile

# for ((i=2; i<=8; i = i+1)); do
# 	export GCC_ARGS=" -D L0=32 -D L1M=${i}"
#
# 	er g++  "src/main.cpp" -std=c++11 -Wno-comment  -Wno-sign-compare \
# 	-O3  -mavx -march=native -Wall -Wno-unused-variable -DNDEBUG  \
# 	$GCC_ARGS -o "./invmatL0-$i" -I. -I./include
#
# 	args="-r$size -i10 -o${testdir}/tile/out_L1-${i}.2.txt"
# 	er ./invmatL0-$i $args
# 	er rm ./invmatL0-$i
# done

mkdir -p ${testdir}/unr

for iunr in 2 4; do
	for junr in 2 4 6 8; do
		export GCC_ARGS=" -D IUNRLL=${iunr} -D JUNRLL=${junr}"

		bin="invmatL0_iunr-${iunr}_junr-${junr}"

		er g++  "src/main.cpp" -std=c++11 -Wno-comment  -Wno-sign-compare \
		-O3  -mavx -march=native -Wall -Wno-unused-variable -DNDEBUG  \
		$GCC_ARGS -o ./$bin -I. -I./include

		args="-r$size -i20 -o${testdir}/unr/out_iunr-${iunr}_junr-${junr}.txt"
		er ./$bin $args
		er rm ./$bin
	done
done
