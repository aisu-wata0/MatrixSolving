#!/bin/bash
#32 33 64 65 128 129 256 257 512 1000 2000
er() { echo "\$ $@" ; "$@" ; }
sizes="32 33 64 65 128 129 256 257 512 1000 2000"

for size in $sizes; do
	for version in "naive_lik" "mult_block_likwid"; do
		er git checkout $version
		er make rebuild
		er mv invmat invmat_${version}
		for LIK_FLAG in "CACHE" "FLOPS_DP" "MEM"; do
			args="-r$size -i10 -olog_${version}_${size}.txt"
			er ./lik.sh "loglik_${LIK_FLAG}_${version}_${size}.txt" ${LIK_FLAG} "./invmat $args"
		done
	done
done

#    Banda de Memória: utilizar o grupo MEM do LIKWID, e apresentar o resultado de "Memory bandwidth [MBytes/s]"; Caso não tenha o grupo MEM, utilize o grupo L3.
#    Cache miss L1: utilizar o grupo CACHE do LIKWID, e apresentar o resultado de "data cache miss ratio". Caso não tenha o cache miss da L1, utilize o cache miss da L2.
#    Operações aritméticas: utilizar o grupo FLOPS_DP ou FLOPS_AVX do LIKWID, e apresentar o resultado de "MFLOP/s"
