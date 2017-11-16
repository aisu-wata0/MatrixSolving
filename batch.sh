#!/bin/bash
# ./batch.sh tests/
# For each size, each version, each flag, in outDir create a directory for each flag, then create logs for each run inside the directory

er() { echo "\$ $@" ; "$@" ; }

outDir=$1
#"32 33 64 65 128 129 256 257 512 1000 2000"
sizes="32 33 64 65 128 129 256 257 512 1000 2000"
#"naive_lik" "mult_block_likwid"
versions="naive_lik mult_block_likwid"

for version in $versions; do
	er git checkout $version
	er make rebuild
	er mv invmat invmat_${version}
	
	for size in $sizes; do
		for LIK_FLAG in "L2CACHE" "FLOPS_DP" "L3"; do
			args="-r$size -i10 -olog-time-${version}-${size}.txt"
			mkdir -p $outDir/$LIK_FLAG
			er ./lik.sh "$outDir/${LIK_FLAG}/liklog-${LIK_FLAG}-${version}-${size}.txt" ${LIK_FLAG} "./invmat_${version} $args"
		done
	done
done

#    Banda de Memória: utilizar o grupo MEM do LIKWID, e apresentar o resultado de "Memory bandwidth [MBytes/s]"; Caso não tenha o grupo MEM, utilize o grupo L3.
#    Cache miss L1: utilizar o grupo CACHE do LIKWID, e apresentar o resultado de "data cache miss ratio". Caso não tenha o cache miss da L1, utilize o cache miss da L2.
#    Operações aritméticas: utilizar o grupo FLOPS_DP ou FLOPS_AVX do LIKWID, e apresentar o resultado de "MFLOP/s"
