# ./plot.sh inDir/ outDir/
# ./plot.sh points/ graphs/
inDir=$1
outDir=$2

mkdir -p $outDir/

# declare a hash
declare -A commandHash
i=0
for filename in $(ls $inDir); do
	op=$(echo $filename | cut -d'-' -f1 | cut -d'.' -f1)
	flag=$(echo $filename | cut -d'-' -f2 | cut -d'.' -f1)
	version=$(echo $filename | cut -d'-' -f3 | cut -d'.' -f1)
	
	if [ "${commandHash[$op-$flag]}" != "" ]; then
		commandHash[$op-$flag]="${commandHash[$op-$flag]}, \"$inDir/$filename\" w l title \"$version\""
		titles="$titles '$inDir/$filename'"
	else
		commandHash[$op-$flag]="\"$inDir/$filename\" w l title \"$version\""
		titles="'$inDir/$filename'"
	fi
	
	echo ${commandHash[$op-$flag]}
	
	i=$(($i+1))
done

for index in ${!commandHash[@]}; do
	command="${commandHash[$index]}"

	op=$(echo $index | cut -d'-' -f1 | cut -d'.' -f1)
	flag=$(echo $index | cut -d'-' -f2 | cut -d'.' -f1)
	
	#gnuplot -c plot.gp "${command}" $outDir/$op-$flag.png
	
	gnuplot -e "command='${command}'; filepath='$outDir/$op-$flag.png'" plot.gp
	gnuplot -e "command='${command}'; filepath='$outDir/$op-$flag.png'" plot.gp
done 