# ./plot.sh inDir/ outDir/
# ./plot.sh points/ graphs/
inDir=$1
outDir=$2

mkdir -p $outDir/

# declare a hash
declare -A commandHash
declare -A commandHash2
i=0
for filename in $(ls $inDir); do
	op=$(echo $filename | cut -d'-' -f1 | cut -d'.' -f1)
	flag=$(echo $filename | cut -d'-' -f2 | cut -d'.' -f1)
	version=$(echo $filename | cut -d'-' -f3 | cut -d'.' -f1)
	
	if [ "${commandHash[$op-$flag]}" != "" ]; then
		#commandHash[$op-$flag]="${commandHash[$op-$flag]}, \"$inDir/$filename\" w l title \"$version\""
		#titles="$titles '$inDir/$filename'"
		
		commandHash2[$op-$flag]="$inDir/$filename"
		title2[$op-$flag]="${version//\_/\\\_}"
		title2[$op-$flag]="${title2[$op-$flag]//\\\_lik*/}"
		echo ${title2[$op-$flag]}
	else
		#commandHash[$op-$flag]="\"$inDir/$filename\" w l title \"$version\""
		#titles="'$inDir/$filename'"
		commandHash[$op-$flag]="$inDir/$filename"
		title[$op-$flag]="${version//\_/\\\_}"
		title[$op-$flag]="${title[$op-$flag]//\\\_lik*/}"
		echo ${title[$op-$flag]}
	fi
	
	echo ${commandHash[$op-$flag]}
	
	i=$(($i+1))
done

for index in ${!commandHash[@]}; do
	command="${commandHash[$index]}"

	op=$(echo $index | cut -d'-' -f1 | cut -d'.' -f1)
	flag=$(echo $index | cut -d'-' -f2 | cut -d'.' -f1)
	
	#gnuplot -e "command='${command}'; filepath='$outDir/$op-$flag.png'" plot.gp
	
	gnuplot -e "datfile1='${commandHash[$index]}'; datfile2='${commandHash2[$index]}'; title1='${title[$index]}'; title2='${title2[$index]}'; filepath='$outDir/$op-$flag.png'" plot.gp
done 
