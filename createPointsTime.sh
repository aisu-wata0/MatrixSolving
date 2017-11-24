# ./createPoints.sh points/ tests/

er() { echo "\$ $@" ; "$@" ; }

outDir=$1
inDir=$2

mkdir -p $outDir

flag="time"
flagDir=$flag

for filename in $(ls $inDir/$flagDir); do
	filepath="$inDir/$flagDir/$filename"
	
	version=$(echo $filename | cut -d'-' -f3 | cut -d'.' -f1)
	
	x=$(echo $filename | cut -d'-' -f4 | cut -d'.' -f1)
	echo $filepath
	ys=$(cat $filepath)
	echo $ys
	
	i=0
	for y in $ys; do
		op="err"
		case "$i" in
			0)  
				;;
			
			1)  op="INV"
				;;
		
			2)  op="RES"
				;;
			
			*)  
				;;
		esac
		if [ $op != "err" ]; then
			echo $x $y >> $outDir/$op-$flag-$version.dat
		fi
		
		i=$(($i+1))
	done
done

for filename in $(ls $outDir); do
	flag=$(echo $filename | cut -d'-' -f2 | cut -d'.' -f1)
	if [ $flag != "time" ]; then
		continue
	fi
	echo $filename
	sort -k1 -n $outDir/$filename > $outDir/$filename.tmp
	mv $outDir/$filename.tmp $outDir/$filename
done