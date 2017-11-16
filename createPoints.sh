# ./createPoints.sh points/ tests/

er() { echo "\$ $@" ; "$@" ; }

outDir=$1
inDir=$2

mkdir -p $outDir

for flagDir in $(ls $inDir); do
	flag=$flagDir
	
	case "$flag" in
		("L2CACHE") 
			value="L2 miss ratio"
			;;
		("L3")
			value="L3 bandwidth [MBytes/s]"
			;;
		("FLOPS_DP")
			value="DP MFLOP/s"
			;;
		(*)
			value="err"
			;;
	esac
	
	echo $flagDir
	echo $value
	
	for filename in $(ls $inDir/$flagDir); do
		filepath="$inDir/$flagDir/$filename"
		
		version=$(echo $filename | cut -d'-' -f3 | cut -d'.' -f1)
		
		x=$(echo $filename | cut -d'-' -f4 | cut -d'.' -f1)
		echo "./likParse.sh \"$filepath\" \"$value\""
		ys=$(./likParse.sh $filepath "$value")
		echo $ys
		
		i=0
		for y in $ys; do
			case "$i" in
				0)  op="INV"
					;;
				
				1)  op="RES"
					;;
			
				2)  op="SUM"
					;;
				
				*)  op="err"
					;;
			esac
			
			echo $x $y >> $outDir/$op-$flag-$version.dat
			i=$(($i+1))
		done
	done
done

for filename in $(ls $outDir); do
	sort -k1 -n $outDir/$filename > $outDir/$filename.tmp
	mv $outDir/$filename.tmp $outDir/$filename
done