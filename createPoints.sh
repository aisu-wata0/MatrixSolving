# ./createPoints.sh points/ tests/*
#                                 * = Every Flag (dir)
outDir=$1

mkdir -p $outDir

agrNum=0
for flagDir in "$@"; do 
	agrNum=$(($agrNum+1))
	if [ $agrNum == "1" ]; then
		continue
	fi
	
	flag=$(echo $flagDir | rev | cut -d'/' -f1 | rev)
	
	case "$flag" in
		"L2CACHE") 
			value="L2 miss ratio"
			;;
		"L3")
			value="L3 bandwidth [MBytes/s]"
			;;
		"FLOPS_DP")
			value="DP MFLOP/s"
			;;
		*)
			value="err"
			;;
	esac
	
	for filename in $(ls $flagDir); do
		x=$(echo $flagDir/$filename | cut -d'-' -f4 | cut -d'.' -f1)
		ys=$(./likParse.sh $flagDir/$filename $value)
		
		i=0
		for y in $ys; do
			case "$i" in
				0)  re="INV"
					;;
				
				1)  re="RES"
					;;
			
				2)  re="SUM"
					;;
				
				*)  re="err"
					;;
			esac
			touch $outDir/$flag.txt
			echo $x $y >> $outDir/$re-$flag-$version.dat
			i=$(($i+1))
		done
	done
done
