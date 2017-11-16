# ./createPoints.sh points/ tests/* 
outDir=$1

mkdir -p $outDir

agrNum=0
for filename in "$@"; do
	agrNum=$(($agrNum+1))
	if [ $agrNum == "1" ]; then
		continue
	fi
	x=$(echo $filename | cut -d'-' -f4 | cut -d'.' -f1)
	ys=$(./missRatio.sh $filename)
	
	i=0
	for y in $ys; do
		case "$i" in
			0)  flag="INV"
				;;
				
			1)  flag="RES"
				;;
			
			2)  flag="SUM"
				;;
				
			*)  flag="err"
				;;
		esac
		
		touch $outDir/$flag.txt
		echo $x $y >> $outDir/$flag.txt
		i=$(($i+1))
	done

done
