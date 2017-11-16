cat $1 | grep "L2 " | sed s/,/\\n/g | grep -A 1 "L2 miss ratio" | sed s/"L2 miss ratio"//g | sed s/--//g | tr '\n\n' ' '

