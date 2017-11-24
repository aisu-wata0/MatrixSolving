set terminal png size 1080,720 enhanced font "Helvetica,23"
print filepath
set output filepath
set logscale y 2

set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

plot datfile1 w l title title1, datfile2 w l title title2
