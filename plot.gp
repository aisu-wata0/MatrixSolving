print command
set terminal png size 400,300 enhanced font "Helvetica,20"
print filepath
set output filepath
set logscale y 2
plot command