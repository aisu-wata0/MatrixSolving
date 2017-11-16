# ./lik.sh "liklog.txt" CACHE "./invmat -r512 -i4"
likwid-perfctr -C 1  -f  -o "$1" -O -g $2 -m "$3"
