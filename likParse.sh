# ./likParse.sh "tests/liklog-L2CACHE-naive_lik-256.txt" "L2 miss ratio"
# ./likParse.sh "tests/liklog-L3-naive_lik-256.txt" "L3 bandwidth [MBytes/s]"
# ./likParse.sh "tests/liklog-FLOPS_DP-naive_lik-256.txt" "DP MFLOP/s"
cat "$1" | 
grep "$2" | 
sed s/,/\\n/g | 
sed s/"$2"//g | 
sed s/--//g | 
tr '\n\n' ' '
