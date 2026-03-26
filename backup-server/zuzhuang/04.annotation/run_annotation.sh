
#去除pep文件内的星号
sed 's/\*//g' sugarglider_row.pep > test.pep
##注释

perl 00.function.annotation.pl -cutf 20 -Interpro -KEGG -Swissprot -Trembl -Cog  -NR  -pep test.pep -nodes node2+node4+node10+node3 -ppn 16
