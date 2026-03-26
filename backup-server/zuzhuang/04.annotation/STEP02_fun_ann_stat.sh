cd  /data01/sunxuebo/project/zuzhuang/04.annotation 
######interpro#####
cat  /data01/sunxuebo/project/zuzhuang/04.annotation/interpro/test.pep.*.fa.tsv  > /data01/sunxuebo/project/zuzhuang/04.annotation/interpro/test.pep.interpro 
perl /public/home/humingliang/software/function_annatation/interpro.change.pl  /public/home/humingliang/software/function_annatation/interpro2go.txt   /data01/sunxuebo/project/zuzhuang/04.annotation/interpro/test.pep.interpro 
cd /data01/sunxuebo/project/zuzhuang/04.annotation && ln -s /data01/sunxuebo/project/zuzhuang/04.annotation/interpro/test.pep.interpro.GO  /data01/sunxuebo/project/zuzhuang/04.annotation/test.pep.GO.result
cd /data01/sunxuebo/project/zuzhuang/04.annotation && ln -s /data01/sunxuebo/project/zuzhuang/04.annotation/interpro/test.pep.interpro.IPR  /data01/sunxuebo/project/zuzhuang/04.annotation/test.pep.InterPro.result

#######kegg####
cat  /data01/sunxuebo/project/zuzhuang/04.annotation/kegg/test.pep.*.fa.kegg.blast  > /data01/sunxuebo/project/zuzhuang/04.annotation/kegg/test.pep.kegg.tmp1 
perl /public/home/humingliang/software/function_annatation/blast_parser.pl  -nohead  /data01/sunxuebo/project/zuzhuang/04.annotation/kegg/test.pep.kegg.tmp1 > /data01/sunxuebo/project/zuzhuang/04.annotation/kegg/test.pep.kegg.tmp2
perl /public/home/humingliang/software/function_annatation/blast_parser.pl  -nohead -tophit 1 -topmatch 1   /data01/sunxuebo/project/zuzhuang/04.annotation/kegg/test.pep.kegg.tmp1 > /data01/sunxuebo/project/zuzhuang/04.annotation/kegg/test.pep.kegg
perl /public/home/humingliang/software/function_annatation/kegg.change.pl /data01/sunxuebo/project/zuzhuang/04.annotation/kegg/test.pep.kegg
cd /data01/sunxuebo/project/zuzhuang/04.annotation && ln -s /data01/sunxuebo/project/zuzhuang/04.annotation/kegg/test.pep.kegg.change  /data01/sunxuebo/project/zuzhuang/04.annotation/test.pep.KEGG.result

#######trembl####
cat  /data01/sunxuebo/project/zuzhuang/04.annotation/trembl/test.pep.*.fa.trembl.blast  > /data01/sunxuebo/project/zuzhuang/04.annotation/trembl/test.pep.trembl.tmp1 
perl /public/home/humingliang/software/function_annatation/blast_parser.pl  -nohead  /data01/sunxuebo/project/zuzhuang/04.annotation/trembl/test.pep.trembl.tmp1 > /data01/sunxuebo/project/zuzhuang/04.annotation/trembl/test.pep.trembl.tmp2
perl /public/home/humingliang/software/function_annatation/blast_parser.pl  -nohead -tophit 1 -topmatch 1   /data01/sunxuebo/project/zuzhuang/04.annotation/trembl/test.pep.trembl.tmp1 > /data01/sunxuebo/project/zuzhuang/04.annotation/trembl/test.pep.trembl
cd /data01/sunxuebo/project/zuzhuang/04.annotation && ln -s /data01/sunxuebo/project/zuzhuang/04.annotation/trembl/test.pep.trembl  /data01/sunxuebo/project/zuzhuang/04.annotation/test.pep.TrEMBL.result

######swissprot#####
cat  /data01/sunxuebo/project/zuzhuang/04.annotation/swissprot/test.pep.*.fa.swissprot.blast  > /data01/sunxuebo/project/zuzhuang/04.annotation/swissprot/test.pep.swissprot.tmp1 
perl /public/home/humingliang/software/function_annatation/blast_parser.pl  -nohead  /data01/sunxuebo/project/zuzhuang/04.annotation/swissprot/test.pep.swissprot.tmp1 > /data01/sunxuebo/project/zuzhuang/04.annotation/swissprot/test.pep.swissprot.tmp2
perl /public/home/humingliang/software/function_annatation/blast_parser.pl  -nohead -tophit 1 -topmatch 1   /data01/sunxuebo/project/zuzhuang/04.annotation/swissprot/test.pep.swissprot.tmp1 > /data01/sunxuebo/project/zuzhuang/04.annotation/swissprot/test.pep.swissprot
cd /data01/sunxuebo/project/zuzhuang/04.annotation && ln -s /data01/sunxuebo/project/zuzhuang/04.annotation/swissprot/test.pep.swissprot  /data01/sunxuebo/project/zuzhuang/04.annotation/test.pep.Swissprot.result

#######cog####
cat  /data01/sunxuebo/project/zuzhuang/04.annotation/cog/test.pep.*.fa.cog.blast  > /data01/sunxuebo/project/zuzhuang/04.annotation/cog/test.pep.cog.tmp1
perl /public/home/humingliang/software/function_annatation/blast_parser.pl  -nohead  /data01/sunxuebo/project/zuzhuang/04.annotation/cog/test.pep.cog.tmp1 > /data01/sunxuebo/project/zuzhuang/04.annotation/cog/test.pep.cog.tmp2
perl /public/home/humingliang/software/function_annatation/blast_parser.pl  -nohead -tophit 1 -topmatch 1   /data01/sunxuebo/project/zuzhuang/04.annotation/cog/test.pep.cog.tmp1 > /data01/sunxuebo/project/zuzhuang/04.annotation/cog/test.pep.cog
cd /data01/sunxuebo/project/zuzhuang/04.annotation && ln -s /data01/sunxuebo/project/zuzhuang/04.annotation/cog/test.pep.cog  /data01/sunxuebo/project/zuzhuang/04.annotation/test.pep.Cog.result

#######nr####
cat  /data01/sunxuebo/project/zuzhuang/04.annotation/nr/test.pep.*.fa.nr.blast  > /data01/sunxuebo/project/zuzhuang/04.annotation/nr/test.pep.nr.tmp1
perl /public/home/humingliang/software/function_annatation/blast_parser.pl  -nohead  /data01/sunxuebo/project/zuzhuang/04.annotation/nr/test.pep.nr.tmp1 > /data01/sunxuebo/project/zuzhuang/04.annotation/nr/test.pep.nr.tmp2
perl /public/home/humingliang/software/function_annatation/blast_parser.pl  -nohead -tophit 1 -topmatch 1   /data01/sunxuebo/project/zuzhuang/04.annotation/nr/test.pep.nr.tmp1 > /data01/sunxuebo/project/zuzhuang/04.annotation/nr/test.pep.nr
cd /data01/sunxuebo/project/zuzhuang/04.annotation && ln -s /data01/sunxuebo/project/zuzhuang/04.annotation/nr/test.pep.nr  /data01/sunxuebo/project/zuzhuang/04.annotation/test.pep.nr.result

perl /public/home/humingliang/software/function_annatation/fun_ann_sta.pl  test.pep
perl /public/home/humingliang/software/function_annatation/stat_fun_ann.pl test.pep.all.function
