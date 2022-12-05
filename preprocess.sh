


rm -r ../smoother_out/radicl.smoother_index

python3 python/main.py indexer init ../smoother_out/radicl ../smoother/Lister427.sizes -d 1000 --test

python3 python/main.py indexer anno ../smoother_out/radicl ../smoother/HGAP3_Tb427v10_merged_2021_06_21.gff3

gdb python3 -ex "run python/main.py indexer repl ../smoother_out/radicl ../smoother_in/anna.sort.test.PRE2 P10_Total"
#python3 python/main.py indexer repl ../smoother_out/radicl ../smoother_in/anna.sort.test.PRE2 P10_Total

#gdb python3 -ex "run python/main.py indexer repl ../smoother_out/hic2 ../smoother_in/claudia.pre1 P10_R1"
#python3 python/main.py indexer repl -q -m ../smoother_out/hic ../smoother_in/claudia.pre1 P10_R1


#gdb python3 -ex "run python/main.py indexer track ../smoother_out/hic ../smoother_in/coverage.tsv.sorted rna_seq"
python3 python/main.py indexer track ../smoother_out/radicl ../smoother_in/coverage.tsv.sorted rna_seq