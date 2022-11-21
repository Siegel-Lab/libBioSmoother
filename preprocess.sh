


rm -r ../smoother_out/hic.smoother_index

python3 python/indexer_parser.py indexer init ../smoother_out/hic ../smoother/Lister427.sizes -d 100 --test

python3 python/indexer_parser.py indexer anno ../smoother_out/hic ../smoother/HGAP3_Tb427v10_merged_2021_06_21.gff3

#python3 python/indexer_parser.py indexer repl ../smoother_out/hic ../smoother_in/anna.sort.test.PRE2 P10_Total
python3 python/indexer_parser.py indexer repl ../smoother_out/hic ../smoother_in/claudia.pre1 P10_R1