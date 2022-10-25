


rm -r ../smoother_out/test.smoother_index

python3 python/indexer_parsers.py indexer init ../smoother_out/test ../smoother/Lister427.sizes -d 1000 --test

python3 python/indexer_parsers.py indexer anno ../smoother_out/test ../smoother/HGAP3_Tb427v10_merged_2021_06_21.gff3

python3 python/indexer_parsers.py indexer repl ../smoother_out/test ../smoother_in/anna.sort.test.PRE2 P10_Total