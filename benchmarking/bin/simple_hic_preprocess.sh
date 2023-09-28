





# align using bwa mem
bwa mem -t 8 -SP ${GENOME_NAME}.fna ${READ_FILE_1} ${READ_FILE_2} | \
    # load the sam fil from bwa mem
    # --drop-sam: drop the sam information
    # --min-mapq 0: do not filter by mapping quality
    # --add-columns mapq,XA: add columns for mapping quality and alternative alignments
    # --walks-policy mask: mask the walks (WW, XX, N*) in the pair_type column
    pairtools parse --drop-sam --min-mapq 0 --add-columns mapq,XA --walks-policy mask -c ${GENOME_NAME}.sizes | \
    # filter out the walks (WW, XX, N*)
    pairtools select '(pair_type!="WW") and (pair_type!="XX") and not wildcard_match(pair_type, "N*")' | \
    # filter out duplicates (sort then dedup)
    pairtools sort --nproc 8 | \
    pairtools dedup | \
    # output as a pair file
    pairtools split --output-pairs ${SAMPLE_NAME}.pairs.gz