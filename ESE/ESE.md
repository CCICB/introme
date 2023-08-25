main
- loop through each variant
- header set as the 6 splicing enhancer motifs
- skip insertions longer than 30bp
- should we skip deletions longer than Xbp? or is it just cus it takes too long right now
- special case for SNV (or maybe case is not unique)
- get score for each of 6 splicing enhancers by diffing SEQ_SCAN of ref vs alt.
    - for fwd and/or rev, for whole length of seq (bar first and last)


    26 mid.large.strand