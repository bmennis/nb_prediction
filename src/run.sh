snakemake -s rules/sf_cov_in_bam.py intersect_wxs_capture \
--rerun-incomplete --drmaa -k -j 2 --latency-wait 40 --greediness 0.7
