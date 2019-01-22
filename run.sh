snakemake -s src/rules/sf_nb_predict_pipeline.py annotateDbnsfp \
--rerun-incomplete -k -j 2 --latency-wait 40 --greediness 0.7
