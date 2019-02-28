snakemake -c  "qsub -cwd -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -l m_mem_free={cluster.m_mem_free} -pe smp {threads}" \
-s src/rules/sf_nb_predict_pipeline.py allCgiScores \
--cluster-config configs/cluster.yaml \
--rerun-incomplete -k -j 2 --latency-wait 90 --greediness 0.7
