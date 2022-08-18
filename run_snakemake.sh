## -pe orte 2
## --keep-going
snakemake --snakefile Snakefile \
          --cluster-config config/cluster.json \
          --cluster "qsub -N '{rule}.{wildcards}' -o 'logs/{rule}/{rule}.{wildcards}.o' -e 'logs/{rule}/{rule}.{wildcards}.o' -V -S /bin/bash -cwd" \
          --jobs 200 \
	  --latency-wait 60 \
	  --rerun-incomplete \
	  --keep-going \
	  --printshellcmds
