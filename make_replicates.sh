#!/bin/bash

for i in {1..3000}; do
    echo Making replicate$i
    sed "s/\REPLICATE/$i/" example.pbs > replicate_$i.pbs
    sbatch replicate_$i.pbs
done



