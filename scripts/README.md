
## Results

Edit paths to output generated from Tools in `scripts/results.py`

```shell
python results.py
```

## Simulations

Edit paths in `scripts/simulations/generate.sh` for
* BWA, Samtools, SVsim, ART executables
* Simulation directory

Requires around 1TB of space for all simulations

```shell
bash simulations/generate.sh
# or SLURM
sbatch simulations/generate.sh
```
