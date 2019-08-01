
## Results

Edit paths to output generated from Tools in `scripts/results.py`

```shell
python3 results.py
```

## Simulations

Edit paths in `scripts/simulations.sh` for
* BWA, Samtools, SVsim, ART executables
* Simulation directory

Requires around 1TB of space for all simulations

```shell
bash simulations.sh
# or using SLURM
sbatch simulations.sh
```
