# Benchmark tooling

This directory contains simulation and validation tooling that is not required
for a normal PolyploidHLA sample run.

- `scripts/`: simulation generation, scoring, and result summarization.
- `slurm/`: local launchers and Slurm jobs for reproducible benchmark matrices.
- Generated reads, truth, logs, runs, and metrics belong in the sibling
  `PolyploidHLA_simulation/` directory, not in this repository.

Run commands from the repository root so paths shown in
[`../docs/SIMULATION_EXPERIMENT_GUIDE.md`](../docs/SIMULATION_EXPERIMENT_GUIDE.md)
resolve consistently.
