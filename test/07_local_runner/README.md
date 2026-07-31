# Sanity check 2: staged runs, run and resumed

Two-stage runs (the second starting from the first's final structure) driven by
both runners. Deliberately tiny — a 10x10 lattice over 6 temperature steps — so
the whole thing takes a couple of seconds.

    python 00_generate_inputs.py
    python 01_run_local.py

`01_run_local.py` prints a line per check, in seven groups.

Groups 1–5 drive `runners.local.run_manifest`:

- a clean run completes both stages of every job;
- a second run skips everything, since each stage wrote its final structure;
- after an interruption is simulated (final structure and the last checkpoints
  deleted), `resume=True` restarts from the last surviving checkpoint and writes
  into the stage's own output directory;
- a second interruption resumes just as cleanly, the offsets accumulating;
- both stages of one job can be interrupted and continued in a single attempt
  without their temporary param files colliding.

Groups 6–7 drive the SLURM path without a cluster, by invoking
`frusa_lattice_mc.runners.task` once per manifest line exactly as an array task
would. One command starts a job that has no checkpoints, continues one that was
interrupted, and no-ops once everything is finished — and the script
`slurm.generate_self_resuming_script` issues one `srun` per array task
rather than one per stage, so it never needs regenerating.

Throughout, the checks enforce the point of the continuation machinery: however
many times a run is interrupted, and whichever runner picks it up, its
checkpoints stay a single gapless `structure_0…structure_{Nt-1}` series in one
directory, with no sibling `structures_*` folder. Continuation *inputs* are
versioned per attempt under `input/staged/_continued/<timestamp>_task<n>/`;
outputs are not.

Note that groups 3–6 delete checkpoints to fake an interruption. Deleting them
*below* a previous attempt's offset — which a real interruption never does —
makes the runner decline to resume and re-run the stage whole, so the fixtures
deliberately keep more each time.

Re-run `00_generate_inputs.py` to start over; it wipes `input/` and `data/`.
