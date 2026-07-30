# Sanity check 2: staged local runs

Two-stage runs (the second resuming from the first) driven by
`runners.local.run_manifest`. Deliberately tiny — a 10x10 lattice over 6
temperature steps — so the whole thing takes a couple of seconds.

    python 00_generate_inputs.py
    python 01_run_local.py

`01_run_local.py` asserts three behaviours and prints a line per check:

- a clean run completes both stages of every job;
- a second run skips everything, since each stage wrote its final structure;
- after an interruption is simulated (final structure and the last checkpoints
  deleted), `resume=True` restarts from the last surviving checkpoint and puts
  the final structure back where the stage promised it.

Re-run `00_generate_inputs.py` to start over; it wipes `input/` and `data/`.
