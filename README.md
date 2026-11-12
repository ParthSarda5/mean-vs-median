# Mean vs. Median

A statistical inference group project comparing the sample mean and median. The completed Monte Carlo study, its code, tests, and results are in [simulation/](simulation/README.md).

## Repository layout

- `simulation/` — simulation code, tests, figures, and results.
- `application/`, `robustness/`, `theory/`, and `report/` — placeholders for future work.
- `requirements.txt` — Python dependencies for the simulation.

## Running the simulation

```bash
pip install -r requirements.txt
cd simulation
python -m pytest -q
python run_all.py --quick
```

See the [simulation README](simulation/README.md) for the study design and detailed results.
