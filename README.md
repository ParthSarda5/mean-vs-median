# Sample Mean vs. Sample Median: Efficiency, Robustness and a Monte Carlo Study

MT3194 *Statistical Inference* group project · **Group 11**: Bhumika, Parth, Satyam, Sriharipriya, Vishwabharathi

> How well do the sample mean and the sample median estimate the centre of a distribution when the data are not normal?
> We derive their variances and compare them with the Cramér–Rao bound, study robustness through influence functions and breakdown points,
> run a Monte Carlo study under heavy tails and outliers, and check the predictions on class data.

## Repository layout

```
mean-vs-median/
├── README.md                  this file
├── requirements.txt.          Python dependencies for the simulation.
├── theory/
├──  robustness/
├──  application/     
├── simulation/
├── report/
│   ├── README.md              design, results, how to run, method notes
│   ├── run_all.py             runs every experiment, then figures and tables
│   ├── mvm/                   the package: distributions, estimators, engine, theory, experiments, plotting, tables
│   ├── tests/                 35 unit and statistical tests
│   └── results/               data/ (CSV), figures/ (PDF + PNG), summary_*.md
├── report/                     .pdf, .tex file 



## Running the simulation

```bash
pip install -r requirements.txt
cd simulation
python -m pytest -q
python run_all.py --quick
```



See the [simulation README](simulation/README.md) for the study design and detailed results.
