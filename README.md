## Simulations of synaptonemal complex partitioning:

1. `copos.py` -- Simulation of continuous signal diffusion and accumulation from crossovers
2. `ssd.py` -- Simulation of discrete SC central elements

`copos.py` -- best run as Jupyter notebook (see copos.ipynb)
Parameters: chromosome length and crossover positions.

This simulation illustrates our prediction of short and long arm identity. A single chromosome has signals emanating in both directions from all specified crossovers. After a certain amount of time has passed, new signals stop emerging and the remaining signals equalize out. The resulting signal height is used to determine, for each crossover, which direction is "short arm" and which direction is "long arm".


`ssd.py` -- can run from notebook or command line
Parameters given in file `ssd.conf`. Each line of `ssd.conf` represents one chromosome that is to be simulated. For each chromosome (each line of the text file), there are at least 2 numbers. The first number is the chromosome length; the second (and more, if present) is the crossover position(s). To specify a non-crossover chromosome, use a crossover position that is greater than the total chromosome length. You may specify as many chromosomes and crossovers as you wish.

The simulation runs as follows: thousands of individual molecules of SYP-1 are created with the following properties:
 1. position (chromosome #, coordinate ∷ -1 for "not localized to a chromosome" e.g. nucleoplasmic)
 2. phosphorylation state (True=phosphorylated, False=non-phosphorylated)
 3. PLK-2-bound state (True=bound, False=unbound) ∷ only phospho-SYP-1 can bind it
Also, a chromosome axis ("sbs" for PLK-2 substrate) (1 location per chromosome position, unmoving) is simulated that can be either phosphorylated or not at each position

at each timestep,
- the SYP-1 molecules move randomly -- they are slowed down by phosphorylated SBSTR8
- those near the crossover have a chance to be phosphorylated
- phosphorylated SYP-1s near the crossover have a chance to pick up PLK-2
- PLK-2 bound to SYP-1 can phosphorylate the substrate if needed
