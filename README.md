# Ancestry-stratified allele frequency trajectories
This repository contains code for the ancestry stratified selection analyses from 
[The Selection Landscape and Genetic Legacy of Ancient Eurasians](
https://doi.org/10.1038/s41586-023-06705-1).

![Figure 4](./figure/Figure_4.png?raw=true)

If you reuse any of this code then please cite the paper:
> Irving-Pease, E.K.&ast;, Refoyo-Martínez, A.&ast;, Barrie, W.&ast; et al. The selection landscape and genetic legacy 
> of ancient Eurasians. Nature 625, 312–320 (2024). https://doi.org/10.1038/s41586-023-06705-1

## Installation
Download the code: 
```bash
git clone git@github.com:ekirving/mesoneo_paper.git && cd mesoneo_paper/
```

The easiest way to install all the dependencies is with the [conda package manager](https://docs.conda.io/en/latest/).

```bash
conda env create --name mesoneo --file environment.yaml
```

Then activate the environment:
```bash
conda activate mesoneo
```

## Running the code

This project contains rules for running tens-of-thousands of selection tests.

```bash
# run all the selection tests using 1000G data  
./run_modern.sh
```

```bash
# run all the selection tests using aDNA time-series
./run_ancient.sh
```

```bash
# run all the selection tests using the ancestry-stratified aDNA time-series
./run_ancestries.sh
```

## Author

Evan K. Irving-Pease, [GLOBE Institute](https://globe.ku.dk/), University of Copenhagen 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
