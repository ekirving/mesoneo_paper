# Ancestry-stratified allele frequency trajectories
This repository contains code for the ancestry stratified selection analyses from 
[The Selection Landscape and Genetic Legacy of Ancient Eurasians](https://).

![Figure 4](./figure/Figure_4.png?raw=true)

If you reuse any of this code then please cite the preprint:
> Irving-Pease, E.K.&ast;, Refoyo-Martínez, A.&ast;, Ingason, A.&ast;, Pearson, A.&ast;, Fischer, A.&ast;, Barrie, 
> W.&ast;, Sjögren, K.-G., Halgren, A.S., Macleod, R., Demeter, F., Henriksen, R.A., Vimala, T., McColl, H., Vaughn, A., 
> Stern, A.J., Speidel, L., Scorrano, G., Ramsøe, A., Schork, A.J., Rosengren, A., Zhao, L., Kristiansen, K., Sudmant, 
> P.H., Lawson, D.J., Durbin, R., Korneliussen, T., Werge, T., Allentoft, M.E., Sikora, M., Nielsen, R., Racimo, F., 
> Willerslev, E., 2022. The Selection Landscape and Genetic Legacy of Ancient Eurasians. *bioRxiv* 2022.09.22.509027. 
> https://www.biorxiv.org/content/10.1101/2022.09.22.509027v1

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
