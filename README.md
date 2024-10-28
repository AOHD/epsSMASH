
# epsSMASH
A bioinformatic tool which predicts biosynthetic gene clusters associated with EPS production.

**NOTE**: This is a work in progress.

## Installation

Clone custom-smash repo to your local machine
`git clone https://github.com/AOHD/epsSMASH.git`

Make a conda environment with the dependencies

`conda create -n antismash_deps`

`conda activate antismash_deps`

`conda install hmmer2 hmmer diamond fasttree prodigal blast muscle glimmerhmm python=3.11 -c bioconda` 


Install custom-smash in the conda environment with pip

`pip install -e .`

The -e makes it so whenever you run custom-smash, any changes you've made in the codebase will be included

Download antismash databases

`download-antismash-databases`

Run custom-smash

`custom-smash --help-showall`

![epsSMASHlogo](customsmash/outputs/html/images/bacteria_epssmash_logo.svg?raw=true "epsSMASH")