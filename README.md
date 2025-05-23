
# epsSMASH
A bioinformatic tool which predicts biosynthetic gene clusters associated with EPS production.


<img src="epssmash/outputs/html/images/bacteria_epssmash_logo.svg" alt="drawing" width="200"/>


**NOTE**: This is a work in progress.

## Installation

Clone epsSMASH to your local machine
`git clone https://github.com/AOHD/epsSMASH.git`

Make a conda environment using the YAML file in `epssmash/config/`

`conda env create -f epssmash/config/config.yaml`


Install custom-smash in the conda environment with pip

`pip install -e .`

The -e makes it so whenever you run custom-smash, any changes you've made in the codebase will be included (So this is only relevant for developers)

Download antismash databases (For Pfam annotation)

`download-antismash-databases`

Run epsSMASH

`epsSMASH --help-showall`

