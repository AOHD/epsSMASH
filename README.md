
# epsSMASH
A bioinformatic tool which predicts biosynthetic gene clusters associated with EPS production.


<img src="epssmash/outputs/html/images/bacteria_epssmash_logo.svg" alt="drawing" width="200"/>


## Installation

Clone epsSMASH to your local machine

`git clone https://github.com/AOHD/epsSMASH.git`

Make a conda environment using the YAML file in `epssmash/config/`

`conda env create -f epssmash/config/config.yaml`

`conda activate antismash_dependencies`

Install epsSMASH in the conda environment with pip

`pip install .`

Download Pfam and clusterblast databases

`download-epssmash-databases`

Run epsSMASH

`epsSMASH --help-showall`

