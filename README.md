
# epsSMASH - the extracellular polymeric substance Secondary Metabolite Analysis SHell

### Installation via Github repository
To install epsSMASH via Github, you first need to clone the repository to your local machine. Then you need to install all the dependencies which epsSMASH relies on, which can be found in a YAML file at `epssmash/config/`. Once the dependencies are installed, epsSMASH can be installed using pip. Databases needed for Clusterblast and Pfam annotation of clusters can then be downloaded using the download-epssmash-databases command. The commands needed to perform the complete install are listed below:
```
git clone https://github.com/AOHD/epsSMASH.git # Clone repo

conda env create -f epssmash/config/config.yaml # Create conda environment with dependencies

conda activate antismash_dependencies # Activate conda environment

pip install . # Install epsSMASH with pip

download-epssmash-databases # Download databases needed for Clusterblast and Pfam annotation

epsSMASH --help-showall # Test installation
```

### Using Docker
[epsSMASH is available on the Docker Hub](https://hub.docker.com/r/antismash/epssmash). Two Docker wrapper scripts can be found at `epssmash/docker/`: `run_epssmash` and `download_epssmash_databases`. Moving these wrappers to your bin folder and running them will allow you to run epsSMASH without going through the installation above. Note that you must specify the input file and output directory when using `run_epssmash`:
```
run_antismash <input file> <output directory> [epsSMASH options]
```
If `download_epssmash_dependencies` is run before `run_epssmash` it will create the epsSMASH Docker image ebfore attempting to download the databases.

