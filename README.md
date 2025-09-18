# snekwrap

directory structure based on CCDS template (https://cookiecutter-data-science.drivendata.org/)

written by Jackson C. Halpin


## overview
A collection of python wrappers to run various bioinformatics/structural biology tools in python


for example,
```python
import snekwrap
snekwrap.colabfold_batch_wrapper(
    input_file_or_directory="input_msa.a3m",
    output_directory="output_colabfold",
    weights="alphafold2_ptm",
    pairmode="unpaired",
    extra_args="--use_amber --num_relax 3",
)
```

The `snekwrap` library is intended to allow you to handle everything in python and mix and match different tools in a single script.

## implemented features
| Feature | Status |
|---------|--------|
| ColabFold | ✅ Implemented |
| USalign | 🚧 in progress |
| TMalign | 🚧 Planned |
| MAFFT | ✅ Implemented |
| MUSCLE | ✅ Implemented |
| Clustal Omega | ✅ Implemented |
| CD-HIT | ✅ Implemented |
| ESMFold | 🚧 in progress |
| proteinMPNN | 🚧 in progress |
| AlphaFold3 | 🚧 Planned |
| database query/download tools | ✅ Implemented |

<!-- | basic sequence manipulation | ✅ Implemented | -->
<!-- | basic structure manipulation | 🚧 planned | -->


## installation

```
git clone https://github.com/jacksonh1/snekwrap.git
cd snekwrap
```

run `make` to see available commands.

to install the project, run:
- `make create_environment` to create a conda environment with the required dependencies.
    - make sure to then activate the environment with `conda activate snekwrap`

if you can't use `make` (e.g. on Windows), you can manually create the conda environment with:
```
conda create -n snekwrap "python~=3.10"
conda env update --name snekwrap --file environment.yml
conda activate snekwrap
```


### configuration
you will need to configure paths to various executables used by the wrappers. there are several options for how to do this:


#### option 1: executables.yaml file
update the `./snekwrap/executables.yaml` file with paths to executable commands.
for example, here's what my `executables.yaml` file looks like:
```yaml
colabfold_batch: "/home/jch/tools/localcolabfold/colabfold-conda/bin/colabfold_batch"
colabfold_data: "/home/jch/tools/localcolabfold/colabfold"
chimerax: "/Applications/ChimeraX-1.9.app/Contents/bin/ChimeraX"
mafft: "mafft"
cd_hit: "cd-hit"
USalign: "USalign"
muscle: "/Users/jackson/tools/muscle/muscle-5.1.0/src/Darwin/muscle"
clustalo: "clustalo"
```

#### option 2: environment variables
alternatively, you can set environment variables for each executable.
for example, in bash:
```bash
export COLABFOLD_BATCH="/home/jch/tools/localcolabfold/colabfold-conda/bin/colabfold_batch"
export COLABFOLD_DATA="/home/jch/tools/localcolabfold/colabfold"
export CHIMERAX="/Applications/ChimeraX-1.9.app/Contents/bin/ChimeraX"
export MAFFT="mafft"
export CD_HIT="cd-hit"
export USALIGN="USalign"
export MUSCLE="/Users/jackson/tools/muscle/muscle-5.1.0/src/Darwin/muscle"
export CLUSTALO="clustalo"
```

#### option 3: .env file
you can also create a `.env` file in the root directory of the project with the following format:
```
COLABFOLD_BATCH=/home/jch/tools/localcolabfold/colabfold-conda/bin/colabfold_batch
COLABFOLD_DATA=/home/jch/tools/localcolabfold/colabfold
CHIMERAX=/Applications/ChimeraX-1.9.app/Contents/bin/ChimeraX
MAFFT=mafft
CD_HIT=cd-hit
USALIGN=USalign
MUSCLE=/Users/jackson/tools/muscle/muscle-5.1.0/src/Darwin/muscle
CLUSTALO=clustalo
```

#### option 4: default paths
if you do not set any of the above, it will assume the executables are in your PATH and will use the default names:
```yaml
colabfold_batch: "colabfold_batch"
colabfold_data: "/path/to/colabfold_data"
chimerax: "chimerax"
mafft: "mafft"
cd_hit: "cd-hit"
USalign: "USalign"
muscle: "muscle"
clustalo: "clustalo"
```
note that you will need to manually set the `colabfold_data` path in your environment if you are using option 4.



#### multiple configuration options
the configuration options can be mixed and matched. The order of precedence is:
1. .env file
2. system level environment variables (e.g. bash export)
3. executables.yaml file
4. default paths



