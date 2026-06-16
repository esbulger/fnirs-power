# Installation

This guide describes how to install and configure `fnirspower` based on the current repository structure and codebase.

The repository provides:
- the fnirspower MATLAB package enabling measurement prediction
- the fnirspower Python package enabling statistical testing
- bundled third-party MATLAB packages
- examples and documentation
- the expected workspace-style directory structure

---

In general, installation consists of:
1. downloading or cloning the repository
2. installing MATLAB and necessary toolboxes
3. creating the Python environment
4. verifying both the MATLAB and Python sides with a small test

---

## 1. Repository structure

The current repository is organized as:

```text
fnirspower/
├─ docs/                #   Installation, Getting Started, and Other Documentation
│
├─ matlab/
│  ├─ +fnirspower/      #   Core MATLAB package code
│  ├─ examples/         #   MATLAB example scripts and workflow drivers
│  └─ thirdparty/       #   Bundled third-party MATLAB dependencies
│
├─ python/
│  ├─ fnirspower_py/    #   Python utilities and helper modules
│  ├─ examples/         #   Python example scripts
│  ├─ data/             #   Data to load (i.e., MATLAB output)
│  ├─ results/          #   Power analysis results and plots
│  └─ jobs/             #   Slurm/job submission scripts
│
└─ workspace/           #   Models, layouts, montages, raw data, and derivative outputs
│  ├─ derivatives/
│  ├─ figures/
│  ├─ layouts/
│  ├─ manuscripts/
│  ├─ models/
│  ├─ montages/
│  ├─ rawdata/
│  └─ simulations/
```

The package is easiest to use when this structure is kept intact.

---

## 2. MATLAB requirements

### MATLAB version
`fnirspower` was developed in **MATLAB R2020b**. It may still function properly with newer MATLAB versions.
You can find MATLAB here: https://www.mathworks.com/products/matlab.html
### MATLAB toolboxes required
The code currently requires the following MATLAB toolboxes:

- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox 
- Optimization Toolbox
- Parallel Computing Toolbox

### Bundled third-party MATLAB dependencies
The current package structure expects bundled third-party folders under:

```text
matlab/thirdparty/
```

These include:
- `iso2mesh`: https://iso2mesh.sourceforge.net/cgi-bin/index.cgi
- `EasyH5`: https://github.com/NeuroJSON/easyh5
- `jsnirfy`: https://www.mathworks.com/matlabcentral/fileexchange/180495-fnirs-snirf-format-reader-writer-for-matlab-octave
- `NIRSToolbox`: https://github.com/huppertt/nirs-toolbox
- `FieldTrip`: https://www.fieldtriptoolbox.org/
- `NIRFASTer`: https://github.com/nirfaster/NIRFASTer

These should remain in the expected repository locations.

### MATLAB path setup
After opening MATLAB, add `fnirspower/matlab` to path and run:

```matlab
P = fnirspower.setup_paths();
```

This should:
- add the parent folder of `+fnirspower`
- add bundled third-party dependencies
- add all other paths according to the repository structure

---


## 3. Python requirements

The Python side of the project is used for:
- cluster-based permutation test power simulations
- Slurm-based batch execution
- result processing and plotting

### Python environment

The Python utilities were developed for **Python 3.9** and currently depend on:

- `numpy`
- `scipy`
- `matplotlib`
- `mne`
- `h5py`

The recommended setup is to install these with `pip` using the provided requirements file:

```bash
python -m venv fnirspower_py_env
source ./fnirspower_py_env/bin/activate
(fnirspower_py_env) [ ~]$ pip install -r python/requirements.txt
```

for other modules, you can use: 

```
(fnirspower_py_env) [ ~]$ pip install <module(s) you want to install>
```

The actual location of the environment may vary. We generally recommend using our repo structure hosted on a PC, 
however, running scripts on a computing cluster typically requires cloning the python folder onto the computing cluster,
and setting up the environment there. Individual implementations may vary. 
All provided code for power analyses is designed to run on a computing cluster.

---

## 4. Required project assets

The repository provides the code and expected structure, but several workflows also require project assets (i.e. data, models, etc.).

### A. Head model
The current codebase is designed to work with the provided head model in the repository/workspace structure. 
This head model was developed using iso2mesh and 3DSlicer. Advanced users may use their own head model, 
but must ensure downstream compatibility.

### B. Forward model files
Some workflows assume the presence of a forward-model MAT file, such as the `*_nirsmodel.mat` 
files produced by the forward-model pipeline.

These are used by workflows including:
- group GLM
- measurement prediction
- variance estimation

Forward model files are included for the example head model and some montages, however, 
custom head models and montages may require the `run_forward_model` pipeline.

### C. Probe geometry / `probeInfo`
The package is designed around **NIRx-compatible probe geometry**.

Forward-model and related workflows require a valid `probeInfo` structure. 
In practice, users should provide a compatible `probeInfo` MAT file with the expected source, detector, and channel definitions.

### D. Layout / montage files
Several plotting steps require:
- a layout MAT file
- a montage name
- or both

The package is currently set up to work with layout files in the project structure. 
Custom layout files are generated using FieldTrip. See `matlab/+fnirspower/+helpers/prepare_layout_NIRx.m` 

### E. Raw subject data

The current code supports explicit SNIRF path lists and is designed around NIRx-compatible data handling. 
This is not required for conducting measurement prediction and power analysis. 
It is required for customizing this pipeline to your own data / experimental setup.

---

## 5. Cluster / Slurm setup

The current Python utilities are designed to support cluster execution through Slurm batch scripts.

Users should also have:
- access to the cluster
- a working Python environment on the cluster
- valid cluster-visible paths for inputs, outputs, and logs
- permission to submit batch jobs

On the cluster, the directory should be a clone or copy of the python folder:

```text
fnirspower/
├─ python/
│  ├─ fnirspower_py/    #   Python functions and utilities
│  ├─ examples/         #   Example scripts for running power analysis
│  ├─ data/             #   Data structures and data to load (e.g., MATLAB output)
│  ├─ results/          #   Power analysis results
│  │ └─ plots/          #   Power plots
│  └─ jobs/             #   Slurm/job submission scripts
```

If access to a computing cluster is not available, the code can be modified to run locally. 
This is not recommended due to long computation times. 
If necessary, please contact the repo maintainer for assistance. 

We provide example SLURM submission scripts to submit job arrays.

---

## 6. Next step

Once installation is complete, proceed to our [getting started guide](docs/getting_started.md).

This guide walks through one full example from MATLAB setup through Python-based power analysis.


___

<br>



### Possible issues

See below for some issues that may occur during installation.

#### MATLAB cannot find `fnirspower`
Use:

```matlab
P = fnirspower.setup_paths();
```

and confirm that the parent folder of `+fnirspower` is on the MATLAB path.

#### FieldTrip causes path conflicts
FieldTrip should generally be added only at the top level rather than recursively adding all subfolders. If there is a pre-existing installation, it should either be removed from the path, or the relevant path files in this repo should be modified.


#### Simulation jobs fail with multiprocessing / pickle errors
If you encounter multiprocessing issues, start with:
- a small test configuration
- a small subject/block grid
- one CPU per task

Then scale up after the test succeeds.

#### Saving MATLAB files is slow
This can happen when large arrays are saved with `-v7.3`. 
In those cases, consider saving only the required outputs rather than entire output structures.

---