# Installation

This guide describes how to install and configure `fnirspower`.

The GitHub repository contains:

- the `fnirspower` MATLAB package for measurement prediction
- the `fnirspower_py` Python package for power analysis
- MATLAB and Python examples
- Slurm job-submission scripts
- documentation.

The required third-party MATLAB dependencies and workspace assets are distributed separately as multipart archive files attached to the corresponding GitHub Release.

---

Installation consists of:

1. Cloning or downloading the source-code repository and assets 
2. Installing MATLAB and the required toolboxes
3. Creating the Python environment
4. Verifying the MATLAB and Python installations

---

## 1. Download the repository and release assets

### A. Download the source code

Clone the repository using Git:

```bash
git clone https://github.com/esbulger/fnirs-power.git
cd <repository>
```

Alternatively, download the source-code archive from GitHub and extract it to a local directory.

### B. Download the required release assets

The third-party MATLAB dependencies and workspace files are not included when the repository is cloned.

Open the GitHub **Releases** page for the version you are installing 
and download every zip file corresponding to a repo folder.

Using WinRAR, 7-Zip, or another compatible archive manager, extract them into the root of the cloned repository.

After extraction, the repository should contain:

```text
fnirspower/
├─ docs/
│
├─ matlab/
│  ├─ +fnirspower/
│  ├─ examples/
│  └─ thirdparty/
│
├─ python/
│  ├─ fnirspower_py/
│  ├─ examples/
│  ├─ data/
│  ├─ results/
│  └─ jobs/
│
└─ workspace/
   ├─ derivatives/
   ├─ figures/
   ├─ layouts/
   ├─ models/
   ├─ rawdata/
   └─ montages/
```

Some workspace directories may be empty; they are present as a placeholder to encourage compatible structuring of files.

---

## 2. MATLAB requirements

### MATLAB version

`fnirspower` was developed in **MATLAB R2020b**. It may also function with newer MATLAB versions.

MATLAB is available from:

<https://www.mathworks.com/products/matlab.html>

### Required MATLAB toolboxes

The current code requires:

- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox

### Third-party MATLAB dependencies

The required third-party MATLAB packages are distributed with the GitHub Release and should be extracted into:

```text
matlab/thirdparty/
```

These include:

- `iso2mesh`
- `EasyH5`
- `jsnirfy`
- `NIRSToolbox`
- `FieldTrip`
- `NIRFASTer`

The supplied versions should remain in their expected locations because `fnirspower.setup_paths` uses this directory structure when configuring MATLAB.

### MATLAB path setup

Open MATLAB and define the location of the cloned repository:

```matlab
repo_root = 'C:\path\to\fnirspower';
```

Add the repository's `matlab` directory to the MATLAB path, then run the package setup function:

```matlab
addpath(fullfile(repo_root, 'matlab'));
P = fnirspower.setup_paths();
```

This should:

- add the parent directory of `+fnirspower`;
- add the supplied third-party MATLAB dependencies; and
- define paths to the workspace, model, layout, montage, and output directories.

To verify the setup, run:

```matlab
which fnirspower.pipeline.run_measurement_prediction
which fnirspower.measpred.compute_subject_ROI
```

MATLAB should return the corresponding files inside the cloned repository.

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

On a computing cluster, create a virtual environment from the repository root:

```bash
python -m venv fnirspower_py_env
```

On Linux, macOS, or a computing cluster:

```bash
source fnirspower_py_env/bin/activate
```

On Windows PowerShell:

```powershell
.\fnirspower_py_env\Scripts\Activate.ps1
```

Install the required packages:

```bash
python -m pip install --upgrade pip
python -m pip install -r python/requirements.txt
```

Verify the Python setup:

```bash
cd python
python -c "import fnirspower_py; import mne; print('Python setup succeeded')"
```

The location of the virtual environment may vary. Running the power-analysis scripts on a computing cluster typically requires cloning the repository, or copying its `python` directory, to a cluster-visible location and creating the Python environment there.

All supplied power-analysis scripts are designed to support computing-cluster execution.

---

## 4. Required project assets

The GitHub Release provides the model, layout, montage, and example workflow assets required by the supplied examples. These files should appear under `workspace/` after the release archive is extracted.

### A. Head model

The current codebase is designed to work with the head model included in the release assets.

This head model was developed using iso2mesh and 3D Slicer. Advanced users may substitute their own head model but must ensure compatibility with the downstream workflows.

### B. Forward-model files

Some workflows require a forward-model MAT file, such as the `*_nirsmodel.mat` files used by the measurement-prediction pipeline.

These files are used by workflows including:

- group GLM;
- measurement prediction; and
- variance estimation.

Forward-model files for the example head model and selected montages are included in the release assets. Custom head models or montages may require running the `run_forward_model` pipeline and exporting the result in the format expected by the downstream functions.

### C. Probe geometry and `probeInfo`

The package is designed around **NIRx-compatible probe geometry**.

Forward-model and related workflows require a valid `probeInfo` structure containing the expected source, detector, and channel definitions.

### D. Layout and montage files

Several plotting steps require:

- a layout MAT file;
- a montage name; or
- both.

The package is configured to use layout and montage files under the workspace structure.

Custom layout files can be generated using FieldTrip. See:

```text
matlab/+fnirspower/+helpers/prepare_layout_NIRx.m
```

### E. Raw subject data

Users adapting the preprocessing, GLM, absorption-estimation, or variance-estimation workflows 
can provide their own compatible data in the `workspace/rawdata` folder.

Raw subject data are not currently included in the public release. They are not required to run the 
supplied measurement-prediction and power-analysis example. Subject data may be made available upon 
request or a future release.

The current code supports explicit SNIRF path lists and is designed around NIRx-compatible data handling.

---

## 5. Cluster and Slurm setup

The Python utilities support cluster execution through Slurm batch scripts.

Users should have:

- access to a Slurm-based computing cluster;
- a working Python environment on the cluster;
- valid cluster-visible paths for inputs, outputs, and logs; and
- permission to submit batch jobs.

On the cluster, the relevant project structure should be available as:

```text
fnirspower/
└─ python/
   ├─ fnirspower_py/    # Python functions and utilities
   ├─ examples/         # Power-analysis example scripts
   ├─ data/             # Input data, including MATLAB output
   ├─ results/          # Power-analysis results
   │  └─ plots/         # Generated power plots
   └─ jobs/             # Slurm job-submission scripts
```

If a computing cluster is not available, the scripts can be adapted to run locally. This is not recommended for large simulations because of their computational cost.

Example Slurm scripts are provided under:

```text
python/jobs/
```

---

## 6. Next step

Once installation is complete, continue to the [getting started guide](getting_started.md).

That guide walks through one complete example from MATLAB setup through Python-based power analysis.

---

## Troubleshooting

### MATLAB cannot find `fnirspower`

Confirm that the repository's `matlab` directory was added before calling the setup function:

```matlab
repo_root = 'C:\path\to\fnirspower';
addpath(fullfile(repo_root, 'matlab'));
P = fnirspower.setup_paths();
```

The path entry must be the parent directory of `+fnirspower`, not the `+fnirspower` directory itself.

### FieldTrip causes path conflicts

FieldTrip should generally be added only at the top level rather than recursively adding all subfolders.

If another FieldTrip installation is already on the MATLAB path, remove it before running `fnirspower.setup_paths`, or modify the relevant path setup code to use the intended installation.

### Python cannot import `fnirspower_py`

Run the example scripts from the repository's `python` directory:

```bash
cd /path/to/fnirspower/python
python examples/sim_cluster_based_power_parallelized.py
```

Confirm that the following directories are present:

```text
python/
├─ fnirspower_py/
└─ examples/
```

You can also test the import directly:

```bash
python -c "import fnirspower_py; print('fnirspower_py import succeeded')"
```

### Saving MATLAB files is slow

Saving large arrays with `-v7.3` can be slow.

Where practical, save only the output variables required by later steps rather than the complete intermediate structures.
