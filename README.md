<h1>  Deconvolution benchmarking </h1>

Source code (Python v3.76.1 and R v4.2.0) to reproduce results described in the article **Performance of tumour microenvironment deconvolution methods in breast cancer using single-cell simulated bulk mixtures** submitted to Nature Communications on 31st July 2022.

- [System requirements](#system-requirements)
  - [Python dependencies](#python-dependencies)
  - [R dependencies](#r-dependencies)
    - [R/3.5.0 dependencies](#r350-dependencies)
    - [R/4.0.2 dependencies](#r402-dependencies)
  - [Singularity dependencies](#singularity-dependencies)
- [Experiment structure](#experiment-structure)
- [Prepare training and test data](#prepare-training-and-test-data)
- [Method execution instructions](#method-execution-instructions)
  - [BisqueRNA (bisque)](#bisquerna-bisque)
  - [BayesPrism (bprism)](#bayesprism-bprism)
  - [CIBERSORTx (CBX)](#cibersortx-cbx)
  - [Cell Population Mapping (CPM)](#cell-population-mapping-cpm)
  - [DWLS](#dwls)
  - [EPIC](#epic)
  - [hspe](#hspe)
  - [MuSiC](#music)
  - [Scaden](#scaden)
  - [Collect test results](#collect-test-results)

---

# System requirements
We generated results for this study using Python/3.6.1, R/3.5.0 on, and R/4.0.2 on 
CenOS Linux 7 (Core). Below are the dependencies for the Python and R languages. 

## Python dependencies
- cuda 10.1
- scaden 1.1.2
- pandas 1.1.5
- numpy 1.19.5
- anndata 0.6.22.post1
- scanpy 1.4.4.post1
- gtfparse 1.2.0
- plotly 4.1.0
- smote-variants 0.4.0

## R dependencies
The method DWLS requires R/3.5.0 environment while all other R-based methods require R/4.0.2 environment. 

### R/3.5.0 dependencies
- quadprog 1.5.5
- reshape 0.8.8
- e1071 1.6.8
- Seurat 2.3.3
- ROCR 1.0.7
- varhandle 2.0.5
- MAST 1.6.0

### R/4.0.2 dependencies
- BisqueRNA 1.0.5
- Biobase 2.50.0
- TED 1.4
- scBio 1.4
- Seurat 4.0.3
- foreach 1.5.1
- doSNOW 1.0.19
- parallel 4.0.2
- raster 3.4.13
- fields 12.5
- limma 3.46.0
- LiblineaR 2.10.12
- sp 1.4.5
- EPIC 1.1.5
- jsonlite 1.7.2
- hspe 0.1
- MuSiC 0.2.0
- xbioc 1.7.2


## Singularity dependencies
The method CIBERSORTx does not provide source code and can only be executed via a Docker container. As Docker requires root access, which could not be run on our infrastructure, we converted the CIBERSORT Docker images to a Singularity image. This requires:
- Docker 20.10.17
- Singularity 3.6.1


# Experiment structure
There are 5 experiments with the overall directory structure as follows:

```
.
├── 01_purity_levels_experiment/
│   ├── include_normal_epithelial/
│   │   └── data/
│   └── exclude_normal_epithelial/
│       └── data/
├── 02_normal_epithelial_lineages_experiment/
│       └── data/
├── 03_immune_lineages_experiment/
│   ├── minor_level/
│       └── data/
│   └── subset_level/
│       └── data/
├── params/
├── src/
├── .gitignore
└── README.md
```

Each experiment is compacted within its own directory with its own data and execution files to run each deconvolution method either as a Jupyter notebook or a `.psb` job.  <br>

Each experiment needs a `data/` directory, which you can download from matching folders on [Google Drive](). each `data/` folder should contain:

```
data/
├── Whole_miniatlas_meta_*.csv
├── Whole_miniatlas_umap.coords.tsv
├── training/ 
    ├── scRNA_ref.h5ad (single-cell reference profiles)
    └── training_sim_mixts.h5ad (training simulated RNA-Seq data, only for Scaden)
├── test/ 
    └── test_sim_mixts.h5ad (test simulated bulk RNA-seq data)
├── pretrained_scaden/ (pretrained Scaden models. Please refer to "Scaden" section)
└── expect_results/ (expected results. Please refer to "Collect test results" section)
    
```
After download, please place the `data/` folders in their corresponding experiments.

---

# Prepare training and test data
Each deconvolution requires either single-cell reference profiles or simulated mixtures as training data to estimate cell populations in the `test_sim_mixts.h5ad` files. To prepare this data, simply run the `0_prepare_method_specific_data.ipynb` notebook in each experiment directory. Once finished, method-specific sub-folders will be created in the `data/` directory:
```
data/
├── bisque/
├── bprism/
├── cbx/
├── cpm/
├── dwls/
├── epic/
├── hspe/
├── music/
├── scaden/
....
```
---
**NOTE:** You will need at least 4 CPUs and 275GiB memory to run this step for each experiment. <br>

# Method execution instructions  
Besides `CIBERSORTx` and `Scaden`, `.psb` run scripts of all methods point to the `R` code in `src/`. `CIBERSORTx` execution is based on the `Singularity`, and `Scaden` execution is based on the bash command `scaden` (please refere to [System requirements](#system-requirements)).

CPUs and memory requirements for each job are already specified in corresponding `.pbs` files.

## [BisqueRNA (bisque)](https://doi.org/10.1038/s41467-020-15816-6)

**Data files:**

- `logged_scRNA_ref.csv`: bisque requires the single-cell reference data to be logged
- `phenotypes.csv`: cell labels. Matches the order in which cells are listed in `logged_scRNA_ref.csv`

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_bisque.pbs` file
- Run method as pbs job

## [BayesPrism (bprism)](https://doi.org/10.1038/s43018-022-00356-3)
**Data files:**

- `scRNA_ref.csv`: single-cell reference profiles. Columns are cell ids and rows are gene symbols
- `single_cell_labels.csv`: single cell labels. Matches the order in which cells are listed in `scRNA_ref.csv`

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_bprism.pbs` file
- (Optional) Increase number of available CPUs in the `PBS -l` option and `N_CPUS` variable.
- Run method as pbs job

## [CIBERSORTx (CBX)](https://doi.org/10.1038/s41587-019-0114-2)
**Data files:**

- `scRNA_ref.txt`: the single-cell reference matrix. Columns are cells and column names are cell types. CBX infers cell types using column names and therefore does not need a single-cell labels files like other methods.
- `test_counts_*_pur_lvl.txt`: as I mentioned above, CBX requires the bulk expressions to be in its directory.<br>

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` in `1_run_cbx.pbs` file
- Copy `test_counts_<PUR_LVL>_pur_lvl.txt` files from `test/` to `cbx/`
- Download the `fractions_latest.sif` file from Google Drive. This is the same container provided at by [CIBERSORTx](https://hub.docker.com/r/cibersortx/fractions). We simply converted the image from Docker to Singularity, as Docker requires sudo access and Singularity does not
- Update Singularity commands
  - Path to the Singularity image `fractions_latest.sif`
  - `--token` with the token requested from [cibersortx.stanford.edu](https://cibersortx.stanford.edu/)

## [Cell Population Mapping (CPM)](https://doi.org/10.1038/s41592-019-0355-5)
**Data:**

We slightly modified the CPM source code to split the execution in half. The first half uses 48 CPUs and 300-500Gi of memory. The second half only uses 4 CPUs and 100Gi of memory. <br><br>

- `scRNA_ref_1330_per_ctype_*.txt`: single-cell reference data. We only use 1,330 or less per cell type.
- `cell_state.csv`: the cell-state space that CPM uses to map how cells transition from one cell type to another. We use UMAP coordinates as cell-state space.
- `single_cell_label.csv`: cell labels. Matches the order in which cells are listed in `scRNA_ref_1330_per_ctype_*.txt`

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_cpm_1.pbs` and `1_run_cpm_2.pbs` files
- Update path to the `1_CPM_functions.R` file in `src/1_run_cpm_1.R` and `src/1_run_cpm_2.R`
- Run `1_run_cpm_1.pbs` first as a pbs job
- Wait for `1_run_cpm_1.pbs` to finish. Then run `src/1_run_cpm_2.pbs` as PBS job

## [DWLS](https://doi.org/10.1038/s41467-019-10802-z)
**Data:**

- `scRNA_ref.txt`: the single-cell reference matrix
- `single_cell_labels.csv`: cell labels. Matches the order in which cells are listed in single-cell reference

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_dwls.pbs` file
- Run method as pbs job

## [EPIC](https://doi.org/10.7554/eLife.26476)
**Data:**

EPIC is the only method that requires a pre-built signature matrix instead of single-cell reference. We provided EPIC with the signature matrix that CBX generated during its execution. This is the reason why you see data files for EPIC in the `cbx_sig_matrix` folder.

- `reference_profiles.csv`: signature matrix that CBX produced
- `marker_gene_symbols.csv`: marker genes among all genes in the signature matrix. EPIC weights these genes higher non non-marker genes. In our study, this is merely a procedural necessity as we list all genes in the signature matrix as marker genes.

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_epic.pbs` file
- Run method as pbs job

## [hspe](https://doi.org/10.1214/20-AOAS1395)
**Data:**

- `scRNA_ref.csv`: the single-cell reference matrix
- `pure_samples.json`: cell labels in a dictionary. Matches the order in which cells are listed in single-cell reference
- `logged_test_counts_*_pur_lvl.txt`: as mentioned above, hspe requires the bulk expressions to be logged. 

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_hspe.pbs` file
- Run method as pbs job

## [MuSiC](https://doi.org/10.1038/s41467-018-08023-x)
**Data:**

- `scRNA_ref.csv`: the single-cell reference matrix
- `pheno.csv`: cell labels and patient ids of cells in the single-cell reference. Matches the order in which cells are listed in single-cell reference

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_music.pbs` file
- Run method as pbs job

## [Scaden](https://doi.org/10.1126/sciadv.aba2619)
**Data:**

- `train_counts.h5ad`: Scaden is the only methods that requires simulated bulk mixtures as training data. It also requires these mixtures to be compressed into an AnnData (.h5ad) object.
- `test_counts_*_pur_lvl.txt`: as I mentioned above, CBX requires the bulk expressions to be in its directory.

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` in `1_run_scaden.pbs` file
- Copy `test_counts_<PUR_LVL>_pur_lvl.txt` files from `test/` folder into `scaden/` folder

**NOTE 1:** <br>
scaden is built using the `tensorflow` library and requires `cuda 10.1` to run. In this study, we installed `scaden` and `cuda 10.1` in 2 separate environments. This is why you see this in `scaden.pbs`
```
source /software/scaden/scaden-0.9.4-venv/bin/activate
module load gpu/cuda/10.1
```
This is not compulsory as long as you can ensure both `scaden` and `cuda 10.1` are loaded before running `scaden`.

**NOTE 2:** <br>
The `Scaden` method is a neural network and therefore produces slightly different results when trained from scratch. However, results from a newly trained model should be similar to results generated by the pre-trained models. 

For testing Scaden pre-trained models:
- Edit the `01_run_scaden.pbs` script to point to `data/pretrained_scaden/` instead of `data/scaden/`
Copy `test_counts_<PUR_LVL>_pur_lvl.txt` files from `test/` folder into `pretrained_scaden/` folder
- Comment out `scaden train` in `01_run_scaden.pbs`. `scaden process` should produce identical processed data files compared to `data/scaden`
- Run `01_run_scaden.pbs`

---

## Collect test results

Once all methods have been executed, you can collect results from all methods using the `2_collate_test_results.ipynb` Jupyter notebook in each experiment directory.

Simply update the `prefix` variable at the top of each notebook with the correct path to the experiment folder and run each notebook from top to bottom. After that, all methods' results will be collated into the `data/results` folder.

In addition, you can also find all results in this study in `data/expected_results/`. Each `expected_results/` sub-folder of each experiment should contain:

```
└── data/
    └── expect_results/
        ├── bisque.csv
        ├── bprism.csv
        ├── cbx.csv
        ├── cpm.csv
        ├── dwls.csv
        ├── epic.csv
        ├── hspe.csv
        ├── music.csv
        ├── scaden.csv
        └── truth.csv
```
