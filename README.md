<h1>  Deconvolution benchmarking </h1>

Source code (Python v3.76.1 and R v4.2.0) to reproduce results descrbied in the article **Performance of tumour microenvironment deconvolution methods in breast cancer using single-cell simulated bulk mixtures** submitted to Nature Communications on 31st July 2022.

---

<h2> Overview </h2> 
There are 5 experiments with the overal directory structure as follows:

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

Each experiment is compacted within its own directory with data files and `.pbs` files to run each method on HPC. With the execption of  `CIBERSORTx` and `Scaden`, `.pbs` files of all other methods point to the corresponding `.R` files in the `src/` folder, which includes their execution code. <br>

Data files in the `data/` folder are organised into methods-specific sub-folders as follows: <br>
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
├── results/
└── test/

```

The `data/test` folder contains the bulk expressions that we deconvolute. `CBX` and `Scaden` require these files to be copied into their folders. `hspe` requires the bulk expressions to be logged and saved in its directory. All others methods require are pointed to the bulk expression files in `data/test` during their execution (see `.pbs` files).

We have compressed and uploaded the training and test data for each experiment to [Figshare](???). You will need to download and extract the `.tar` files in the matching experiment folder, e.g. `exclude_normal_epithelial.tar` should be extracted inside `01_purity_levels_experiment/exclude_normal_epithelial/data/`. This will give you all the required training & test data organised in the structure described above.

The `data/results/` folder of each experiment contains the results of all methods. To collect methods' results, refer to the section **Collect test results**.

---

<h2> Method execution instructions  </h2> 

<h2> bisque </h2>

**Data files:**

- `logged_scRNA_ref.csv`: bisque requires the single-cell reference data to be logged
- `phenotypes.csv`: cell labels. Matches the order in which cells are listed in `logged_scRNA_ref.csv`

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_bisque.pbs` file
- Run method as pbs job

<h2> BayesPrism (bprism) </h2>

**Data files:**

- `scRNA_ref.csv`: single-cell reference profiles. Columns are cell ids and rows are gene symbols
- `single_cell_labels.csv`: single cell labels. Matches the order in which cells are listed in `scRNA_ref.csv`

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_bprism.pbs` file
- (Optional) Increase number of available CPUs in the `PBS -l` option and `N_CPUS` variable.
- Run method as pbs job

<h2> CIBERSORTx (CBX) </h2>

**Data files:**

- `scRNA_ref.txt`: the single-cell reference matrix. Columns are cells and column names are cell types. CBX infers cell types using column names and therefore does not need a single-cell labels files like other methods.
- `test_counts_*_pur_lvl.txt`: as I mentioned above, CBX requires the bulk expressions to be in its directory.<br>

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` in `1_run_cbx.pbs` file
- Make sure `test_counts_<PUR_LVL>_pur_lvl.txt` files are copied from `test/` folder into `cbx/` folder
- Download CIBERSORTx Docker images from [Docker hub](https://hub.docker.com/r/cibersortx/fractions) and convert it to Singularity image
- Update Singularity commands
  - Path to Singularity image
  - `--token` with the token requested from [cibersortx.stanford.edu](https://cibersortx.stanford.edu/)

*NOTE*: CIBERSORTx requires a token to run, which can be requested at <i>cibersortx.stanford.edu</i>

<h2> CPM </h2>

**Data:**

We slightly modified the CPM source code to split the execution in half. The first half uses 48 CPUs and 300-500Gi of memory. The second half only uses 4 CPUs and 100Gi of memory. <br><br>

- `scRNA_ref_1330_per_ctype_*.txt`: single-cell reference data. We only use 1,330 or less per cell type.
- `cell_state.csv`: the cell-state space that CPM uses to map how cells transition from one cell type to another. We use UMAP coordinates as cell-state space.
- `single_cell_label.csv`: cell labels. Matches the order in which cells are listed in `scRNA_ref_1330_per_ctype_*.txt`

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_cpm_1.pbs` and `1_run_cpm_2.pbs` files
- Update path to the `1_CPM_functions.R` file in `1_run_cpm_1.pbs`
- Run `1_run_cpm_1.pbs` first as a pbs job
- Wait for `1_run_cpm_1.pbs` to finish. Then run `1_run_cpm_2.pbs` as PBS job

<h2> DWLS </h2>

**Data:**

- `scRNA_ref.txt`: the single-cell reference matrix
- `single_cell_labels.csv`: cell labels. Matches the order in which cells are listed in single-cell reference

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_dwls.pbs` file
- Run method as pbs job

<h2> EPIC </h2>

**Data:**

EPIC is the only method that requires a pre-built signature matrix instead of single-cell reference. We provided EPIC with the signature matrix that CBX generated during its execution. This is the reason why you see data files for EPIC in the `cbx_sig_matrix` folder.

- `reference_profiles.csv`: signature matrix that CBX produced
- `marker_gene_symbols.csv`: marker genes among all genes in the signature matrix. EPIC weights these genes higher non non-marker genes. In our study, this is merely a procedural necessity as we list all genes in the signature matrix as marker genes.

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_epic.pbs` file
- Run method as pbs job

<h2> hspe </h2>

**Data:**

- `scRNA_ref.csv`: the single-cell reference matrix
- `pure_samples.json`: cell labels in a dictionary. Matches the order in which cells are listed in single-cell reference
- `logged_test_counts_*_pur_lvl.txt`: as mentioned above, hspe requires the bulk expressions to be logged. 

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_hspe.pbs` file
- Run method as pbs job

<h2> MuSiC </h2>

**Data:**

- `scRNA_ref.csv`: the single-cell reference matrix
- `pheno.csv`: cell labels and patient ids of cells in the single-cell reference. Matches the order in which cells are listed in single-cell reference

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` and `SRC_DIR` paths in `1_run_music.pbs` file
- Run method as pbs job

<h2> Scaden </h2>

**Data:**

- `train_counts_kondrashova.h5ad`: Scaden is the only methods that requires simulated bulk mixtures as training data. It also requires these mixtures to be compressed into an AnnData (.h5ad) object.
- `test_counts_*_pur_lvl.txt`: as I mentioned above, CBX requires the bulk expressions to be in its directory.

**Running:**
- Replace user email in `PBS -M` option
- Modify `WORK_DIR` in `1_run_scaden.pbs` file
- Make sure `test_counts_<PUR_LVL>_pur_lvl.txt` files are copied from `test/` folder into `scaden/` folder

---

<h2> Collect test results </h2> 

Once all methods have been executed, you can collect results from all methods using the `2_collate_test_results.ipynb` Jupyter notebook in each experiment directory.

Simply update the `prefix` variable at the top of each notebook with the correct path to the experiment folder and run each notebook from top to bottom. After that, all methods' results will be collated into the `data/results` folder.
