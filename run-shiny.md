## Running Locally in R

1. Open R inside the repository path and install/load devtools to install the package dependencies:
```R
devtools::install(dependencies = TRUE)
```

2. Load the package and launch the app:
```R
library(seqpac)
run_seqpac_app()
```

## seqpac dashboard running guidelines (on demodata)

### Step 1: Landing Page & Loading Data
1. Open your browser to `http://localhost:3838`. You will see the **Home** (Launch Page) with the workflow overview.
2. Click the **Get Started** button. This will automatically redirect you to the **Load / Create PAC** tab.
3. On the **Example Data** sub-tab (left panel), click the **Load Drosophila Dataset** button.
4. You should see the **PAC Object Overview** populate with:
   - **Sequences (Rows)**: 9,131
   - **Samples (Columns)**: 9
   - **Annotations columns**: `Biotypes_mis0`, `Biotypes_mis3`
   - **Phenotypic metadata**: `stage`, `batch`, `sample`
5. Click through the three main-panel tabs (**Pheno Table**, **Counts Preview**, **Anno Preview**) to inspect the loaded datasets.

### Step 2: Filtering & Normalization
1. Click the **Filter & Normalize** tab in the main top navbar.
2. Keep the default settings:
   - **Size Range**: 20 to 30
   - **Min Counts (Threshold)**: 5
   - **Min Coverage (% of samples)**: 20%
   - **Normalizations**: Check **CPM** and **VST**
3. Click the **Apply Filter & Normalize** button.
4. Verify the outputs in the main panel:
   - **Sequence Counts Summary**: It will show that the filters retained **472 of 9131 sequences** (~5.17%).
   - **Available Normalizations**: Shows `cpm` and `vst` are now available.
   - **Sequence Length Distribution**: The histogram will refresh showing the distribution of the filtered sequences (mostly peak at 21–23 nt).

### Step 3: Annotation Explorer
1. Click the **Annotation Explorer** tab.
2. Under **Select Annotation Column to Analyze**, choose `Biotypes_mis0` (or `Biotypes_mis3`).
3. The **Annotation Statistics** table on the left will show the breakdown (e.g., how many sequences are annotated as `miRNA`, `tRNA`, `rRNA`, etc.).
4. The main panel table displays the complete `Anno` metadata matrix mapping sequences to their biotypes. You can search or filter this table interactively.

### Step 4: Post-Filtering Analysis
Click the **Post-Filtering Analysis** tab. This tab has four sub-tabs:

#### A. PCA (Principal Component Analysis)
1. Go to the **Principal Component Analysis** sub-tab.
2. Set **Color by Group** to `stage`.
3. Check the **Show Sample Labels** box.
4. Click **Run PCA**.
5. You will see a 2D scatter plot showing how your samples (e.g. Stage1, Stage3, Stage5) group together.

#### B. DESeq2 (Differential Expression)
1. Go to the **Differential Expression (DESeq2)** sub-tab.
2. In the design formula box, type `~ stage` (or keep the default).
3. Set **Primary Factor** to `stage`.
4. Click **Run DESeq2 Analysis**.
5. After a few seconds, the main panel will populate with a searchable interactive table of differential expression results (`baseMean`, `log2FoldChange`, `pvalue`, `padj`, etc.).

#### C. Size & Nucleotide Bias
1. Go to the **Size & Nucleotide Bias** sub-tab.
2. Set **Annotation Column** to `Biotypes_mis0`.
3. Choose **Size Distribution** and click **Generate Plots**. You will see biotype-specific sequence size histograms.
4. Switch to **Nucleotide Bias**, set position to `1`, and click **Generate Plots** to inspect the sequence starting nucleotide bias (e.g., high Uracil/T bias at position 1).

#### D. Composition (Bar & Pie)
1. Go to the **Composition (Bar & Pie)** sub-tab.
2. Set **Annotation Column** to `Biotypes_mis0`.
3. Select **Stacked Bar** and click **Generate Composition Chart** to see the relative abundance (percentages) of different small RNA biotypes across your samples.