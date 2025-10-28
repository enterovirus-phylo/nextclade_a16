# Nextclade Workflow for &lt;your virus&gt;

This repository provides a robust, reproducible workflow for building a custom [Nextclade](https://github.com/nextstrain/nextclade) dataset for &lt;your virus&gt;. It enables you to generate reference and annotation files, download and process sequence data, infer an ancestral root, and create all files needed for Nextclade analyses and visualization.

---

## Folder Structure

Follow the [Nextclade example workflow](https://github.com/nextstrain/nextclade_data/tree/master/docs/example-workflow) or use the structure below:

```bash
mkdir -p dataset data ingest resources results scripts
```

---

## Workflow Overview

This workflow includes several modular steps:

1. **Reference Generation**  
   Extracts reference and annotation files from GenBank.
2. **Dataset Ingest**  
   Downloads and processes sequences and metadata from NCBI Virus.
3. **Phylogenetic Root Inference (optional)**  
   Infers a dataset-specific ancestral root sequence to use as a reference in Nextclade, improving mutation and clade assignments.
4. **Augur Phylogenetics & Nextclade Preparation**  
   Builds phylogenetic trees, prepares sequence alignments, and generates all required files for Nextclade and Auspice.
5. **Visualization & Analysis**  
   Enables both command-line and web-based Nextclade analyses, including local dataset hosting.

---

## Setup Instructions

### 1. Generate Reference Files

Run the script to extract the reference FASTA and genome annotation from GenBank:

```bash
python3 scripts/generate_from_genbank.py --reference "<accession_id>" --output-dir dataset/
```

During execution, follow the prompts for CDS annotation selection.

**Outputs:**
- `dataset/reference.fasta`
- `dataset/genome_annotation.gff3`

---

### 2. Configure `pathogen.json`

Edit `pathogen.json` to:
- Reference your generated files (`reference.fasta`, `genome_annotation.gff3`)
- Update metadata and QC settings as needed

> [!WARNING]  
> If QC is not set, Nextclade will skip quality checks.

See the [Nextclade pathogen config documentation](https://docs.nextstrain.org/projects/nextclade/en/latest/user/input-files/05-pathogen-config.html) for details.

---

### 3. Prepare GenBank Reference

Copy your GenBank file to `resources/reference.gb`.  
Edit protein names and features if necessary for your use case.

---

### 4. Update the `Snakefile`

- Adjust the workflow parameters and file paths as needed for your dataset.
- Ensure required files are available:
  - `data/sequences.fasta`
  - `data/metadata.tsv`
  - `resources/auspice_config.json`

Sequences and metadata can be downloaded automatically via the ingest process (see below).

---

## Subprocesses

### Ingest

Automates downloading of EV sequences and metadata from NCBI Virus.  
See [ingest/README.md](ingest/README.md) for specifics.

**Required packages:**  
`csvtk, nextclade, tsv-utils, seqkit, zip, unzip, entrez-direct, ncbi-datasets-cli` (installable via conda-forge/bioconda)

---

### Inferred Root (Optional but Recommended)

The `inferred-root/` directory contains a reproducible pipeline to infer a dataset-specific ancestral root, which can be used as a reference sequence in Nextclade. This enhances mutation and clade call accuracy for your dataset.

- **See:** [`inferred-root/README.md`](inferred-root/README.md) for details.
- To enable, set `ANCESTRAL_ROOT_INFERRENCE = True` in your config and run with  
  `--config root_inference_confirmed=true`.
- Without confirmation, the workflow will halt and display an opt-in message.

> [!NOTE]  
> To skip the inferred root step, leave `ANCESTRAL_ROOT_INFERRENCE = False`.

---

### 5. Configure `auspice_config.json`

Edit `resources/auspice_config.json` to:
- Match the **RefSeq node** to your virus reference  
- Update the **title**, **build_url**, and **maintainers**  

This enables three reference options for mutation calling:  

- **Reference** → RefSeq or commonly used reference sequence
- **Static Ancestral Root** → Ancestral sequence from `inferred-root` workflow; run only once
- **Tree root** → Tree-based root (changes with every re-run)  

*Example reference selector UI:*

![Options](./image.png)

---

## Running the Workflow

To generate the Auspice JSON and a Nextclade example dataset:

```bash
snakemake --cores 9 all --config root_inference_confirmed=true
```

This will:
- Build the reference tree and produce the Nextclade dataset in `dataset/`
- Run Nextclade on the example sequences in `out-dataset/sequences.fasta`
- Output results to `test_out/` (alignment, translations, summary TSV)

### Labeling Mutations of Interest
To label mutations of interest, execute the `mutLabels` rule as a standalone instance. They will be added to the `out-dataset/pathogen.json` file.

---

## Visualizing Your Custom Nextclade Dataset

To use the dataset in Nextclade Web, serve it locally:

```bash
serve --cors out-dataset -l 3000
```

Then open:

```
https://master.clades.nextstrain.org/?dataset-url=http://localhost:3000
```

- Click "Load example", then "Run"
- Consider reducing "Max. nucleotide markers" to 500 under "Settings" → "Sequence view" to optimize performance

---

## Author & Contact

- Maintainers: Nadia Neuner-Jehle, Alejandra Gonzalez Sanchez, Emma B. Hodcroft ([hodcroftlab](https://github.com/hodcroftlab))
- For questions or suggestions, please [open an issue](https://github.com/enterovirus-phylo/dataset-template-inferred-root/issues/new) or email: eve-group[at]swisstph.ch

---

## Troubleshooting and Further Help

- For issues, see the [official Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/index.html#).
- For details on the inferred root workflow, see [`inferred-root/README.md`](inferred-root/README.md).

---

This guide provides a structured, scalable approach to building and using high-quality Nextclade datasets for your specific Enterovirus — and can be adapted for other enterovirus types as well.
