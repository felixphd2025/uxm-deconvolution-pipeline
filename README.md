# UXM Deconvolution Pipeline

Nextflow pipeline for cell type deconvolution of Oxford Nanopore methylation data using UXM and wgbs_tools.

## Overview

This pipeline processes Oxford Nanopore BAM files to determine cell type composition through methylation-based deconvolution.

### Pipeline Steps

1. **MERGE_BAMS** - Merge multiple BAM files per sample
2. **BAM_TO_PAT** - Convert BAM to PAT format using wgbs_tools
3. **UXM_DECONV_U25** - Deconvolve with U25 atlas (25 tissue types)
4. **UXM_DECONV_U250** - Deconvolve with U250 atlas (250 tissue types)
5. **SUMMARIZE_RESULTS** - Combine and visualize results

## Requirements

### Software
- Nextflow (v25.04+)
- samtools
- Python 3.8+
- pandas

### Tools (must be installed locally)
- [wgbs_tools](https://github.com/nloyfer/wgbs_tools)
- [UXM_deconv](https://github.com/nloyfer/UXM_deconv)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/YOUR_USERNAME/uxm-deconvolution-pipeline.git
cd uxm-deconvolution-pipeline
```

2. Install dependencies:
```bash
# Install wgbs_tools
git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools
python setup.py
wgbstools init_genome hg38
wgbstools set_default_ref --name hg38

# Install UXM
git clone https://github.com/nloyfer/UXM_deconv.git
```

## Usage

### 1. Prepare sample sheet

Create `samples.tsv` with your samples:
```tsv
sample_id	bam_directory
sample1	/path/to/sample1/bam_pass
sample2	/path/to/sample2/bam_pass
```

### 2. Configure paths

Edit `nextflow.config` to set paths to your tools:
```groovy
params {
    atlas_u25 = "/path/to/UXM_deconv/supplemental/Atlas.U25.l4.hg38.tsv"
    atlas_u250 = "/path/to/UXM_deconv/supplemental/Atlas.U250.l4.hg38.full.tsv"
}
```

### 3. Run pipeline
```bash
./run_uxm.sh
```

Or with Nextflow directly:
```bash
nextflow run main.nf
```

### 4. Monitor progress
```bash
# Watch logs
tail -f .nextflow.log

# Check results
ls -lh results/
```

## Output Structure
```
results/
├── sample1/
│   ├── merged/
│   │   ├── sample1_merged.bam
│   │   └── sample1_merged.bam.bai
│   ├── pat/
│   │   └── sample1_merged.pat.gz
│   └── deconvolution/
│       ├── sample1_deconv_U25.csv
│       └── sample1_deconv_U250.csv
├── sample2/
│   └── ...
└── summary/
    ├── all_samples_U25.csv
    ├── all_samples_U250.csv
    └── summary_report.txt
```

## Configuration

### Hardware requirements

- **CPU**: 32+ cores recommended
- **RAM**: 128+ GB for large BAM files
- **Storage**: ~200GB per sample for intermediate files

### Performance

Processing time per sample (with 2,000-7,000 split BAMs):
- Merge: 1-5 hours
- BAM to PAT: 30 min - 1 hour
- Deconvolution: 10-30 minutes
- **Total: 2-6 hours per sample**

## Results Interpretation

### U25 Atlas (25 tissue types)
- Faster, cleaner visualization
- Good for identifying major cell populations
- Recommended for initial analysis

### U250 Atlas (250 tissue types)
- More detailed cell type resolution
- Useful for in-depth characterization
- May show more cross-reactivity

### Key cell types for brain tumors
- **Neuron** - Neuronal cells
- **Oligodend** - Oligodendrocytes (glial cells)
- **Endothel** - Endothelial cells (vasculature)
- **Blood-Mono+Macro** - Tumor-associated macrophages
- **Blood-NK/T/B** - Immune infiltration

## Citation

If you use this pipeline, please cite:

- UXM: [Loyfer et al., Nature Biotechnology, 2023](https://doi.org/10.1038/s41587-022-01465-8)
- wgbs_tools: [Loyfer et al., Bioinformatics, 2020](https://doi.org/10.1093/bioinformatics/btaa859)

## License

MIT License

## Author

Felix Adams  
University of Manchester  
Contact: felix.adams@postgrad.manchester.ac.uk

## Acknowledgments

- Developed for pediatric low-grade glioma (pLGG) research
- Based on Oxford Nanopore methylation sequencing data
