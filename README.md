# Unhearing the expected: subtractive prediction error emerges in the auditory thalamus

This repository contains the experimental paradigm and analysis code for an fMRI study investigating how the brain detects the absence of expected sensory stimuli.

> **Note:** This paper is currently under review. Citation information will be added upon publication.

## Overview

How does the brain detect what is absent? This study proposes that the detection of absence relies on a distinct form of prediction error dedicated to reducing beliefs about stimulus occurrence, termed **subtractive prediction error**. Using high-resolution fMRI of the human auditory system, we show that:

- Subtractive prediction error is encoded in the auditory thalamus (medial geniculate body; MGB) and primary auditory cortex (aHG)
- Unlike feature prediction error, subtractive prediction error is not encoded in the auditory midbrain (inferior colliculus; IC)
- These findings suggest that absence-related error signals emerge at the thalamic level and are supported by a distinct neural circuit

### Experimental Paradigm

The study used an auditory omission paradigm where participants listened to sequences of pure tones with occasional omissions. Participants (N=15) attended to sequences of up to 8 repetitions of pure tones (500ms ISI) and reported the position of omissions when present. Omissions occurred at positions 4, 5, or 6 with probability 6/7. Data were acquired using high-resolution partial-coverage fMRI optimized for subcortical auditory structures.

## Repository Structure

```
.
├── presentation/          # Experimental paradigm (MATLAB/Psychtoolbox)
│   └── omres.m           # Main experiment script
├── results/              # Preprocessed data for reproducing figures
│   ├── betas_t1w.pickle      # Beta estimates from first-level GLM
│   ├── funcloc_t1w.pickle    # Functional localizer contrasts
│   ├── noise_t1w.pickle      # Residual variance estimates
│   └── zscores_t1w.pickle    # Z-scored beta estimates
├── analysis.py           # Second-level analyses and figure generation
├── nitabs3.py            # Generic fMRI preprocessing and analysis library
├── omresconfig.py        # Study-specific configuration
└── pipeline.py           # First-level analysis pipeline
```

## Requirements

### For running the experiment (presentation/)
- MATLAB
- Psychophysics Toolbox

### For data analysis
- Python 3.7+
- Required Python packages:
  - numpy
  - scipy
  - pandas
  - matplotlib
  - seaborn
  - nibabel
  - nipype

### For preprocessing (first-level analysis)
- fMRIPrep 23.2.2+ (for preprocessing raw fMRI data)
- SPM12 (for GLM estimation)
- FreeSurfer 7.3.2+ (for cortical surface reconstruction)

## Usage

### Reproducing Figures and Tables

If you only want to reproduce the figures and statistical tables from the preprocessed results:

```bash
python analysis.py
```

This will generate all figures and supplementary tables in the `figures/` directory:
- `fig2-omres.pdf` - Adaptation and omission responses across regions
- `fig3-rois.pdf` - Negative prediction error effects
- `fig4-emergence.pdf` - Hierarchical emergence of prediction errors
- `figS2-allregs.pdf` - All regressors across ROIs
- `figS3-participants.pdf` - Individual participant results
- `figS4-residuals.pdf` - Residuals-weighted analyses
- `tabS1-omres.tex` - Statistics for adaptation and omission effects
- `tabS2-negpe.tex` - Statistics for negative prediction error
- `tabS3-emergence.tex` - Statistics for hierarchical emergence
- `tabS4-residuals.tex` - Statistics for residuals-weighted analyses

### First-Level Analysis Pipeline

To run the complete first-level analysis pipeline on raw data (requires raw fMRI data):

1. Configure paths in `omresconfig.py` to point to your data directories
2. Run the pipeline:

```bash
python pipeline.py
```

The pipeline performs:
1. **Preprocessing** (via fMRIPrep):
   - Motion correction
   - Coregistration to T1w
   - Normalization to MNI space
   - Confound extraction (motion, CSF, white matter)

2. **ROI Definition**:
   - Subcortical ROIs (IC, MGB): Anatomical atlas + functional localizer
   - Cortical ROIs (aHG): FreeSurfer parcellation

3. **GLM Estimation** (via SPM12):
   - Main experiment: 6 regressors (std0, std1, std2, om4, om5, om6)
   - Functional localizer: sound vs. silence contrast
   - Parametric modulation for repeated standards (std1, std2)

4. **Extraction**:
   - Beta estimates and z-scores for each condition and ROI
   - Functional localizer contrasts
   - Residual variance estimates

## Data Availability

### Preprocessed Results (Local)

The `results/` folder contains preprocessed derivatives necessary to reproduce all figures and statistics:
- Beta estimates from first-level GLM analyses
- Functional localizer contrasts
- Z-scored estimates
- Residual variance estimates

### First-Level Results (Zenodo)

First-level GLM results (beta images, contrast images, ROIs) for all subjects are publicly available on Zenodo:

**DOI:** [10.5281/zenodo.17603000](https://doi.org/10.5281/zenodo.17603000)
**URL:** https://zenodo.org/records/17603000

The Zenodo repository contains the complete derivatives folder for all 15 subjects, including:
- Beta images for all conditions in T1w and MNI space
- Contrast and t-statistic images
- Residual images
- ROI masks (anatomical and functional)
- JSON metadata files

#### Using the Zenodo Data

The Zenodo data contains all first-level GLM results needed to run the extraction and second-level analyses (lines 39-42 of `pipeline.py`).

After downloading from Zenodo, update the file paths to match your local system:

```python
import nitabs3

# Relocate koffer to your local path
nitabs3.relocate_koffer(
    '/path/in/downloaded/json',  # Path embedded in downloaded JSON files
    '/your/local/path/koffer',   # Where you actually saved the data
    mode='verify'                # Update paths without moving files
)
```

Then run the extraction steps in `pipeline.py` (lines 39-42) to reproduce the `results/` files, followed by `analysis.py` to generate all figures and tables.

### Raw Data (Unavailable)

Raw fMRI data (structural and functional images) are not publicly available because they contain personally identifiable biometric information and we do not have permission from participants to share them.

## Methods Summary

- **Participants:** 15 healthy adults (10 female, age 19-34)
- **Acquisition:** 3T Siemens Prisma, partial-coverage EPI (24 slices, 1.5mm isotropic, TR=2100ms)
- **Preprocessing:** fMRIPrep 23.2.2
- **GLM Estimation:** SPM12, including motion, CSF, and white matter nuisance regressors
- **ROI Definition:**
  - IC & MGB: Anatomical atlas (Sitek et al. 2019) + functional localizer
  - aHG: FreeSurfer Destrieux atlas
- **Statistics:** Non-parametric ranksum tests, Holm-Bonferroni correction for multiple comparisons

## Key Contrasts

1. **Adaptation:** std0 > std1 (response reduction to repeated standards)
2. **Omission response:** om4 > std1 (response to omission vs. expected sound)
3. **Subtractive prediction error:** om4 > om6 (position 4 omission > position 6 omission)

## Project Structure Details

### nitabs3.py
Generic library for fMRI data analysis providing:
- Interface to fMRIPrep preprocessing
- GLM estimation via SPM12/Nipype
- ROI extraction utilities
- Beta estimate and contrast extraction

### omresconfig.py
Study-specific configuration including:
- Data paths and directory structure
- ROI definitions (anatomical and functional)
- GLM design specifications for main task and localizer
- Event file parsing functions

### pipeline.py
Orchestrates the complete first-level analysis:
- Calls fMRIPrep for preprocessing
- Creates ROIs using anatomical and functional criteria
- Estimates GLMs for main task and localizer
- Extracts beta estimates, z-scores, and residual variance

### analysis.py
Implements all second-level (group) analyses:
- Loads preprocessed results from `results/`
- Statistical tests across participants and ROIs
- Generates all manuscript figures and tables
- Applies multiple comparison corrections

## Citation

This work is currently under review. Citation information will be provided upon publication.

**Preregistration:** The study was preregistered at OSF (https://osf.io/86a4x).

## Contact

For questions about the code or methods, please contact Alejandro Tabas.

## License

Code is provided for reproducibility and research purposes. Please cite the associated paper when using this code.
