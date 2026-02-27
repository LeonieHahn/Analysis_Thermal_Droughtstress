## Capturing early signs of drought stress in tree saplings to support irrigation monitoring at drought-prone reforestation sites in Central Europe with a low-cost thermal camera

This repository contains the code for reproducing the data analysis, results and figures from Hahn et al. (2026).

### Workflow and execution order
Please run the scripts **in the order indicated by their numbering**.

##### Python Scripts

The following notebooks must be executed in **Python**:

- `03_Image_Registration.ipynb`
- `05_Segment_Trees_Negative_Samples_SAM2.ipynb`

Use the provided conda environment:

```bash
conda env create -f Analysis_Thermal_Droughtstress.yml
conda activate Analysis_Thermal_Droughtstress
```

##### R scripts
All remaining scripts are written in **R**.

### Citation
This repository is archived on Zenodo
