# 📈 AP Intensity Profile Extractor – Cortical Plane Confocal Images

This MATLAB script extracts and analyses anterior-posterior (AP) intensity profiles from cortical plane confocal images of the cortex or membrane.

The images processed by this script represent cortical bands that are 100 pixels high but may vary in length due to differences in cortex coverage across embryos.

---

## 📂 Script: `CorticalPlaneProfExtractor.m`

### Overview of Functionality

This script performs the following steps:

1. **Input**
   - Requires a path to the directory containing cortical band images.

2. **Data Organisation**
   - Loads image matrices and organises them into a data structure named `cortex`, based on the folder hierarchy.
   - The folder structure can be customised depending on experimental design (e.g., genotypes, RNAi conditions).

3. **AP Intensity Profile Extraction**
   - Moves a sliding window of a defined width (`ROI_width`) across the image band, advancing one pixel at a time.
   - At each position, computes the average intensity of pixels within the window.
   - The resulting AP intensity profile is stored within the `cortex` data structure.

4. **Standardisation via Percentile Binning**
   - To account for variability in cortex length between embryos, each AP profile is divided into equal-sized bins based on percentiles (e.g., 0–100% of cortex length).
   - The binned profiles are saved as `percentile_cortex` within `cortex`.

5. **Grouping by Experimental Condition**
   - Constructs a second data structure, `s2`, that groups `percentile_cortex` profiles according to experimental condition.
   - Computes the mean, standard deviation (stdev), and standard error of the mean (SEM) for each condition.
   - All results are stored in `s2`.

6. **HiLo Ratio Filtering and Plotting**
   - Applies a HiLo (High/Low) ratio filter to normalise each AP intensity profile.
   - For each experimental condition, calculates the mean, stdev, and SEM of the HiLo-normalised profiles.
   - Plots for each condition the mean ± SEM of the HiLo-normalised AP intensity profiles.

---

## 📌 Notes

- Cortical bands are tif images processed in FIJI
- The script assumes that images are aligned along the AP axis and that their vertical size (height) is fixed at 100 pixels.
- Adapt the folder-reading section to match your specific experiment naming or structuring conventions.
- Suitable for batch processing of embryos and comparative analysis between experimental conditions.

---

## 🧪 Requirements

- MATLAB (R2020 or later recommended)
- Image Processing Toolbox (for image loading and manipulation)

---

## 📁 Example Folder Structure

/Data
├── Condition1
│ ├── embryo1.tif
│ └── embryo2.tif
└── Condition2
├── embryo1.tif
└── embryo2.tif