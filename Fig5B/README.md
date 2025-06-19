# 🔬 Foci Property Extraction from Time-Lapse Microscopy Data

This MATLAB script automates the extraction and analysis of **foci properties** (such as size, circularity, and intensity ratio) from time-lapse `.tif` images acquired in multiple experimental conditions. It is intended for use with spinning disk confocal microscopy data of *C. elegans* embryos and supports structured batch processing of different strains or conditions.

---

## 📁 Folder Structure

The script assumes the following folder organisation:

/ParentFolder
├── Strain1
│ ├── image1.tif
│ ├── image2.tif
│ └── /masks/
│ ├── mask1.tif
│ └── mask2.tif
├── Strain2
│ └── ...
└── ...


Each strain or condition has its own folder containing `.tif` time-lapse image stacks and a subfolder named `/masks/` with corresponding binary masks for each image.

---

## ⚙️ Script Overview

The script performs the following steps for each strain folder:

1. **Initialisation**
   - Loops through each subdirectory (strain or condition).
   - Prepares a structure (`allData`) to store analysis results.

2. **Image and Mask Handling**
   - Loads each time-lapse `.tif` stack and its corresponding binary mask.
   - Applies the mask to restrict analysis to the relevant cell region.

3. **Foci Detection and Measurement**
   - At selected timepoints (frames 10, 20, 30, 40, 50), thresholding is applied using `multithresh`.
   - Binary image processing (e.g. `bwareaopen`, `imopen`) cleans up noise and isolates foci.
   - Connected components (foci) are detected and labelled.
   - For each focus, the script measures:
     - **Area**
     - **Circularity**
     - **Raw intensity**
     - **Background intensity**
     - **Intensity ratio (focus/cytoplasm)**

4. **Data Aggregation**
   - Foci data is concatenated across all frames and stored in:
     - `frameFociProps` (per image stack)
     - `allFociProps` (per strain)

5. **Visualisation**
   - Optionally displays segmentation overlays and foci scatter plots.
   - Plots foci size distributions using `KDFplot`.

6. **Statistics and Export**
   - Calculates and stores:
     - Mean and standard deviation of foci size and intensity ratio
     - Number of foci detected
     - Mean number of foci per frame
   - Saves all results as an `.xlsx` file per strain.

7. **Group Comparison**
   - Uses `barErrorSigJitter` to compare:
     - Foci sizes
     - Filtered foci sizes (e.g. size < 3000 & circularity > 0.5)
     - Intensity ratios
   - Plots include mean ± error bars (likely SEM or 95% CI)

---

## 📦 Output

- **Excel spreadsheets** per strain with foci statistics
- **MATLAB figures** showing:
  - Foci size distributions
  - Group comparisons across strains
- **`allData` structure** in memory containing all measured properties

---

## 🧪 Requirements

- MATLAB (recommended: R2020 or later)
- Image Processing Toolbox
- Custom functions (must be available in MATLAB path):
  - `imEx1.m` – for reading image files
  - `KDFplot.m` – for kernel density function plots
  - `barErrorSigJitter.m` – for bar plots with jitter and significance testing

---

## 📌 Notes

- The script uses hardcoded paths (e.g. for saving Excel files). These should be edited to suit your file system.
- The analysis skips the first two entries in the `dir` output (`.` and `..`) by starting at `s = 3`.
- Circularity and intensity ratio thresholds can be adjusted for custom filtering.
- Make sure mask filenames partially match their corresponding image filenames.

---

## ✏️ Customisation Tips

- Change `tp=1:5` if you wish to analyse different or more frames.
- Modify `th=multithresh(...)` if your images require a different thresholding strategy.
- Edit the statistical plotting section to include additional comparisons or normalisations.

---

## 📫 Questions or Issues?

If you encounter problems or would like to request improvements, please open an issue on this repository.

