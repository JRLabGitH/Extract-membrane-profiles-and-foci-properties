# MembraneCrossProfAndPeakVals.m

This script generates membrane-crossing intensity profiles and extracts peak membrane fluorescence values from straightened cortical images of *C. elegans* zygotes.

---

## 🔍 Overview

The script reads in 60-pixel-wide straighten stripes encompassing the zygote cell membrane, drawn clockwise starting and finishing at the anterior pole. It then:

- Extracts membrane-crossing intensity profiles across the anterior and posterior cortex, spanning from the cell exterior into the cytoplasm.
- Calculates peak membrane fluorescence intensity values at both the anterior and posterior profiles.

---

## 🗂 Input Requirements

- A **directory path** pointing to the folder containing straightened cortical images of zygotes.
- Images must follow a folder hierarchy that can be adapted to reflect experimental design (e.g. genotypes, RNAi conditions, etc.).

---

## 🧬 Processing Steps

1. **Load Images**  
   Loads image matrices from the given path and organises them into a structured array called `cortex` based on the folder hierarchy.

2. **Parse Folder Hierarchy**  
   Automatically assigns metadata (e.g. condition or genotype) based on folder structure. This can be adjusted for different experimental designs.

3. **Region Selection**  
   - Two **60 × 60 pixel anterior regions** are selected, starting 60 pixels from each end of the stripe.  
   - One **60 × 120 pixel posterior region** is selected at the centre of the stripe.  
   These regions are stored in the `cortex` structure.

4. **Normalisation**  
   Pixel intensities within each region are normalised by dividing each value by a local cytoplasmic intensity, calculated as the mean of the bottom 10 pixel rows. Normalised data are stored in `cortex`.

5. **Profile Generation**  
   Normalised regions are averaged along the x-axis to generate one-dimensional profiles crossing the membrane (from outside the embryo, through the membrane, and into the cytoplasm).  
   - Anterior regions are concatenated into a 60X120 area before averaging.  
   - All profiles are stored in `cortex`.

6. **Peak Detection and Representation**  
   The script identifies the first intensity peak ~30 pixels from the top of each profile.  
   - Anterior and posterior peak positions and values for each embryo are stored in `cortex` (for easy access anterior and posterior values are also stored in 'peakv' matrix).  
   - Use these data to export into GraphPad Prism for plotting.

7. **Profile Alignment and Plotting**  
   - Profiles were aligned to their peak intensity by padding with NaNs so that peak positions match across embryos.  
   - Centred profiles are stored in the matrix `ant_prof_center`.  
   - Profiles are then averaged by condition. Their mean, standard deviation, number of profiles, and peak position (for alignment to zero) are stored in the `prfCent` data structure.  
   - Use these data to export into GraphPad Prism for plotting.

---

## 📁 Output

- `cortex`: Structured array containing all processed regions, normalised data, membrane-crossing profiles, and peak intensity values.
- `prfCent`: Data structure with mean, standard deviation, and alignment info for averaged, centred anterior profiles.
- `ant_prof_center`: Matrix containing individual centred profiles.

---

## 🛠 Notes

- This script assumes that anterior is located at both ends of the straighten stripe, while posterior lies in the centre.
- Some manual adjustments may be required if images contain background noise or broken cortex traces.

---

## 🧑‍💻 Author

Originally written for the analysis in **Figure S2**. Adapt and extend as needed for similar membrane profile analyses.
