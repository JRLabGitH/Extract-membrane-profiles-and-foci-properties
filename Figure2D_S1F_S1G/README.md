Embryo Image Analysis
This MATLAB script processes embryo microscopy images to compute normalized profiles and mean anterior values. The results are saved in a summary table.

Steps Overview
Step 1: Find the Images
Searches the current directory for .tif image files.
Step 2: Loop Through and Manipulate Images
Opens each image and determines its halfway point.
Splits the image into two halves:
The first half is stored as norm_ratiol.
The second half is flipped horizontally and stored as norm_ratior.
Stores these halves in structures DL and DR with filenames prefixed by "L_" and "R_", respectively.
Concatenates DL and DR into structure D.
Step 3: Pre-allocate Structures
Initializes norm_profs to store normalized profiles.
Pre-allocates anterior_means as a column vector of NaNs.
Pre-allocates embryo_names as a column vector of strings.
Step 4: Generate Profiles
For each profile:
Extracts the cortex section and averages its top 10 values.
Computes the cytoplasm mean value.
Normalizes the cortex mean by the cytoplasm mean.
Stores the normalized profiles and corresponding filenames.
Step 5: Compute Mean Anterior Values
If norm_prof length ≥ 100: computes mean for indices 20-100.
If length is between 20 and 100: computes mean for indices 20 to end.
If length < 20: sets value to NaN.
Step 6: Merge and Save Data
Creates a EmbryoSummary table with EmbryoName and AnteriorMean.
Saves the table as EmbryoAnteriorMeans.mat.
Usage
Place this script in a directory containing .tif files.
Run the script in MATLAB.
The output file, EmbryoAnteriorMeans.mat, contains the summary table.
Requirements
MATLAB (R2017a or later recommended).
Image files in .tif format.
Contributing
Contributions are welcome! Open an issue or submit a pull request for improvements.

License
This project is licensed under the MIT License. (Adjust according to your preference.)
