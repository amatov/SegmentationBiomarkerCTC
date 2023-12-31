### I designed the whole algorithm for the segmentation of prostate-specific membrane antigen labeling in ciruculating tumor cells of metastatic patients and wrote code in Matlab to demonstrate a proof of princip of this analysis 

### The final Matlab implementation was done by Shayan Modiri (the group of Mubarak Shah)

#### The initial image segmentation step is accomplished by stationary wavelet transform, which identifies bright pixel clusters in noisy images; the seeding step. Active contour, next, identifies precisely the edges of the image features based on the seeds. Watershed transformation of the seeding step image is overlaid, with reversed intensities, on the active contour image. Logical conjunction / “and” of the active contour image and the watershed image identifies the bright areas and their exact borders.  

#### For an image with an example, see Figure 2 of the DoD application here: https://researchgate.net/publication/374544332_Automated_Enumeration_and_Analysis_of_Circulating_Tumor_Cells_from_Peripheral_Blood_of_Metastatic_Cancer_Patients_for_Diagnosis_and_Refinement_of_Therapy_2012-2015
