### I designed the whole algorithm for the segmentation of prostate-specific membrane antigen labeling in ciruculating tumor cells of metastatic patients and wrote code in Matlab to demonstrate a proof of principle of this analysis 

### The final Matlab implementation was done by Shayan Modiri (the group of Mubarak Shah, PhD)

#### The initial image segmentation step is accomplished by stationary wavelet transform, which identifies bright pixel clusters in noisy images; the seeding step. Active contour, next, identifies precisely the edges of the image features based on the seeds. Watershed transformation of the seeding step image is overlaid, with reversed intensities, on the active contour image. Logical conjunction / “and” of the active contour image and the watershed image identifies the bright areas and their exact borders.  

#### For an image with an example, see Figure 8 here: https://www.researchgate.net/publication/374059796_Personalized_Therapy_for_Sensitization_of_Resistant_Tumors_and_Identification_of_Putative_Targets_within_the_Microtubule_Transcriptome_2007_-_2017
