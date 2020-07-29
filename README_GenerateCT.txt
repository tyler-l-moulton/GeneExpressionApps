Generate qPCR DATA -- README 

SUMMARY:

This app generates data for simulated qPCR experiments. Colors on plots represent different subjects or individuals in the experiment. The simulated data are based on a logistic function of the form Fluorescence = Alpha / (1 + Beta * e^(-Gamma*cycle)). The user may adjust parameters alpha, beta, and gamma to alter the shape of the curve. The user may also establish the number of PCR cycles, replicate trials, number of subjects, and may translate the curve along the x-axis by specifying "Shift CT". The user controls variability of the data. Simulated data are previewed at the bottom of the app page and can be downloaded in a ".csv" format. Plots may be downloaded as an image by right-clicking (or cmnd+Click). Only handles one gene and one treatment at a time




INPUT:

LOGISTIC CURVE PARAMETERS

Alpha			Fluorescence = Alpha / (1 + Beta * e^(-Gamma*cycle))

Beta			Fluorescence = Alpha / (1 + Beta * e^(-Gamma*cycle))

Gamma			Fluorescence = Alpha / (1 + Beta * e^(-Gamma*cycle))


qPCR PARAMETERS

Number of Cycles	Specify number of cycles for qPCR assay

qPCR Replicates		Specify number of replicates per subject

Number of Subjects	Specify number of subjects in assay

Shift CT		Translate curve along the x-axis in units of "cycle"



ADJUST VARIABILITY

Noise			Implements a random jittering function to simulate the varition in 
			fluorescence readings (random variation along y-axis). 
			NOTE: This jittering does NOT correspond linearly to fluorescence units.

SD of CT among Subjects	Controls variability in the value of CT along the x-axis, in "cycle" units.
			This randomization is implemented before the "Noise" randomization.

Treatment		Specify a name for the experimental treatment or species (dependent variable)
			for comparisons of relative expression.

Gene			Specify a name for the gene (dependent variable) for comparisons of
			relative expression.

BUTTONS

Update Data		This action button implements the data generation routine with the variables
			specified above. NOTE: Randomization algorithms are implemented again.

Download Data		Download a ".csv" file with the name specified in "File Name" window.
			NOTE: Do NOT include "." or file extensions in "File Name".


OUTPUT

PLOTS			2 plots are simultaneously created and stacked to form one downloadable image

qPCR Fluorescence	Curves for each subject and replicate, fluorescence vs. cycle number. 
			Colors code for different subjects.

log10 Transformed	Curves for each subject and replicate, log10(fluor) vs. cycle number.
qPCR Fluorescence	Colors code for different subjects.

DATA			Data table is displayed at the bottom of the app page. Displays the first 
			100 rows of the simulated data with the following headings:
			"subject", "qPCR.Trial", "cycle", "log10.fluor", "fluor", "treatment", "gene"




By Tyler Moulton (5/18/2020)

