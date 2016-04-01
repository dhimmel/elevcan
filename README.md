# A data exploration of elevation, oxygen & lung cancer

This repository hosts the analysis for:

> Simeonov KP, Himmelstein DS (2015) [**Lung cancer incidence decreases with elevation: evidence for oxygen as an inhaled carcinogen**](). *PeerJ* 2:e705 DOI: 10.7717/peerj.705

###Directory Structure

+ `code/`: All code that produced final project output. Run `Rscript code/run.R` (from the project's root directory) to perform the analysis starting from the county-level dataset.

+ `code/data-creation/`: contains code that was used to create the county-level dataset including resource parsing, data integration, and population-weighted elevaton computation.

+ `output/`: output from running the analysis including the stdout log, underlying figure data, and persistent copy of the final R session.

+ `tables/`: Latex tables generated from the analysis. Additional manual editing was performed on some tables.

+ `figures/`: Figures as pdfs generated from the analysis. Additional manual editing was performed for some figures.

+ `data/`: text files of the county-level dataset.

+ `manual/`: Contains manually editted files (files that are not produced by running the analysis). Final versions of the figures are available as pdfs. Saved sessions from propreitary software used for image processing are also included.

