# original-algorithm-paper-corrections

This repository contains corrected supplemental files for the paper "Automated identification of implausible values in growth data from pediatric electronic health records" in the Journal of the American Medical Informatics Association, Volume 24, Issue 6, November 2017, Pages 1080–1087 by Daymont et al. 

https://academic.oup.com/jamia/article/24/6/1080/3767271

The JAMIA paper describes the development of what has been updated and is now referred to as the growthcleanr pediatrics algorithm. The JAMIA publication included five supplemental files/folders, four of which have required updates because of small errors. Further details about the supplemental files and updates are below.

Supplemental File 1: Updated code available on this repository. The original Stata code. There was an error in the EWMA code (5 was added to the difference in age in days before the absolute value was taken, rather than after which was correct), which was corrected. There was also an error in step 12ei regarding dealing with what were then referred to as "duplicates" (now referred to as same-day-extraneous measurements) for a subject with no nonduplicate values and no values for the other parameter. An alternate method of addressing this issue is presented in the updated Stata code available on this repository. In testing, these errors had minimal impact on data cleaning results. Please note that there is one aspect of the Stata code (relating to macros for dealing with both growth parameters at one time) that still works on Stata on Macs but does not work on Stata on Windows. Importantly, note that there have been other updates to the algorithm that are available only in the R code (see Supplemental File 2 below). The edits to the do file only address the two errors that were identified.

Supplemental File 2: Updated code available at https://github.com/carriedaymont/growthcleanr. The supplemental file was the original R code. The main growthcleanr github repository now has updated R code. Errors were corrected in step 12eii and step 16; both of which had minimal effect on results of tested datasets. Additionally, although the pediatric algorithm is largely the same, changes have been made to improve performance on large datasets, an adult algorithm has been added, and other tweaks have been made. Additional documentation is available on the main growthcleanr repository, and any further updates will also be available there. 

Supplemental File 3: Updated files (clean and tracked changes) available on this repository. This is the original algorithm specifications or "English instructions" that were used to translate the algorithm from Stata into R. The corrected method for addressing subjects with no nonduplicate values and no measurements for the other parameter in step 12ei are described in this updated document. As with supplemental file 1, the additional updates that have been made to the R code are not addressed in this supplemental file. 

Supplemental File 4: Updated files available at https://github.com/carriedaymont/growthcleanr. There were errors in three total rows of the two WHO height velocity tables. Corrected versions are now available at the growthcleanr repository. These tables affect one step of the algorithm, and we identified very small differences when growthcleanr was run with the original vs updated tables. 

Supplemental File 5: No updated files. This document described one component of validation of the algorithm. 


