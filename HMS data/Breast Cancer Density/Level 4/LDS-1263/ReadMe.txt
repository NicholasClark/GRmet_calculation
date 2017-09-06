################################################################
Dataset Information
################################################################

Dataset Name:
Breast Cancer Cell Density: Imaging assay of the density- and context-dependence of small molecule perturbagen response in breast cancer cell lines. Dataset 3 of 3: Calculated dose response metrics.

Dataset Description:
The density- and context-dependent sensitivities of 6 cell lines plated at 6 different cell densities were assayed by measuring cell counts 72 hours after treatment with each of 12 small molecule inhibitors across a 9-point dose range. This dataset presents a comparison of the dose response metrics calculated using the normalized and non-normalized cell count data reported in HMS LINCS Datasets #20256 and 20257.

--Files in Package:
20258.txt
Cell_Line_Metadata.txt
MetadataPropertyMapping.txt
ReadMe.txt
Small_Molecule_Metadata.txt

################################################################
Center-specific Information
################################################################

Center-specific Name:
HMS_LINCS

Center-specific Dataset ID:
20258

Center-specific Dataset Link:
http://lincs.hms.harvard.edu/db/datasets/20258/

################################################################
Assay Information
################################################################

Assay Protocol:
1. For each cell line, cells were plated at the indicated density in 60 uL of growth media in quadruplicate 384-well plates. The cells then were grown for 24 hours. For growth conditions, see PMC3845839.<br />
2. Three plates per cell line were treated with the indicated small molecules by direct dispensing of DMSO stock solutions to the indicated concentrations into 60 µL of media using an HP D300 Digital Dispenser.<br />
3. Immediately after treatment of these three plates, the untreated fourth plate for each cell line was fixed and stained by adding 20 ¼L of fixative containing formaldehyde solution and Hoechst 33342 at a final concentration of 3% formaldehyde and 250 ng/mL Hoechst. This plate was sealed and kept at 4°C until day 3 of the assay.<br />
4. On day 3 the treated plates were fixed, stained, and sealed as described in step 3.<br />
5. After 1 hour all plates were scanned with a PE Operetta high-throughput plate scanner.<br />
6. Nuclei were counted using Columbus software.<br />
7. Nuclei counts were normalized to DMSO-treated controls at the same plating density and on the same plate to yield relative cell count and normalized growth rate inhibition (GR) values for each plate (technical replicate) for each cell line/density/small molecule/concentration combination (HMS LINCS Datasets #20256). Relative cell count is the measured cell count for a given treatment divided by the 50%-trimmed mean of the cell count of the DMSO-treated control wells on the same plate. Normalized growth rate inhibition values were calculated according to the following formula: 2^[log2(x(c)/x0)/log2(xctrl/x0)]-1 where x(c) is the measured cell count after a given treatment, x0 is the 50%-trimmed mean of the cell count from the day 0 untreated plate grown in parallel until the time of treatment, and xctrl is the 50%-trimmed mean of the cell count of the DMSO-treated control wells on the same treated plate.<br />
8. Second, the results from the three treated plates (technical replicates) were averaged to yield the mean relative cell count and the mean normalized growth rate inhibition (GR) value for each cell line/density/small molecule/concentration combination (HMS LINCS Datasets #20257).<br />
9. For each cell line/density/small molecule combination, mean normalized growth rate inhibition (GR) values across all tested concentrations were fitted to a sigmoidal curve to extract the GR<sub>50</sub>, GR<sub>max</sub>, AUC, EC<sub>50</sub>, GR<sub>inf</sub>, and Hill coefficient (HMS LINCS Dataset #20258). For comparison purposes only, mean relative cell counts across the tested concentrations for each cell line/density/small molecule combination also were fitted to a sigmoidal curve to extract the IC<sub>50</sub>, E<sub>max</sub>, AUC, EC<sub>50</sub>, E<sub>inf</sub>, and Hill coefficient (HMS LINCS Dataset #20258).

Date Updated:
2016-08-19

Date Retrieved from Center:
11/13/2015

################################################################
Metadata Information
################################################################

Metadata information regarding the entities used in the experiments is included in the accompanied metadata. A metadata file per entity category is included in the package. For example, the metadata for all the cell lines that were used in the dataset are included in the Cell_Lines_Metadata.txt file.
Descriptions for each metadata field can be found here: http://www.lincsproject.org/data/data-standards/
[/generic/datapointFile]
[/generic/reagents_studied]
