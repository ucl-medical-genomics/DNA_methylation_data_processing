***

# DNA Methylation Data Preprocessing Summary

***Simone Ecker***
*[s.ecker@ucl.ac.uk](mailto:s.ecker@ucl.ac.uk)*\
***Medical Genomics, UCL Cancer Institute, London, UK***

*Last updated: 13 September 2019*

***

## 1 Main R packages

The main R packages used for DNA methylation data preprocessing are:

*   `Minfi`<sup>1</sup>
*   `ChAMP`<sup>2,3</sup>
*   `RnBeads`<sup>4</sup>

## 2 Data import

IDAT files obtained from Illumina’s iScan are directly imported into R as raw data using `minfi`’s<sup>1,5</sup> read function. The corresponding sample sheet and, if necessary, a specific array manifest file are also read in. Data are then converted into an `RGChannelSet` containing the raw intensities in the green and red channels as well as the intensities of internal control probes.

## 3 Quality control

### 3.1 Initial quality control

For quality control in `minfi`, the `RGChannelSet` is converted into a `MethylSet` to convert the red and green intensities to methylated and unmethylated signals according to the microarray’s probe design. This does not perform any normalization. Then, the following quality control plots are produced:

*   Log median intensity in methylated and unmethylated channels via `getMeth()` and `getQC()`
*   Densities: `densityPlot()`
*   Density bean plots: `densitybeanPlot()`

### 3.2 Quality assessment using control probes

`Minfi` can produce several control probes plots for quality assessment: 

*   Control probes plot for bisulfite conversion (for probe type I and II): `controlStripPlot()`
*   All these plots can be generated along with further informative control probes plots in one step using the function `qcReport()`

In addition, quality checks as implemented in Illumina’s BeadArray Controls Reporter<sup>6</sup> can be performed with an R package `ewastools` published by Heiss _et al._<sup>7</sup> which evaluates 17 control metrics calculated from control probes. The following metrics are calculated and reported back in numbers and figures can also be generated (most of them are also shown in a different style in `minfi`’s qcReport):

*   Restoration
*   Staining green
*   Staining red
*   Extension green
*   Extension red
*   Hybridization high/medium
*   Hybridization medium/low
*   Target removal 1
*   Target removal 2
*   Bisulfite conversion I green
*   Bisulfite conversion I red
*   Bisulfite conversion II
*   Specificity I green
*   Specificity I red
*   Specificity II
*   Non-polymorphic green
*   Non-polymorphic red

Although less informative than the above, control probes plots as provided by Illumina GenomeStudio<sup>8</sup> can be obtained with the R package `HumMethQCReport`<sup>9</sup> (currently not available for EPIC). It generates the following figures: 

*   Barplots of Illumina internal controls
*   Intensity at high and low Beta-values
*   Percentage of non-detected CpGs
*   Average detection p-value

### 3.3 Identification of sample mismatches

The `ewastools` package further provides the option to identify potential sample mix-ups using high-frequency SNP probes. RnBeads<sup>4</sup> also provides this functionality. For the same purpose of detecting possible mismatches, the sex of the sample donor can be determined based on the methylation signal at sex chromosomes. Sample contamination can also be assessed with `ewastools`.

### 3.4 Quality assessment based on all cg probes used for downstream analyses

A `RatioSet` object can be generated to store actual methylation values instead of the methylated and unmethylated signals (a copy number matrix can also be obtained if required). The `RatioSet` is generated with the function `ratioConvert()`. Subsequently, both a Beta- and a M-value matrix are obtained from the `RatioSet` via `getBeta()` and `getM()`. Alternatively, the methylation value matrices can also be obtained from an `RGChannelset` or `MethylSet`.

The following plots for quality control are produced from both the Beta- and M-value matrix: 

*   Sample-wise density plots colored by phenotype of interest 
*   Group-wise density plots (CpG-wise mean per group)
*   Principal component analysis (PCA) using all CpGs, colored by phenotype of interest
*   Multidimensional scaling (MDS) using all CpGs, colored by phenotype of interest
*   Hierarchical clusterings (complete, average, median, wardD, centroid and single) performed on Euclidean distances, including a heatmap and dendrograms, colored by phenotype of interest
*   PCA, MDS and clusterings as above, but colored by all additional covariates to be assessed, in particular, technical (or biological) variables that could potentially lead to batch effects
*   Singular Value Decomposition (SVD) plots examining associations between principal components and all technical and biological covariates available
*   Age distribution per phenotype of interest (overlayed histograms or density plots)

## 4 Filtering

### 4.1 Sample filtering

Based on the above quality assessment, samples that turn out to be outliers corresponding to multiple measurements and plots might be removed from the dataset for downstream analyses.

### 4.2 Probe filtering

Probes meeting the following criteria are removed from all samples:

*   Detection p-values > 0.01 in ≥ 1 sample or more stringent thresholds as recently suggested<sup>10,11</sup>
*   Bead count < 3 in ≥ 5% of the samples
*   Ambiguous genomic locations (cross-hybridization) as provided by Zhou _et al._<sup>12</sup>
*   SNPs of the corresponding population provided by Zhou _et al.<sup>12_</sup>
*   Non-cg probes
*   Sex chromosomes (unless they are needed for specific analyses)
*   Non-variable probes (optional, only if required for specific reasons)

## 5 Normalization

Bisulfite conversion of DNA and other experimental steps introduce assay variability and batch effects. Technical variation should be reduced as much as possible. Illumina 450K and EPIC BeadChips use two fluorescent dyes (Cy3-green and Cy5-red) and two different probe types (Infinium I and Infinium II), leading to technical biases which should be corrected. Background noise also needs to be removed. To this aim, negative control probes and out-of-band (OOB) probes are present on Illumina’s arrays. 

Many different normalization methods are available and the normalization method to be applied depends on each dataset, the research question and the downstream analyses to be performed. Usually, suitable normalization methods are performed (on either Beta- or M-values as required by the method), and their results are compared to each other. Then, the quality assessment steps listed in section 3.4 are repeated to help to choose the most appropriate method(s) for the specific dataset and research question. Individual steps from different normalization methods can also be combined and performed sequentially to obtain optimal results. 

For example, normalization methods that perform (only) specific tasks such as dye bias, probe type bias or background correction may be followed by normalization methods that remove technical variance across arrays. Several normalization methods combine the correction of the above-mentioned technical biases and subsequent global normalization including both within- and across-sample normalization to remove unwanted technical variance. Across-sample normalization may not always be suitable<sup>13</sup>. A few methods that only apply within-sample normalization are available such as for example PBC and ssNoob. These normalization methods are particularly useful for datasets where the global distribution of DNA methylation levels is expected to differ significantly across samples<sup>13</sup> and for projects where the normalization of single samples needs to be reproducible independently from the rest of the samples, or where new additional samples will be added to the dataset over time<sup>5</sup>.

*   **Illumina<sup>1</sup>** (available in `minfi`):  \
reverse engineered Illumina GenomeStudio<sup>8</sup> control normalization with background subtraction
*   **Lumi<sup>14</sup>** (available in `methylumi`): \
background correction, variance stabilization and normalization
*   **PBC<sup>15</sup>** (<ins>P</ins>eak-<ins>B</ins>ased <ins>C</ins>orrection, available in `ChAMP`): \
intra-array probe type bias correction based on global profiles
*   **Quantile<sup>1,16</sup>** (available in `minfi`): \
global quantile normalization performing both sample normalization and probe type correction
*   **Quantile stratified<sup>1,16</sup>** (available in `minfi`): \
as above, but with quantile normalization performed within genomic region strata
*   **SWAN<sup>17</sup>** (<ins>S</ins>ubset-quantile <ins>W</ins>ithin <ins>A</ins>rray <ins>N</ins>ormalization, available in `ChAMP`): \
reduces technical variation within and between arrays
*   **BMIQ<sup>18</sup>** (<ins>B</ins>eta-<ins>MI</ins>xture <ins>Q</ins>uantile normalization, available in `ChAMP`): \
intra-array normalization which adjusts for probe type bias
*   **Noob<sup>19</sup>**(<ins>N</ins>ormal-exponential <ins>O</ins>ut-<ins>O</ins>f-<ins>B</ins>and, available in `minfi`): \
background correction method with dye-bias normalization
*   **Data-driven normalization<sup>20</sup>** (available in `wateRmelon`): \
background correction, between-array quantile normalization of methylated and unmethylated signal intensities and the two probe types separately (for ‘dasen’ normalization, for example)
*   **FunNorm<sup>21</sup>** (<ins>Fun</ins>ctional <ins>Norm</ins>alization, available in `ChAMP`): \
between-array quantile normalization of methylated and unmethylated signal intensities and the two probe types separately, uses control probes to remove unwanted technical variation
*   **ENmix<sup>22</sup>**(<ins>E</ins>xponential-<ins>N</ins>ormal <ins>mix</ins>ture signal intensity background correction, `ENmix `package):** \
**background correction, probe type (RCP<sup>23</sup>, <ins>R</ins>egression on <ins>C</ins>orrelated <ins>P</ins>robes) and dye bias correction (RELIC<sup>24</sup>, <ins>Re</ins>gression on <ins>L</ins>ogarithm of <ins>I</ins>nternal <ins>C</ins>ontrol probes), inter-array normalization
*   **ssNoob<sup>5</sup>** (<ins>s</ins>inge-<ins>s</ins>ample <ins>N</ins>ormal-exponential <ins>o</ins>ut-<ins>o</ins>f-<ins>b</ins>and, available in `minfi`): \
intra-array normalization suitable for incremental preprocessing of individual samples and when integrating data from multiple generations of Illumina microarrays
*   **SeSAMe<sup>25</sup>**(<ins>SE</ins>nsible <ins>S</ins>tep-wise <ins>A</ins>nalysis of <ins>Me</ins>thylation data, `SeSAMe` package): \
probe quality masking with `pOOBAH` (<ins>p</ins>-value with <ins>O</ins>ut-<ins>O</ins>f-<ins>B</ins>and probes for <ins>A</ins>rray <ins>H</ins>ybridization), bleed-through correction in background subtraction, non-linear dye bias correction, control for bisulfite conversion, reduction of inter- and intra-array technical variation

In our experience, SWAN, BMIQ and Funnorm together with ENmix's probe type correction give good results with successful reduction of technical variation and correction for background and probe type bias. Dye-bias normalization can be achieved with Noob, for example.

The mentioned methods provide appropriate results for analyses of global DNA methylation and variability patterns, rank-based approaches, site- and region-based analysis of differential mean and variability of DNA methylation, and epigenetic clock analyses. Probe type bias correction is particularly important for the assessment of global DNA methylation patterns, clusterings, rank-based approaches and region-based analyses. The best possible removal of all technical variation is crucial under all circumstances, but even more critical for analyses of biological variation such as epigenetic drift and differential variability. 

SWAN works very well when no global shifts in DNA methylation are expected while BMIQ is better suited when significant global differences are present between the phenotypes of interest (see above). This can be the case when comparing normal samples with samples of a cancer type in which global hypomethylation occurs, for example. The quantro method<sup>13</sup> can aid the decision by assessing such global differences.

Quantile and quantile stratified normalization are less well suited for DNA methylation analyses and can introduce additional bias. Illumina’s GenomeStudio normalization method does not reduce technical variation<sup>26,27</sup>. It can even introduce additional variability<sup>28</sup>. FunNorm and BMIQ can also lead to biased distributions, in particular for DNA methylation M-values. SWAN, BMIQ and PBC have been shown to outperform other methods<sup>26,28,29</sup>. PBC can lead to discontinuities in probe type II density distributions under certain circumstances<sup>29</sup> and might thus not always be the best suitable approach. For EPIC, we have achieved very good results with Funnorm in combination with ENMix's probe type bias correction. All principally suitable methods should be applied to the data, and their results should be assessed carefully (see above).

## 6 Batch-effect correction

Useful methods to detect batch effects and biological confounders are **Surrogate Variable Analysis (SVA)** implemented in the `sva` package<sup>30</sup> and **Singular Value Decomposition (SVD)<sup>31</sup>** available in `ChAMP`. These methods assess the impact of significant components of variation present in a dataset. Unknown sources of variation can be detected, and importantly, the contribution of all available technical and biological covariates can and should be assessed. One 450K or EPIC BeadChip accommodates multiple samples (n=12 and n=8, respectively) which can cause a batch effect associated with array ID (‘Sentrix ID’). The arrays’ position on the BeadChip (row and column, ‘Sentrix Position’) also impacts the data. If any significant batch effects or unwanted confounders are observed, the data can be corrected for these effects by the application of **ComBat<sup>32</sup>** as long as they are not directly confounded with the phenotype of interest.

ComBat is the most widely used method for batch effect correction. It is implemented in `sva`, `minfi`, `ChAMP` and `RnBeads` and usually works best on DNA methylation M-values. However, sometimes better results are achieved when applying ComBat to Beta-values instead. Thus, it can be useful to perform ComBat on both M-values and Beta-values and assess the results with the help of SVA/SVD and the quality assessment plots described in section 3.4 to compare them to each other.

## 7 Subpopulation correction

When working with composite samples such as whole blood, peripheral blood mononuclear cells (PBMCs), or brain tissue, cell composition differences across samples may need to be corrected. The most common reference-based method to assess and correct for subpopulation differences is **Houseman’s algorithm<sup>33</sup>** in combination with a blood cell reference dataset of the six most common white blood cell types<sup>34</sup>, available in `minfi`, `ChAMP` and `RnBeads`. Reference-free methods also exist<sup>35–39</sup>. The results of the correction should be assessed visually by plotting the cell type distributions per sample, mean cell type distributions per group, and by repeating the quality assessment plots listed in section 3.4. 

Two very interesting new methods to identify the cell types driving differential DNA methylation signal in composite samples such as whole blood have recently been developed: **CellDMC** (DMC = <ins>D</ins>ifferentially <ins>M</ins>ethylated <ins>C</ins>ytosines)<sup>40</sup> and **TCA** (<ins>T</ins>ensor <ins>C</ins>omposition <ins>A</ins>nalysis)<sup>41</sup>.

## 8 Beta- and M-value conversion

As mentioned, some methods require Beta-values as input. M-values can be converted into Beta-values and vice versa using the following functions implemented in `RnBeads:`

```
  beta2mval <- function(betas, epsilon = 0.00001) {
    if (!is.numeric(betas)) {
      stop("invalid value for betas")
    }
    if (!(is.numeric(epsilon) && length(epsilon) == 1 && (!is.na(epsilon)))) {
      stop("invalid value for epsilon")
    }
    if (epsilon < 0 || epsilon > 0.5) {
      stop("invalid value for epsilon; expected 0 <= epsilon <= 0.5")
    }
    betas[betas < epsilon] <- epsilon
    betas[betas > (1 - epsilon)] <- 1 - epsilon
    return(log2(betas / (1 - betas)))
  }

  mval2beta <- function(mvals) {
    if (!is.numeric(mvals)) {
      stop("invalid value for mvals")
    }
    intmf<-2^(mvals)
    return(intmf/(intmf+1))
  }

```

## 9 Additional notes

M-values are used for all downstream statistical analyses such as differential mean methylation or analyses of variability and epigenetic drift. Beta-values are only used for the visualization of DNA methylation patterns. 

An exception are epigenetic clocks to perform epigenetic age estimations<sup>42–46</sup>. These methods require Beta-values as input. When working with whole blood data or PBMCs, the epigenetic clocks usually work better on data that were not corrected for subpopulation effects.

## References

1.	Aryee, M. J. _et al._ Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. _Bioinformatics_ **30,** 1363–9 (2014).

2.	Morris, T. J. _et al._ ChAMP: 450k Chip Analysis Methylation Pipeline. _Bioinformatics_ **30,** 428–430 (2013).

3.	Tian, Y. _et al._ ChAMP: updated methylation analysis pipeline for Illumina BeadChips. _Bioinformatics_ 2014–6 (2017). doi:10.1093/bioinformatics/btx513

4.	Assenov, Y. _et al._ _Comprehensive analysis of DNA methylation data with RnBeads_. _Nat Methods_ **11,** (2014).

5.	Fortin, J. P., Triche, T. J. & Hansen, K. D. Preprocessing, normalization and integration of the Illumina HumanMethylationEPIC array with minfi. _Bioinformatics_ **33,** 558–60 (2017).

6.	Illumina Inc. BeadArray Controls Reporter: Software Guide. (2015). doi:10.17485/IJST/2016/V9I47/94892

7.	Heiss, J. A. & Just, A. C. Identifying mislabeled and contaminated DNA methylation microarray data: An extended quality control toolset with examples from GEO. _Clin Epigenetics_ **10,** (2018).

8.	Illumina Inc. _GenomeStudio Methylation Module v1.8 User Guide_. _Illumina_ (2008).

9.	Mancuso, F. M., Montfort, M., Carreras, A., Alibés, A. & Roma, G. HumMeth27QCReport: An R package for quality control and primary analysis of Illumina Infinium methylation data. _BMC Res Notes_ **4,** 546 (2011).

10.	Lehne, B. _et al._ A coherent approach for analysis of the Illumina HumanMethylation450 BeadChip improves data quality and performance in epigenome-wide association studies. _Genome Biol_ **16,** 37 (2015).

11.	Heiss, J. A. & Just, A. C. Guidance on filtering DNA methylation microarray probes by detection p-values. _bioRxiv_ (2018). doi:10.1101/245217

12.	Zhou, W., Laird, P. W. & Shen, H. Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. _Nucleic Acids Res_ **45,** e22 (2017).

13.	Hicks, S. & Irizarry, R. Quantro: a Data-Driven Approach To Guide the Choice of an Appropriate Normalization Method. _Genome Biol_ **16,** 117 (2015).

14.	Du, P., Kibbe, W. A. & Lin, S. M. lumi: A pipeline for processing Illumina microarray. _Bioinformatics_ **24,** 1547–8 (2008).

15.	Dedeurwaerder, S. _et al._ Evaluation of the Infinium Methylation 450K technology. _Epigenomics_ **3,** 771–784 (2011).

16.	Touleimat, N. & Tost, J. Complete pipeline for Infinium® Human Methylation 450K BeadChip data processing using subset quantile normalization for accurate DNA methylation estimation. _Epigenomics_ **4,** 325–41 (2012).

17.	Maksimovic, J., Gordon, L. & Oshlack, A. SWAN: Subset-quantile Within Array Normalization for Illumina Infinium HumanMethylation450 BeadChips. _Genome Biol_ **13,** R44 (2012).

18.	Teschendorff, A. E. _et al._ A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450 k DNA methylation data. _Bioinformatics_ **29,** 189–96 (2013).

19.	Triche, T. J., Weisenberger, D. J., Van Den Berg, D., Laird, P. W. & Siegmund, K. D. Low-level processing of Illumina Infinium DNA Methylation BeadArrays. _Nucleic Acids Res_ **41,** e90 (2013).

20.	Pidsley, R. _et al._ A data-driven approach to preprocessing Illumina 450K methylation array data. _BMC Genomics_ **14,** 293 (2013).

21.	Fortin, J.-P. _et al._ Functional normalization of 450k methylation array data improves replication in large cancer studies. _Genome Biol_ **15,** 503 (2014).

22.	Xu, Z., Niu, L., Li, L. & Taylor, J. A. ENmix: A novel background correction method for Illumina HumanMethylation450 BeadChip. _Nucleic Acids Res_ **44,** e20 (2016).

23.	Niu, L., Xu, Z. & Taylor, J. A. RCP: A novel probe design bias correction method for Illumina Methylation BeadChip. _Bioinformatics_ **32,** 2659–63 (2016).

24.	Xu, Z., Langie, S. A. S., De Boever, P., Taylor, J. A. & Niu, L. RELIC: A novel dye-bias correction method for Illumina Methylation BeadChip. _BMC Genomics_ **18,** 4 (2017).

25.	Zhou, W., Triche, T. J., Laird, P. W. & Shen, H. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. _Nucleic Acids Res_ **46,** e123 (2018).

26.	Marabita, F. _et al._ An evaluation of analysis pipelines for DNA methylation profiling using the illumina humanmethylation450 BeadChip platform. _Epigenetics_ **8,** 333–46 (2013).

27.	Yousefi, P. _et al._ Considerations for normalization of DNA methylation data by Illumina 450K BeadChip assay in population studies. _Epigenetics_ **8,** 1141–52 (2013).

28.	Wu, M. C. _et al._ A systematic assessment of normalization approaches for the Infinium 450K methylation platform. _Epigenetics_ **9,** 318–29 (2014).

29.	Wang, T. _et al._ A systematic study of normalization methods for Infinium 450K methylation data using whole-genome bisulfite sequencing data. _Epigenetics_ **10,** 662–9 (2015).

30.	Leek, J. T., Johnson, W. E., Parker, H. S., Jaffe, A. E. & Storey, J. D. The SVA package for removing batch effects and other unwanted variation in high-throughput experiments. _Bioinformatics_ **28,** 882–3 (2012).

31.	Teschendorff, A. E. _et al._ An epigenetic signature in peripheral blood predicts active ovarian cancer. _PLoS One_ **4,** e8274 (2009).

32.	Johnson, W. E., Li, C. & Rabinovic, A. Adjusting batch effects in microarray expression data using empirical Bayes methods. _Biostatistics_ **8,** 118–27 (2007).

33.	Houseman, E. A. _et al._ DNA methylation arrays as surrogate measures of cell mixture distribution. _BMC Bioinformatics_ **13,** (2012).

34.	Reinius, L. E. _et al._ Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. _PLoS One_ **7,** (2012).

35.	Houseman, E. A. _et al._ Reference-free deconvolution of DNA methylation data and mediation by cell composition effects. _BMC Bioinformatics_ **17,** 259 (2016).

36.	Kaushal, A. _et al._ Comparison of different cell type correction methods for genome-scale epigenetics studies. _BMC Bioinformatics_ **18,** 216 (2017).

37.	Zou, J., Lippert, C., Heckerman, D., Aryee, M. & Listgarten, J. Epigenome-wide association studies without the need for cell-type composition. _Nat Methods_ **11,** 309–11 (2014).

38.	Rahmani, E. _et al._ Sparse PCA corrects for cell type heterogeneity in epigenome-wide association studies. _Nat Methods_ **13,** 443–445 (2016).

39.	Rahmani, E. _et al._ Correcting for cell-type heterogeneity in DNA methylation: a comprehensive evaluation. _Nat Methods_ **14,** 218–219 (2017).

40.	Zheng, S. C., Breeze, C. E., Beck, S. & Teschendorff, A. Identification of differentially methylated cell-types in Epigenome-Wide Association Studies. _Nat Genet Methods_ **15,** 1059–66 (2018).

41.	Rahmani, E. _et al._ Cell-type-specific resolution epigenetics without the need for cell sorting or single-cell biology. _bioRxiv_ (2018). doi:10.1101/437368

42.	Horvath, S. DNA methylation age of human tissues and cell types. _Genome Biol_ **14,** R115 (2013).

43.	Hannum, G. _et al._ Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. _Mol Cell_ **49,** 359–67 (2013).

44.	Levine, M. E. _et al._ An epigenetic biomarker of aging for lifespan and healthspan. _Aging (Albany NY)_ **10,** 573–91 (2018).

45.	Horvath, S. _et al._ Epigenetic clock for skin and blood cells applied to Hutchinson Gilford Progeria Syndrome and ex vivo studies. _Aging (Albany NY)_ **10,** 1758–75 (2018).

46.	Lu, A. T. _et al._ DNA methylation GrimAge strongly predicts lifespan and healthspan. _Aging (Albany NY)_ (2018).
