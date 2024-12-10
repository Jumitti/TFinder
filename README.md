<picture>
    <img
        src="./img/banners_TFinder.jpg">
</picture>

# TFinder 🧬🔍 [![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://tfinder-ipmc.streamlit.app/)

## Overview

TFinder is a Python easy-to-use web tool for identifying Transcription Factor Binding Sites (TFBS) and Individual Motif (IM). Using the NCBI API, it can easily extract either the promoter or terminal regions of a gene through a simple query of NCBI gene name or ID. It enables simultaneous analysis across five different species for an unlimited number of genes. The tool searches for TFBS and IM in different formats, including IUPAC codes and JASPAR entries. Moreover, TFinder also allows the generation and use of a Position Weight Matrix (PWM). Finally, the data may be recovered in a tabular form and a graph showing the relevance of the TFBSs and IMs as well as its location relative to the Transcription Start Site (TSS) or gene end. The results may be sent by email to the user facilitating the ulterior analysis and data sharing.

TFinder is written in Python and is freely available on GitHub under the MIT license: https://github.com/Jumitti/TFinder and can be accessed as a web application implemented in Streamlit at https://tfinder-ipmc.streamlit.app.

DOI: [https://doi.org/10.21203/rs.3.rs-3782387/v1](https://doi.org/10.21203/rs.3.rs-3782387/v1)

<picture>
    <img
        src="./img/tfinder_schema.png">
</picture>

## Description
Transcription factors (TFs) are proteins that bind to DNA to regulate gene expression. They specifically recognize a nucleotide sequence called a transcription factor binding site (TFBS) in the promoter and terminator regions of genes. The search of these TFBSs is an empirical discipline in the field of genomics that concerns a key step before TFBS functional validation by gel shift assays (EMSA) and chromatin immunoprecipitation (ChIP) that allow the examination of the interaction between a TF and DNA (Jayaram, Usvyat and R. Martin 2016)


The in-silico research of TFBS can be tedious and time-consuming at various stages, especially for novices in the discipline. Thus, first, it is necessary to retrieve the promoter or terminator nucleotide sequence of a gene. This step may be achieved by the utilization of several databases such as NCBI, UCSC and Ensembl, but they are not intuitive and user-friendly. Next, after identifying the promoter sequence of interest, one may use TF databases such as JASPAR (Castro-Mondragon et al. 2022) and TRANSFAC (Matys 2006), but they have their limitations. For example, these platforms do not allow the search of TFBS from an unreferenced TF and may be subject to a fee. Other tools such as PROMO (Farre 2003), TFBIND (Tsunoda and Takagi 1999), TFsitescan make it possible to find all the TFs binding to a nucleotide sequence; nevertheless, they all use JASPAR and TRANSFAC databases and do not allow use of personal TF and their TFBS. Moreover, these tools are rather archaic and not very user-friendly. There is only MEME that allows research with of your “own” TFBS (Bailey et al. 2015). MEME has a large tool library but is a niche software suite. FIMO is their most similar tool to TFinder (Grant, Bailey and Noble 2011).

TFinder is an ultra-intuitive, easy-to-use and fast analysis open source and free tool that allows both the retrieval and search of TFBS in a unique site. TFinder allows the analysis of an unlimited number of genes; the selection of up to five different species (human, mouse, rat, drosophila, zebrafish); the choice and examination of either promoter or terminator gene regions; the configuration of an upstream downstream window of sequence analysis and the search of TFBS in different formats including IUPAC code, a JASPAR ID or a Position Weight Matrix. TFinder, searches for TFBS on the sense and antisense strand but also considers the search with the complementary forms. The software takes care of everything in record time.

## How to install/use

**No installation is required**. You can access it by clicking here [![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://tfinder-ipmc.streamlit.app/)

A beta version of TFinder exists here. [![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://tfinder-jti99236fxoncgfqnhrcrw.streamlit.app/)


## Functions

### Browser compatibility

- Opera GX
- Chrome (also Chromium)
- Safari
- Edge
- Mozilla
- Phone

# **Gene Regulatory Regions Extractor**

### **Description**
Gene Regulatory Regions Extractor is a tool designed to streamline the extraction of regulatory regions (promoters and terminators) for multiple species and transcript variants. The tool supports customization to extract regions based on gene name, transcript ID, or ENTREZ_GENE_ID, while providing flexibility in choosing genome assemblies.

## **Features**

### **Supported Species**
- **Human** 🙋🏼‍♂️  
- **Mouse** 🖱  
- **Rat** 🐀  
- **Drosophila** 🦟  
- **Zebrafish** 🐟  

### **Extraction Capabilities**
- **Flexible Input Options**:  
  - Extract regulatory regions using:
    - **ENTREZ_GENE_ID**
    - **Gene Name**
    - **Transcript ID** (e.g., NM, XM, NR, XR)
- **Transcript Variants**:  
  - Retrieve regions for specific transcript IDs.
  - Option to extract all transcript variants for a single gene simultaneously.

### **Genome Assembly Selection**
Choose between the current or previous genome assemblies for greater flexibility:

| **Species**     | **Current Genome Assembly** | **Previous Genome Assembly** |
|------------------|-----------------------------|--------------------------------|
| **Human**       | GrCh38                      | GrCh37                         |
| **Mouse**       | GRCm39                      | GRCm38                         |
| **Rat**         | GRCr8                       | mRatBN7.2                      |
| **Drosophila**  | BDGP6                       | BDGP5                          |
| **Zebrafish**   | GRCz11                      | GRCz10                         |

### **Customizable Parameters**
- Set **Upstream** and **Downstream** regions relative to:
  - **Transcription Start Site (TSS)**  
  - **Gene End**
- Fine-tune the length of promoter and terminator regions.

### **Advanced Mode**
The **Advanced Mode** allows for:
- Simultaneous extraction of promoter and terminator regions for multiple species.
- Application of different parameters to each extraction task.

Voici une version adaptée pour votre `README.md` :

---

## **Understanding the FASTA Format for Regulatory Regions**

The FASTA format is used for representing nucleotide sequences, such as promoters and terminators. Below is a guide to help users understand the format and integrate their custom sequences.

### **FASTA Example**
```plaintext
>NM_004562 PRKN | Homo sapiens chromosome 6, GRCh38.p14 Primary Assembly NC_000006.12 | Strand: minus | Promoter | TSS (on chromosome): 162727766 | TSS (on sequence): 2000
TATGAATACAGGTTTAGGAAAAAACAGAAAAGAACCCCAACCAGTAAAAAAAAAATTAAAGTATAACATTAAAAAACATCAAAATTGTAAATATTGTGTAGAAGAAAAACTAAATGATTAACCTGAATGG...
```

### **Key Components**
1. **Header Line (`>`):** Contains metadata about the sequence.
   - `NM_004562 PRKN`: Transcript ID and gene name.
   - `Homo sapiens chromosome 6`: Species and chromosome.
   - `GRCh38.p14 Primary Assembly NC_000006.12`: Genome assembly version and chromosome accession.
   - `Strand: minus`: Indicates the strand (`plus` or `minus`).
   - `Promoter`: Specifies the region type (e.g., Promoter or Terminator).
   - `TSS (on chromosome): 162727766`: Transcription start site (TSS) on the genome.
   - `TSS (on sequence): 2000`: TSS relative to the sequence provided.

2. **Sequence Line:** The nucleotide sequence in uppercase letters (`A`, `T`, `G`, `C`).

---

### **Custom FASTA Integration**

#### **Adding Your Own Sequences**
You can upload custom FASTA files with headers that include metadata. The application automatically extracts the following parameters:
- **Gene name**: Derived from the first part of the header.
- **Species**: Detected from predefined species names (*Homo sapiens*, *Mus musculus*, etc.).
- **Region Type**: Identified as `Promoter` or `Terminator`.
- **Strand**: Extracted as `plus` or `minus`.
- **TSS (on chromosome)**: Parsed from the header if present.

#### **Manually Adding Parameters**
If your header lacks sufficient metadata, you can manually include fields in the following format:
```plaintext
>YourGeneID YourGeneName | YourSpecies | Strand: [plus/minus] | [Promoter/Terminator] | TSS (on chromosome): [TSS_position] | TSS (on sequence): [TSS_on_sequence]
```

If fields are missing, default values will be applied:
- **Species**: `"n.d"` (not determined).
- **Region Type**: `"n.d"`.
- **Strand**: `"n.d"`.
- **TSS (on chromosome)**: `"n.d"`.
- **TSS (on sequence)**: `0`.

---

### **How Metadata is Processed**

1. **Species Detection:** Compares species in the header to the supported list:
   - *Homo sapiens*
   - *Mus musculus*
   - *Rattus norvegicus*
   - *Drosophila melanogaster*
   - *Danio rerio*

   If no match is found, it defaults to `"n.d"`.

2. **Region Type Detection:** Searches for `Promoter` or `Terminator`. Defaults to `"n.d"` if absent.

3. **Strand Validation:** Ensures the strand is valid (`plus` or `minus`). Defaults to `"n.d"` otherwise.

4. **TSS Extraction:** Extracts TSS coordinates if present in the header. Defaults to `0` if missing.

---

### **Custom FASTA Example**
```plaintext
>NM_123456 MyGene | Mus musculus | Strand: plus | Terminator | TSS (on chromosome): 123456789 | TSS (on sequence): 2000
```

This format ensures compatibility with the application and proper parameter extraction.

### **Recommendation**
- Include all relevant metadata in the header for precise parameter assignment.
- Use consistent formatting to facilitate parsing by the tool.

--- 

Voici une version complète et adaptée pour votre `README.md` :  

---

## **Individual Motif Finder**

The **Individual Motif Finder** tool allows for the identification of specific DNA motifs across multiple sequences. It offers flexible input formats, customizable parameters, and robust statistical analysis to ensure accurate motif detection.

---

### **Features**

- **Multiple Input Options:**  
  - Accepts multiple DNA sequences in FASTA format.  
  - Supports motif detection using:
    - **IUPAC codes** (e.g., `Y = C/T`).
    - **JASPAR transcription factor IDs** via the [JASPAR API](https://jaspar.genereg.net/api/v1/docs/).  
    - **Custom Position Weight Matrices (PWM)** generated from user-provided sequences or uploaded directly.  

---

### **Advanced Analysis Tools**
- Detects motif occurrences (e.g., transcription factor binding sites, enzyme restriction sites, or custom patterns).
- Calculates distances of detected motifs to **Transcription Start Sites (TSS)** or **Gene Ends** and **on Chromosome**.
- Outputs results with **Relative Scores**, **Scores**, **Adjusted Scores**, and **P-values**.
    - **Relative Score** is computed by comparing the motif score to the maximum and minimum values from a reference PWM matrix.
    - **Score** is calculated as below
    - **Score Adjusted** accounts for the background nucleotide frequencies from the sequence being studied.
    - The **Relative Score Adjusted** normalizes this score using the adjusted scores.
    - The **p-value** is calculated by generating 1,000,000 random sequences of the same length as the element being analyzed, based on the nucleotide frequencies (either fixed or derived from the sequence) and compares the relative score of the motif to the scores from these random sequences.
- **Thresholds** for these calculations are automatically set for accurate analysis, though advanced users can manually adjust certain parameters for customization.
- Generates **interactive graphs** and allows users to:
  - Download results in `.xlsx` format.
  - Share results via email.

---

### **Motif Scoring**

#### **Score and Adjusted Score Calculations**
- The **Score** is calculated using a **log-odds PSSM** (Position-Specific Scoring Matrix) generated from the PWM.  
- **Score (and Relative Score)**: Assuming equal base distribution (0.25/base)
- **Adjusted Score (Score Adj.) (and Relative Score Adj.):** Accounts for the nucleotide frequencies of the analyzed sequence rather than assuming equal base distribution (0.25/base).  
- **Customization:**  
  - The pseudocount for converting PWM to log-odds PSSM is handled automatically using the [Perl TFBS module](https://biopython.org/docs/dev/Tutorial/chapter_motifs.html#:~:text=Choice%20of%20pseudocounts%3A) but can be modified in the advanced settings.  
  - Background nucleotide frequencies (0.25/base by default) can also be adjusted.

#### **Score Formula:**
The score for a k-mer sequence is computed as:  
<p align="center">
  <img src="https://latex.codecogs.com/svg.image?{\color{cyan}\text{Score} = \sum_{i=1}^{L} \log_2 \left( \frac{P_{\text{PWM}}(a_i, i)}{P_{\text{bg}}(a_i)} \right)
}" alt="relscore equation">
</p>

Where:
- (P_(PWM)(a_i, i)): Probability of nucleotide (a_i) at position (i) based on the PWM (Position Weight Matrix).
- (P_(bg)(a_i)): Background nucleotide probability for (a_i).
- (L): Length of the sequence being analyzed.


#### **Relative Score Formula:**
The **Relative Score** normalizes the raw score to a scale of 0 to 1:
<p align="center">
  <img src="https://latex.codecogs.com/svg.image?{\color{cyan}\text{Relative Score} = \frac{\text{Score of the element found} - \text{Minimum score of the reference matrix}}{\text{Maximum score of the reference matrix} - \text{Minimum score of the reference matrix}}
}" alt="relscore equation">
</p>

---

### **P-value Calculation**

#### **P-value Definition:**
To estimate the accuracy of motif detection, a p-value is calculated by generating 1,000,000 random sequences of the same length as the detected motif:  
- The **p-value** is the proportion of random sequences with a relative score greater than or equal to the detected motif's relative score.  

#### **P-value Formula:**
<p align="center">
  <img src="https://latex.codecogs.com/svg.image?{\color{cyan}\text{p-value} = \frac{\text{Number of random k-mers with Relative Score} \geq \text{Relative Score of the element found}}{\text{Total number of random k-mers}}
}" alt="relscore equation">
</p>

#### **Customizable Background Frequencies:**
- Users can select between:  
  - **Fixed nucleotide frequencies** (default: 0.25/base).  
  - **Nucleotide frequencies derived from the analyzed sequence.**  

---

### **Visualization and Output**
- **Interactive Graphs:** View motif occurrences and distances to TSS or Gene End.  
- **Download Options:** Save results as `.xlsx` files for further analysis.  
- **Email Integration:** Export results directly via email.
- 
---

![graph_webui](https://raw.githubusercontent.com/Jumitti/TFinder/main/img/promtermoriginal.png)
![graph_webui](https://raw.githubusercontent.com/Jumitti/TFinder/main/img/bsfMS.png)
![graph_webui](https://raw.githubusercontent.com/Jumitti/TFinder/main/img/Graph%20WebUI.png)

## More

Report an issue/bug 🆘 ➡️ [Click here](https://github.com/Jumitti/TFinder/issues/new/choose)

Banner was generated with [Adobe Firefly](https://firefly.adobe.com/inspire/images)

Artwork made by [Minniti Pauline](https://minnitidesign.fr/)

## Credit & Licence & Citation

Copyright (c) 2023 Minniti Julien.

This software is distributed under an MIT licence. Please consult the LICENSE file for more details.

PREPRINT: [https://www.researchsquare.com/article/rs-3782387/v1](https://www.researchsquare.com/article/rs-3782387/v1)

