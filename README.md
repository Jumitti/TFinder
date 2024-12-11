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

## Features - FULL DOCUMENTATION [HERE](https://jumitti.notion.site/tfinder?pvs=4)

### Promoter and Terminator Extractor

- Support gene name (for Human, Mouse, Rat, Drosophila, Zebrafish), Gene ID and transcripts (NM_, XM_, NR_, XR_)
- Support Current genome version and previous (see [documentation](https://www.notion.so/jumitti/Promoter-and-Terminator-Extraction-15973c6e8b9c8016ad21d5f070eac965?pvs=4#15973c6e8b9c80efae51fc76c56bbc06))
- Support extraction of all transcripts variants
- Extract promoter and/or terminator of gene by setting upstream and downstream from TSS/gene end
- Extract unlimited number of gene
- Advance mode for settings (ex: extract promoter and terminator at the same time)
- Speed: 12 sequence / min

### Motif Finder

- Support IUPAC individual motif
- Support JASPAR ID
- Support Position Weight Matrix
- Create a PWM from FASTA sequence
- Calcul of Score (Weight) and normalization for better interpretation
- Automatic Threshold function to Score normalized
- (Pseudocount, background nucleotide frequencies... can be set by users)
- p-value dependant from nucleotide frequencies from the sequence or imposed
- Speed: 80kb/sec (without p-value)

## How to install/use

**No installation is required**. You can access it by clicking here [![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://tfinder-ipmc.streamlit.app/)

A beta version of TFinder exists here. [![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://tfinder-jti99236fxoncgfqnhrcrw.streamlit.app/)

## Image

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

