# Copyright (c) 2023 Minniti Julien

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of TFinder and associated documentation files, to deal
# in TFinder without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of TFinder, and to permit persons to whom TFinder is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of TFinder.

# TFINDER IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH TFINDER OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import streamlit as st
from PIL import Image
import datetime

def home_page():
    st.header('TFinder 🧬🔎')
    colhome1, colhome2 = st.columns ([1.5, 0.5])
    with colhome1:
        st.image('https://raw.githubusercontent.com/Jumitti/TFinder/main/img/Graph%20WebUI.png', caption='Just a test')
    with colhome2:
        st.markdown('**Overview**')
        st.markdown('<div style="text-align: justify;">TFinder is a Python easy-to-use web tool for identification of putative Transcription Factor Binding Sites (TFBS) in a sequence. It allows extracting directly the promoter or terminal regions of a gene via the NCBI API for 5 different species, with no limit on the number of genes. The reference pattern (ex: a TFBS) accepts both IUPAC codes and JASPAR entries. It is also possible to generate and to use a Position Weight Matrix (PWM). Finally, the data are presented in tabular form, along with a graph showing the relevance of the TFBSs found as a function of their relative position on the sequence.</div>', unsafe_allow_html=True)
    st.divider()
    st.markdown('<div style="text-align: justify;"><p style="text-indent: 2em;">Transcription Factors (TFs) are proteins that bind to DNA to regulate gene expression. They specifically recognize a nucleotide sequence called a Transcription Factor Binding Site (TFBS) present in the 5’ end (most of the time proximal and core promoter regions) and sometimes 3’ end of the regulatory sequences of a gene. The identification   and validation of these TFBSs is an empirical discipline in the field of genomics. To validate TF/DNA interaction , Electro-Mobility Shift Assay (EMSA),ChIP, DNase I footprinting/protection and Systematic Evolution of Ligands by Exponential enrichment (SELEX)  are the most used techniques [1]. However, both techniques require approximative localization of the TFBS, which necessitates in silico studies.</p></div>', unsafe_allow_html=True)
    st.markdown("<div style='text-align: justify;'><p style='text-indent: 2em;'>In silico research can be tedious and time-consuming at various stages, especially for novices. Firstly, databases such as NCBI, UCSC and Ensembl provide access to the regulatory sequences of a gene. However, they are not intuitive and user-friendly. Even If we have the desired nucleotide sequence, TF databases such as JASPAR [2] and TRANSFAC [3] have their limitations. It is important to recognize that there is a certain ease of use with these databases. However, it is impossible to use an unreferenced pattern of a new TF does not present in a database. In addition, TRANSFAC is an expensive tool and if it is used occasionally, it is not worth the cost. Other tools such as PROMO [4], TFBIND [5], TFsitescan make it possible to find all the TFs binding to a nucleotide sequence;  never the less, they all use JASPAR and TRANSFAC databases and do not allow use of personal TF and their TFBS. Moreover, these tools are rather archaic and not very user-friendly. There is only MEME that allows research with of your “own” TFBS [6]. MEME has a large tool library but is a niche software suite. FIMO is their most similar tool to TFinder [7].</p></div>", unsafe_allow_html=True)
    st.markdown('')
    st.markdown('<div style="text-align: justify;"><p style="text-indent: 2em;">A big problem has emerged over time. Many people are unfamiliar with "classic" bioinformatics tools such as NCBI, UCSC and Ensembl tools. It is therefore already very complicated for them to understand these sites and retrieve the information they are looking for. We often must search manually (“Ctlr-F”) to find a specific TFBS in a promoter region, especially if it is not referenced in databases. Most of the tools propose a search on the sense and antisense strand but do not consider the complementary forms. Although there are many complementary tools and databases, only one software offers to use your own TFBS. The vast majority of the tools are free but the most successful: TRANSFAC, is expensive.</p></div>', unsafe_allow_html=True)
    st.markdown('<div style="text-align: justify;"><p style="text-indent: 2em;">This is why faced with all these shortcomings, we decided to build a software that meets our needs. First, quickly extract many promoter/terminators regions only with the name of the gene. Then in a second step add the possibility of searching for an unreferenced TFBS and generating a PWM from it. Our Software, TFinder, searches for TFBS on the sense and antisense strand but also considers the search with the complementary forms. In addition, we have authorized the use of TFBS referenced in JASPAR. Our priority for this software, TFinder, is that it is Open source and free.</p></div>', unsafe_allow_html=True)