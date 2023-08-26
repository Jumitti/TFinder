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
    st.divider()
    st.markdown(
        "<h3 style='text-align: center; color: black;'>TFinder: a Python web tool for predicting Transcription Factor Binding Sites</h1>",
        unsafe_allow_html=True)
    with open('./clock.time', 'r') as f:
        last_updated_on = f.readlines()[0]
    st.caption(last_updated_on)
    st.markdown('')
    st.image('https://raw.githubusercontent.com/Jumitti/TFinder/Beta/img/tfinder_schema.png')
    st.markdown('')
    st.markdown('**Overview**')
    st.markdown(
        '<div style="text-align: justify;"> TFinder is a Python easy-to-use web tool for identifying Transcription Factor Binding Sites (TFBS). It allows easy extraction of either the promoter or terminal regions of a gene by simple inquiry of unique NCBI API keys. It allows the simultaneous analysis of five different species of an unlimited number of genes. The tool allows the search of TFBS in different formats, including IUPAC codes and JASPAR entries. Moreover, TFinder also allows the generation and use a Position Weight Matrix (PWM). Finally, the data may be recovered in a tabular form and a graph showing the relevance of the TFBSs as well as its location relative to the transcription start site (TSS) or gene end. </div>',
        unsafe_allow_html=True)
    st.divider()
    st.markdown(
        '<div style="text-align: justify;"><p style="text-indent: 2em;">Transcription factors (TFs) are proteins that bind to DNA to regulate gene expression. They specifically recognize a nucleotide sequence called a transcription factor binding site (TFBS) in the promoter and terminator regions of genes. The search of these TFBSs is an empirical discipline in the field of genomics that concerns a key step before TFBS functional validation by gel shift assays (EMSA) and chromatin immunoprecipitation (ChIP) that allow the examination of the interaction between a TF and DNA. </p></div>',
        unsafe_allow_html=True)
    st.markdown(
        "<div style='text-align: justify;'><p style='text-indent: 2em;'>The in-silico research of TFBS can be tedious and time-consuming at various stages, especially for novices in the discipline. Thus, first, it is necessary to retrieve the promoter or terminator nucleotide sequence of a gene. This step may be achieved by the utilization of several databases such as NCBI, UCSC and Ensembl, but they are not intuitive and user-friendly. Next, after identifying the promoter sequence of interest, one may use TF databases such as JASPAR and TRANSFAC, but they have their limitations. For example, these platforms do not allow the search of TFBS from an unreferenced TF and may be subject to a fee. Other tools such as PROMO, TFBIND, TFsitescan make it possible to find all the TFs binding to a nucleotide sequence; nevertheless, they all use JASPAR and TRANSFAC databases and do not allow use of personal TF and their TFBS. Moreover, these tools are rather archaic and not very user-friendly. There is only MEME that allows research with of your “own” TFBS. MEME has a large tool library but is a niche software suite. FIMO is their most similar tool to TFinder.</p></div>",
        unsafe_allow_html=True)
    st.markdown(
        '<div style="text-align: justify;"><p style="text-indent: 2em;">TFinder is an ultra-intuitive, easy-to-use and fast analysis open source and free tool that allows both the retrieval and search of TFBS in a unique site. TFinder allows the analysis of an unlimited number of genes; the selection of up to five different species (human, mouse, rat, drosophila, zebrafish); the choice and examination of either promoter or terminator gene regions; the configuration of an upstream downstream window of sequence analysis and the search of TFBS in different formats including IUPAC code, a JASPAR ID or a Position Weight Matrix. TFinder, searches for TFBS on the sense and antisense strand but also considers the search with the complementary forms. The software takes care of everything in record time.</p></div>',
        unsafe_allow_html=True)
