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
        '<div style="text-align: justify;">TFinder is an easy-to-use Python web portal allowing the identification of Individual Motifs (IM) such as Transcription Factor Binding Sites (TFBS). Using the NCBI API, TFinder extracts either promoter or the gene terminal regions through a simple query based on NCBI gene name or ID. It enables simultaneous analysis across five different species for an unlimited number of genes. TFinder searches for TFBS and IM in different formats, including IUPAC codes and JASPAR entries. Moreover, TFinder also allows the generation and use of a Position Weight Matrix (PWM). Finally, the data can be recovered in a tabular form and a graph showing the relevance of the TFBSs and IMs as well as its location relative to the Transcription Start Site (TSS) or gene end. The results are then sent by email to the user facilitating the subsequent analysis and data analysis sharing.</div>',
        unsafe_allow_html=True)
    st.divider()
    st.markdown(
        '<div style="text-align: justify;"><p style="text-indent: 2em;">A DNA Individual Motif (IM) is a short pattern conserved between species that can be bind by proteins like Transcription Factors (TFs) enabling gene regulation. They specifically recognize a nucleotide IM sequence called Transcription Factor Binding Site (TFBS) either in gene promoter or terminator regions. Searching of TFBSs is an empirical discipline of genomics which is a key step prior to TFBS functional validation either by gel shift assays (EMSA) or by chromatin immunoprecipitation (ChIP). Both techniques allow the examination of the interaction between a TF and DNA sequence (Jayaram, Usvyat and R. Martin 2016). </p></div>',
        unsafe_allow_html=True)
    st.markdown(
        "<div style='text-align: justify;'><p style='text-indent: 2em;'>The in-silico research of IM can be tedious and time-consuming at various stages, especially for academics or biologist not familiar with bioinformatics. Thus, it is first necessary to retrieve the regulatory region sequence (promoter or terminator). This step may be achieved by the utilization of several databases such as NCBI, UCSC or Ensembl, but they are not intuitive and user-friendly. Next, after identifying the regulatory region sequence, one may use TF databases such as JASPAR (Castro-Mondragon et al. 2022) and TRANSFAC (Matys 2006), but they have their limitations. For instance, these platforms do not allow the search of TFBS for an unreferenced TF and may be subject to fees. Other tools such as PROMO (Farre 2003), TFBIND (Tsunoda and Takagi 1999) and TFsitescan allow searching multiple TFBSs in a unique nucleotide sequence; nevertheless, they all use JASPAR and TRANSFAC databases and do not allow a custom IM or un-referenced TFBS. Finally, a few web tools like FiMO, a module of MEME Suite (Grant, Bailey and Noble 2011; Bailey et al. 2015), allow an unreferenced TFBS or an IM. Of note all the above cited tools are rather archaic, not user-friendly and do not allow the retrieval of regulatory regions prior to motif finding.</p></div>",
        unsafe_allow_html=True)
    st.markdown(
        '<div style="text-align: justify;"><p style="text-indent: 2em;">TFinder is an intuitive, easy-to-use, fast analysis open source and free software that allows both the retrieval of sequences and the search of IM in a unique web application. TFinder allows (1) the analysis of an unlimited number of genes in a record time; (2) the selection up to five different species (human, mouse, rat, drosophila, zebrafish); (3) the choice and examination of either promoter and/or terminator gene regions; (4) the search of IM/TFBS in different formats (IUPAC code, JASPAR ID or a Position Weight Matrix (PWM)); (5) the search of IM/TFBS on the sense and antisense strand but also considers with the complementary forms and (6) the export of the resulting analysis by email.</p></div>',
        unsafe_allow_html=True)
