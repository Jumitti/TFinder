# Copyright (c) 2023 Minniti Julien

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import streamlit as st

def resource_page():
    st.header('Introduction')
    st.write("TFinder is a Python easy-to-use web tool for identification of putative Transcription Factor Binding Sites (TFBS) in a sequence. It allows extracting directly the promoter or terminal regions of a gene via the NCBI API for 5 different species, with no limit on the number of genes. The reference pattern (ex: a TFBS) accepts both IUPAC codes and JASPAR entries. It is also possible to generate and to use a Position Weight Matrix (PWM). Finally, the data are presented in tabular form, along with a graph showing the relevance of the TFBSs found as a function of their relative position on the sequence.In this document, we will detail each part of TFinder, the method used and the resulting advantages and limitations. The software has 2 main modules with different sub-modules necessary for its operation. We will not go into the specific details of the underlying code but will explain the principles and processes. We will first look at how to retrieve a nucleotide sequence on NCBI, then how to find a specific pattern in a nucleotide sequence.")
    st.header('Extraction of gene regulatory sequences (promoter/terminator)')
    st.write("TFinder makes it easy to extract gene regulatory regions from its name or its Gene ID. TFinder accepts gene names and Gene IDs (**Fig.1 Step 1.1 and Fig.2**). Since the ID gene is already species dependent, it is not necessary to configure the species. This also implies that you are not limited to the 5 species offered. However, if you use the gene name, then the species will be required (**Fig.1 Step 1.1 and 1.2**). TFinder allows mixing of gene name and gene ID. Please select the desired species. You can therefore easily compare the same regulatory region of 2 or more different species with the gene ID for the same gene")
    st.write("Transcription Factors are proteins that bind to DNA to regulate gene expression. They specifically recognize a nucleotide sequence called a Transcription Factor Binding Site present in the 5’ end (most of the time proximal and core promoter regions) and sometimes 3’ end of the regulatory sequences of a gene (terminator regions) (**Fig.3**). TFinder allows the extraction of these 2 regions. NCBI does not allow direct extraction of regions external to genes. You have to do an extraction directly on the chromosome. We therefore need to know on which chromosome is the gene of interest, where it begins and ends. The API makes it possible to recover this information and it is also present in the GeneBank of the gene. Once the information has been retrieved, the direction of the gene can be determined (if start coordinate < end coordinate then sense strand, otherwise anti-sense strand). This makes it possible to transform the extracted region always in the 5' to 3' direction of the gene. With the start and end coordinates, we can set the window to extract with the upstream/downstream information. Thanks to the ACCESION of the chromosome by setting the coordinates calculated with the upstream/downstream, we can extract the requested region. okook You can set the upstream and the downstream with the sliders. (Fig.1 Step 1.4 and Fig.3). TFinder can recognize the direction of the gene to always present it in the 5' to 3' direction. After extraction, the sequences are in FASTA format. The NCBI API does not allow the possibility to extract a nucleotide sequence external to a gene. It is then necessary to find the chromosomal coordinates of the gene. The coordinates represent the start and end of the gene on the chromosome. The beginning and the end of the gene are accessible in the GeneBank of the gene in the ACCESSION section (Fig.4). We also find the accession number of the corresponding chromosome that we will use later to extract our region. Therefore, we can know the direction of the genes and use the NCBI API to retrieve a chunk of the gene using upstream and downstream coordinates based on the beginning of the gene or the end.")
    st.write("Limit: the NCBI API does not allow to extract regions external to genes in a simple way. We have to make a request for a piece of chromosome. But the coordinates are dependent on the requested gene. NCBI does not have coordinates for different transcripts of the same gene. We decided to display the coordinates of the TSS or the gene end in the FASTA in order to find more easily which transcript it corresponds to. You can hover your mouse over the NCBI genetic map to see approximate coordinates on the chromosome and the different transcripts below (Fig.5)")
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Figure 6"])
    with tab1:
        st.image('https://raw.githubusercontent.com/Jumitti/TFinder/main/img/NCBI%20API.png', caption='Figure 1: Screenshot of TFinder promoter and terminator extraction tool')
    with tab2:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/NCBI%20gene%20ID.png?raw=true', caption='Figure 2: NCBI gene name (upper red rectangle) and Gene ID (lower red rectangle)')
    with tab3:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/whatisagene.png?raw=true', caption='Figure 3: What is a gene ?')
    with tab4:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/GeneBank.png?raw=true', caption='Figure 4: GeneBank of a gene')
    with tab5:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/coordinates.png?raw=true', caption='Figure 5: Chromosic coordinates on a genetic map from NCBI')
    with tab6:
        st.code('>PRKN | Homo sapiens | NC_000006.12 | Promoter | TSS (on chromosome): 162727765\nATGAATACAGGTTTAGGAAAAAACAGAAAAGAACCCCAACCAGTAAAAAAAAAATTAAAGTATAACATTAAAAAACATCAAAATTGTAAATATTGTGTAGAAGAAAAACTAAATGATTAACCTGAATGGTTATGGTATTGCTGATAAATGCATCATCTTGA\n\n>APP | Homo sapiens | NC_000021.9 | Terminaotr | Gene end (on chromosome): 26171127\nACGCCATTCTCCTGCCTCAGCCTCCCCAGTAGCTGGGACTACAGGCGCCCGCCACGACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTGTTGATCTCCTGACCTCGTGATCCGCCCGCCTCAGCCTCCCAA')
    st.header('Transcription Factors Binding Site')