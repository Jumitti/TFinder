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

TEST = "TFinder makes it easy to extract gene regulatory regions by simply providing the gene name or its Gene ID " \
       "(<strong>Fig.1 Step 1.1, Fig.2, Fig.3</strong>). We have added an option to check if the gene is accessible for each species (\"Check genes availability\" button Fig.2). Since the ID gene already takes in account the species, it is not necessary to configure the species analyzed. This also implies that you are not limited to the five species proposed by the program. However, if you use the gene name, then the species will be required (<strong>Fig.1 Step 1.1 and 1.2</strong>). TFinder allows mixing of gene name and gene ID. Please select the desired species. You can therefore easily compare the same regulatory region of 2 or more different species with the gene ID for the same gene."
def resource_page():
    st.divider()
    st.markdown("<h3 style='text-align: center; color: black;'>Introduction</h1>", unsafe_allow_html=True)
    st.markdown('')
    st.markdown('<div style="text-align: justify;"><p style="text-indent: 2em;">TFinder is a Python easy-to-use web tool for the identification of putative Transcription Factor Binding Sites (TFBS) in a sequence. It allows the extraction of the promoter or terminal regions of an unlimited number of genes via the NCBI API of up to five different species. The reference pattern (ex: a TFBS) accepts both IUPAC codes and JASPAR entries. It is also possible to generate and to use a Position Weight Matrix (PWM). Finally, the data may be recovered in either a tabular or graphic formats, showing the relevance score of the TFBSs found as a function of their relative position in the sequence. In this document, we will detail each part of TFinder, the methodology used and the resulting advantages and limitations of the software. TFinder is composed of two main modules that are structured in different sub-modules necessary for its functioning. We will not go into the specific details of the underlying code but will explain the principles and processes. We will first describe how to retrieve a nucleotide sequence on NCBI, then how to find a specific pattern in a nucleotide sequence.</p></div>', unsafe_allow_html=True)
    st.divider()
    st.markdown("<h3 style='text-align: center; color: black;'>Regulatory Regions Extractor (promoter/terminator) 🧬</h1>", unsafe_allow_html=True)
    st.markdown('')
    st.markdown(f'<div style="text-align: justify;"><p style="text-indent: 2em;">{TEST}</p></div>', unsafe_allow_html=True)
    st.markdown("<div style='text-align: justify;'><p style='text-indent: 2em;'>Transcription Factors are proteins that bind to DNA to regulate gene expression. They specifically recognize a nucleotide sequence called a Transcription Factor Binding Site present in the 5’ end (most of the time proximal and core promoter regions) and sometimes 3’ end of the regulatory sequences of a gene (terminator regions) (<strong>Fig.4</strong>). TFinder allows the extraction of these 2 regions. NCBI does not allow direct extraction of regions external to genes. You have to do an extraction directly on the chromosome. We therefore need to know on which chromosome is located the gene of interest, where it begins and ends. The API makes it possible to recover this information as well as GeneBank of the gene (see <strong>Fig.5</strong>). With the API we can find the start and end coordinates of the gene on the chromosome. By having these coordinates, it is possible to know the meaning. Chromosomal coordinates are in the 5` to 3` direction. This means that the first nucleotide is 1 and the last 100000 (for example). If a gene starts at chromosomal coordinates 200 and ends at 500, it is in the 5`-> 3` direction. If a gene starts at 500 but ends at 200, it's the other way around. This makes it possible to transform the extracted region always in the 5` to 3` direction of the gene. After obtaining the start and end coordinates, we can choose to extract the upstream/downstream information by setting the corresponding window (<strong>Fig.4</strong>). Thanks to the ACCESION of the chromosome by setting the coordinates calculated with the upstream/downstream button, one can extract the requested region.</p></div>", unsafe_allow_html=True) 
    st.markdown('<div style="text-align: justify;"><p style="text-indent: 2em;">An "Advanced" mode allows you to choose more specifically what you want to extract (<strong>Fig.6</strong>). With this mode, you can select the region and species you want to recover for each gene. You can make a multiple selection. Of course, it is not possible to choose the species when using gene IDs, only for the gene names. However, in all cases, the region can be selected. The extracted sequences are converted to FASTA format (<strong>Fig.7</strong>).</p></div>', unsafe_allow_html=True) 
    st.markdown('<div style="text-align: justify;"><p style="text-indent: 2em;"><strong>Note 1</strong>: the NCBI API does not allow to extract regions external to genes in a simple way. It is necessary to request the information concerning a piece of the chromosome. But the coordinates are dependent on the requested gene. NCBI does not provide coordinates for different transcripts of the same gene. TFinder displays the coordinates of the TSS or the gene end in the FASTA format in order to easier the transcript identification. You can hover your mouse over the NCBI genetic map to access the coordinates on the chromosome and the different transcripts below (<strong>Fig.8</strong>).</p></div>', unsafe_allow_html=True) 
    st.markdown('<div style="text-align: justify;"><p style="text-indent: 2em;"><strong>Note 2</strong>: If you want to use your own FASTAs, you can. However, verify if they contain the TSS or the end of the gene at the same distance from the beginning of their sequence. Otherwise, you assume the inconsistency of the Relative Position. For the rest, TFinder recognizes if it is a promoter or a terminator if the name of the FASTA contains the "promoter" or "terminator" assignments.</p></div>', unsafe_allow_html=True) 
    tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs(["Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Figure 6", "Figure 7", "Figure 8"])
    with tab1:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/promtermoriginal.png?raw=true', caption='Figure 1: Screenshot of TFinder promoter and terminator extraction tool')
    with tab2:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/promtermcheckgene.png?raw=true', caption='Figure 2: Check genes avaibility')
    with tab3:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/NCBI%20gene%20ID.png?raw=true', caption='Figure 3: NCBI gene name (upper red rectangle) and Gene ID (lower red rectangle)')
    with tab4:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/whatisagene.png?raw=true', caption='Figure 4: What is a gene ?')
    with tab6:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/promtermadvance.png?raw=true', caption='Figure 6: Advance mode')
    with tab5:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/GeneBank.png?raw=true', caption='Figure 5: GeneBank of a gene')
    with tab8:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/coordinates.png?raw=true', caption='Figure 8: Chromosic coordinates on a genetic map from NCBI')
    with tab7:
        st.code('>PRKN | Homo sapiens | NC_000006.12 | Promoter | TSS (on chromosome): 162727765\nATGAATACAGGTTTAGGAAAAAACAGAAAAGAACCCCAACCAGTAAAAAAAAAATTAAAGTATAACATTAAAAAACATCAAAATTGTAAATATTGTGTAGAAGAAAAACTAAATGATTAACCTGAATGGTTATGGTATTGCTGATAAATGCATCATCTTGA\n\n>APP | Homo sapiens | NC_000021.9 | Terminaotr | Gene end (on chromosome): 26171127\nACGCCATTCTCCTGCCTCAGCCTCCCCAGTAGCTGGGACTACAGGCGCCCGCCACGACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTGTTGATCTCCTGACCTCGTGATCCGCCCGCCTCAGCCTCCCAA')
        
    st.divider()
    st.markdown("<h3 style='text-align: center; color: black;'>Individual Motif Finder 🔎</h1>", unsafe_allow_html=True)
    st.markdown('')
    st.markdown('<div style="text-align: justify;"><p style="text-indent: 2em;">TFinder allows you to search for specific patterns in the desired sequences. The sequences must be put in <strong>Step 2.1</strong> (<strong>Fig.1</strong>). The FASTA format is authorized for multiple sequences.', unsafe_allow_html=True)
    st.markdown('<div style="text-align: justify;"><p style="text-indent: 2em;">The pattern to search for can be added in <strong>Step 2.2</strong>. A manual pattern can be used (<strong>Fig.1</strong>). TFinder supports the IUPAC format for your specific pattern (<strong>Fig.1</strong> Step 2.3 and <strong>Fig.2</strong>). A Position Weight Matrix (PWM) is generated to be reused by the user if necessary, as well as a weblogo. You can also use a JASPAR_ID (<strong>Fig.3</strong> and <strong>Fig.4</strong>) or use a PWM (<strong>Fig.5</strong>). If you decide to use a PWM, you can create it from sequences in FASTA format of the same length, or use a PWM generated with our tool.</p></div>', unsafe_allow_html=True) 
    st.markdown("<div style='text-align: justify;'><p style='text-indent: 2em;'><strong>Step 2.4</strong> allows setting the TSS/gene end. It corresponds to the distance of the TSS/gene end from the beginning of the sequence where you are searching. It calculates the Relative Position (Rel Position) at the TSS/gene end for better visualization. If you don't know, put 0 or don't match this column in the results. Otherwise, put the Upstream value that you used previously in positive (example: if you have an upstream of -2000, put in this box 2000). The Relative Score Threshold eliminates the patterns found with a Relative Score (Rel Score) lower than the threshold (<strong>Step 2.5</strong>).</p></div>", unsafe_allow_html=True) 
    st.latex(r'''Relative  \space  Score = \frac {Score - Min \space Score \space PWM}{Max \space Score \space PWM - Min \space Score \space PWM}''')
    st.markdown('<div style="text-align: justify;"><p style="text-indent: 2em;">For each k-mer of the PWM’s length (nucleotide sequence of PWM’s length), a score is calculated by summing the corresponding frequencies of each nucleotide at each position. In order to refine the score calculation, the relative score is used, and calculated as described in the above formula. It allows the normalization with respect to the maximum PWM score of the reference TFBS while subtracting the minimum authorized score. The Relative Score represents the similarity between k-mer and the PWM. Thus, the closer the Relative Score is to 1, the more likely the TFBS identified is not a false positive. The p-value calculation is an information provided by TFinder that requires more processing time. In the p-value calculation, a 1000000 k-mer random sequence is generated taking in consideration the percentage of A, T, G and C in the analyzed target sequence and the length of the PWM. For each k-mer randomly generated, the Relative Score is calculated as described above. The p-value represents the probability that a random sequence has a relative score greater than or equal to the Relative Score of a TFBS. Next, the p-value is presented:</p></div>', unsafe_allow_html=True) 
    st.latex(r'''p \space value = \frac {Nb \space Rel \space Score \space random \space kmer \geq Rel \space Score \space TFBS}{Nb \space random \space kmer}''')
    ttab1, ttab2, ttab3, ttab4, ttab5 = st.tabs(["Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5"])
    with ttab1:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/bsfMS.png?raw=true', caption='Figure 1: Screenshot of TFinder Binding Site Finder tools (option Manual Sequence)')
    with ttab2:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/IUPAC.png?raw=true', caption='Figure 2: IUPAC code')
    with ttab3:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/bsfJI.png?raw=true', caption='Figure 3: Screenshot of TFinder Binding Site Finder tools (option JASPAR)')
    with ttab4:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/JASPAR%20ID.png?raw=true', caption='Figure 4: JASPAR_ID')
    with ttab5:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/bsfM.png?raw=true', caption='Figure 5: Screenshot of TFinder Binding Site Finder tools (option PWM)')

    st.divider()
    st.markdown("<h3 style='text-align: center; color: black;'>References 🔎</h1>", unsafe_allow_html=True)
    st.write("Bailey TL, Johnson J, Grant CE et al. The MEME Suite. Nucleic Acids Res 2015;43:W39–49.", unsafe_allow_html=True)
    st.markdown(
        "Castro-Mondragon JA, Riudavets-Puig R, Rauluseviciute I et al. JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles. Nucleic Acids Research 2022;50:D165–73.", unsafe_allow_html=True)
    st.markdown(
        "Farre D. Identification of patterns in biological sequences at the ALGGEN server: PROMO and MALGEN. Nucleic Acids Research 2003;31:3651–3.", unsafe_allow_html=True)
    st.markdown(
        "Grant CE, Bailey TL, Noble WS. FIMO: scanning for occurrences of a given motif. Bioinformatics 2011;27:1017–8.", unsafe_allow_html=True)
    st.markdown(
        "Jayaram N, Usvyat D, R. Martin AC. Evaluating tools for transcription factor binding site prediction. BMC Bioinformatics 2016;17:547.", unsafe_allow_html=True)
    st.markdown(
        "Khan A, Fornes O, Stigliani A et al. JASPAR 2018: update of the open-access database of transcription factor binding profiles and its web framework. Nucleic Acids Research 2018;46:D260–6.", unsafe_allow_html=True)
    st.markdown(
        "Khan A, Mathelier A. JASPAR RESTful API: accessing JASPAR data from any programming language. Wren J (ed.). Bioinformatics 2018;34:1612–4.", unsafe_allow_html=True)
    st.markdown(
        "Matys V. TRANSFAC(R) and its module TRANSCompel(R): transcriptional gene regulation in eukaryotes. Nucleic Acids Research 2006;34:D108–10.", unsafe_allow_html=True)
    st.markdown(
        "Sayers EW, Bolton EE, Brister JR et al. Database resources of the national center for biotechnology information. Nucleic Acids Research 2022;50:D20<br>–6.", unsafe_allow_html=True)
    st.markdown(
        "Tareen A, Kinney JB. Logomaker: beautiful sequence logos in Python. Valencia A (ed.). Bioinformatics 2020;36:2272–4.", unsafe_allow_html=True)
    st.markdown(
        "Touzet H, Varré J-S. Efficient and accurate P-value computation for Position Weight Matrices. Algorithms Mol Biol 2007;2:15.", unsafe_allow_html=True)
    st.markdown(
        "Tsunoda T, Takagi T. Estimating transcription factor bindability on DNA. Bioinformatics 1999;15:622–30.", unsafe_allow_html=True)
    st.markdown(
        "Turatsinze J-V, Thomas-Chollier M, Defrance M et al. Using RSAT to scan genome sequences for transcription factor binding sites and cis-regulatory modules. Nat Protoc 2008;3:1578–88.", unsafe_allow_html=True)
    st.markdown(
        "Wasserman WW, Sandelin A. Applied bioinformatics for the identification of regulatory elements. Nat Rev Genet 2004;5:276–87.", unsafe_allow_html=True)

    