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


def tfinder_api():
    st.divider()
    st.markdown("<h3 style='text-align: center; color: black;'>Documentation for TFinder python package</h1>",
                unsafe_allow_html=True)
    st.markdown('How to use:')

    st.divider()
    st.markdown("<h3 style='text-align: center; color: black;'>DNA Region Extractor</h1>",
                unsafe_allow_html=True)

    with st.expander('Analyse gene availability'):
        st.markdown('Some genes do not have the same name in different species. It can also happen that the gene ID is incorrect.')
        analyse_gene = 'NCBIdna.analyse_gene(gene_id)'
        st.code(analyse_gene)

        st.divider()
        st.markdown('**Parameter**')
        st.markdown('**gene_id** (list):')
        st.markdown('Support only one by one gene_id. See example. Analyse if ID GENE is valid or if NAME GENE exist for Human, Mouse, Rat, Drosophila, Zebrafish')

        st.divider()
        st.markdown('**Return**')
        st.markdown('(list)')
        st.markdown('Order of results:')
        st.markdown('gene_id | ID | Human | Mouse | Rat | Drosophila | Zebrafish')

        st.divider()
        st.markdown('**Example**')
        st.markdown('Code:')
        example = '''gene_ids = ['4843', 'PRKN']
results = []

for gene_id in gene_ids:
    results.append(NCBIdna.analyse_gene(gene_id))

print(results)'''
        st.code(example)
        st.markdown('Result:')
        st.code("['4843', '✅', 'n.d', 'n.d', 'n.d', 'n.d', 'n.d'], ['PRKN', 'n.d', '✅', '✅', '✅', '❌', '✅']")

    with st.expander('Extraction of regulatory gene sequences'):
        st.markdown('Recovery of promoter or terminal sequences of a given gene for 5 different species.')
        find_sequences = 'NCBIdna(gene_id, prom_term, upstream, downstream, species=None).find_sequences()'
        st.code(find_sequences)

        st.divider()
        st.markdown('**Parameters**')
        st.markdown('**gene_id** (list):')
        st.markdown('It should be **gene name** or **gene ID**. Support only one by one gene_id.')
        st.markdown('**prom_term** ("promoter" or "terminator"):')
        st.markdown('Choose between "promoter" or "terminator". It defines origin of extraction (TSS for promoter and gene end for terminator). Chromosomal coordinates are display in results')
        st.markdown('**upstream** (integer):')
        st.markdown('Upstream defines the number of base pairs (bp) before the TSS/gene end ')
        st.markdown('**downstream** (integer):')
        st.markdown('Downstream defines the number of base pairs (bp) after the TSS/gene end"')
        st.markdown('**species** (str):')
        st.markdown('Any species of https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/. **Note:** if gene ID is use, species is useless')

        st.divider()
        st.markdown('**Return**')
        st.markdown('(str)')
        st.markdown('FASTA output with gene name, species, accession reference of chromosome used, '
                    'regulatory region (prom/term), chromosomal coordinate of TSS or gene end, distance of TSS or gene '
                    'end on the sequence and sequence extracted ')

        st.divider()
        st.markdown('**Example**')
        exemple = '''gene_ids = ['4843', 'PRKN']
upstream = 1000
downstream = 500

results = []
for gene_id in gene_ids:
    if gene_id.isdigit():
        prom_term = 'Promoter'
        results.append(NCBIdna(gene_id, prom_term, upstream, downstream).find_sequences())
    else:
        prom_term = 'Terminator'
        species = 'Macaca mulatta'
        results.append(NCBIdna(gene_id, prom_term, upstream, downstream,species).find_sequences())

print(results)'''
        st.code(exemple)
        st.markdown("**Result**")
        result = '''>NOS2 | Homo sapiens | NC_000017.11 | Promoter | Gene end (on chromosome): 27756765 | Gene end (on sequence): 1000
TTTGAGAGGCTGAAGTGGGCAGATCACTTGAGCTTCAGAGTTCGAGACCAGCATGGACAACATGGTGAAACCCAGTCTCTACCAAAAACACAAAAATATTAGCTGGGTGTGGTGGTGCATGCCTGTAGTCCCAGCTACTCAGGAGGCTGAGGTGGGAGGATCGCTTGAGCCTGGGAGGCAGAAGTTGCAATGAGCAGAGATCGTGCCACTCCGCTCCAGTCTTGGTGACAGAATGAGACTCCATCTCAAAAATAAATAAATAAATAAAATAAATGAAATGAAATTATAAGAAATTACCACTTTTTCATGTAAGAAGTGATCATTTCCATTATAAGGGAAGGAATTTAATCCTACCTGCCATTCCACCAAAGCTTACCTAGTGCTAAAGGATGAGGTGTTAGTAAGACCAACATCTCAGAGGCCTCTCTGTGCCAATAGCCTTCCTTCCTTTCCCTTCCAAAAACCTCAAGTGACTAGTTCAGAGGCCTGTCTGGAATAATGGCATCATCTAATATCACTGGCCTTCTGGAACCTGGGCATTTTCCAGTGTGTTCCATACTGTCAATATTCCCCCAGCTTCCTGGACTCCTGTCACAAGCTGGAAAAGTGAGAGGATGGACAGGGATTAACCAGAGAGCTCCCTGCTGAGGAAAAAATCTCCCAGATGCTGAAAGTGAGGCCATGTGGCTTGGCCAAATAAAACCTGGCTCCGTGGTGCCTCTATCTTAGCAGCCACCCTGCTGATGAACTGCCACCTTGGACTTGGGACCAGAAAGAGGTGGGTTGGGTGAAGAGGCACCACACAGAGTGATGTAACAGCAAGATCAGGTCACCCACAGGCCCTGGCAGTCACAGTCATAAATTAGCTAACTGTACACAAGCTGGGGACACTCCCTTTGGAAACCAAAAAAAAAAAAAAAAAAAAGAGACCTTTATGCAAAAACAACTCTCTGGATGGCATGGGGTGAGTATAAATACTTCTTGGCTGCCAGTGTGTTCATAACTTTGTAGCGAGTCGAAAACTGAGGCTCCGGCCGCAGAGAACTCAGCCTCATTCCTGCTTTAAAATCTCTCGGCCACCTTTGATGAGGGGACTGGGCAGTTCTAGACAGTCCCGAAGTTCTCAAGGCACAGGTCTCTTCCTGGTTTGACTGTCCTTACCCCGGGGAGGCAGTGCAGCCAGCTGCAAGGTGAGTTGCCTTCATTTCTGGGGAAGCGGCTGTTTTGAGAGGGTTTTGTTTCTTCCTCTTTGAGAAGGCTCAGAAAATTGTGGGAATTTTCTGCCTACAGAGAGAAGGTGTTGGAAAGTCTAGGTAAAAAAATGCCGACATGAGGTTTGGTCATTTGAACATGGCATCTTGGTCAGATTTCTTTTCTGCAAATATTAGCTGTGGTTGTTACATGACAAGGAAAAACTTTCTGAAAAGTTCAAAATGGGAGTTTGTATCTGCTAAGAATTTTTTTTAACAGAGATTATCTTATCAGTCCTTATAACATGTAA

>PRKN | Macaca mulatta | NC_041757.1 | Terminator | Gene end (on chromosome): 9397206 | Gene end (on sequence): 1000
TGTCTGAAGATCTGCCGACTCCAGCATGGCCCCATGGTGACAACAGACCTGCGACAGGAAGCCCAAAGCTCACATCGAAATGGTAGAGAGATCAAAGTCTCTATCGTAAGGGAAAAAAAGAGAGGTCACAGGCATGAGCCCCGGCAGCCAGTGGCTTGTGTCCACACGGAGTCCAGACCCTGATCAGGCCCTGACTTAGTATCGCTGGCAATCCCACTAAATTGCAGTTCCTTACACTGGCCCGATGCAACAAATCAGGTGGCTCCCTTCTGTCACAGGGACCACACAGTGTTTCCCATCATCCATAGCTTTCTTCCTGATGATGTTTGCATGATTGCGCCTTCCCAGCCTGCATGCTGCATTGGGCTTGCAGTGCCTGAACAAGGTTTGCTCCCAAGAGCTCAGGCACCCTAGGAACCCCTGTTAGACTATTAGACTGTCCAGCTTGGTCTCCTTCCCTTTCTTGGTGGTGGTCTTTTCCCTTTCCAGAATAGAACAGTGATTCTTAAAATAAGTTCGAGCAGGCTGGGCACGGTGGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGTGGGTGGATCACGAGGTCAGGAGTTCAAGACCAGCCTGGCCAAGACAGCGAAACCCCATCTCTACTAAAAATACAAAGTTAGCTGGGCGCGGTGGCAGGCGCCTGTAATCCCAGCTTCTCGGGAGGCTGAGGCAGGAGAATCACTTGAACCCGGGGGGCAGAGGTTGCAGTGAGCCGAGATCACGCCACTGAACTCCAGCCTGGGCAACAGAGTAAGACTCTGTCTCAAAAAAAGCAAAATCAGCATCCACTACACATGAAAACGAATCACAGTATTATTTGCACAGGAAGGGTGTAACAAAATATGAATGTATCAAAAAATAGAAATAAAGTCTTTGCAGAAAAAATCTATTTTCTCTGACATGTGTTGAGATTATCTGACAACTCTAAGATTGTACTTAAATTGTCAATAAAAGCATCAAAAGAGCTTCTGAGTCCTGTCTTTGAAGATACATTTTCAATCGATGTGAAGGGTGAGGGAGGGAATCACAGGCTGTGTTACTCCATGGACCTCACTTTTGGCCTGGAGGAATTGCCTGGTGATGAAATGGCCATTTTCATATCGACTCCAAGCTCCAGATTAACCCTGGCTGTTTCTCTTTGTTAAATAAGAGCAGTTTGTTAAGGCTATCAAAATTTGAAAGGGTTCCAAAAAGCAAACCCTCATTCAAAACAATCTAAAAGTTACTTTCTTGCCTTACTAGCCAACTTTTCTAATAAGTGTTTGAAAGAGAAGTTTGCATTGGAAAAGACTATGTATAAGGATTTTTTCTCTCATGTCATTCTAGAAAATTTCAGTGACTTTGCATTGCTGTCAATGTGGGTGTGTTGACGCTGCATTTTCAAAGAACTGCTTTAGTGGCCACGCATTGTTTATTTCACAATGAGCCCATGATGGAGGCATCTGGAATGGAACATTCCGAGAT
'''
        st.code(result)

        with st.expander('Convert gene name to gene ID'):
            st.markdown('Convert gene name to gene ID')
            convert_gene_to_entrez_id = 'NCBIdna.convert_gene_to_entrez_id(gene_name, species)'
            st.code(convert_gene_to_entrez_id)

            st.divider()
            st.markdown('**Parameters**')
            st.markdown('**gene_name** (list):')
            st.markdown('It should be **gene name** or **gene ID**. Support only one by one gene_id.')
            st.markdown('**species** (str):')
            st.markdown(
                'Any species of https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/.')

            st.divider()
            st.markdown('**Return**')
            st.markdown('(str)')
            st.markdown('Gene ID of the gene name for desire species')

            st.divider()
            st.markdown('**Example**')
            exemple = '''gene_ids = ['APP', 'PRKN']
        species_list = ['Macaca mulatta', 'Human']

        results = []
        for gene_id in gene_ids:
            for species in species_list:
                results.append(NCBIdna.convert_gene_to_entrez_id(gene_id, species))

        print(results)'''
            st.code(exemple)
            st.markdown("**Result**")
            result = """'100427716', '351', '710899', '5071'"""
            st.code(result)
