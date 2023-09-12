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
        analyse_gene = 'NCBIdna(gene_id).analyse_gene()'
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
    results.append(NCBIdna(gene_id).analyse_gene())

print(results)'''
        st.code(example)
        st.markdown('Result')
        st.code("['4843', '✅', 'n.d', 'n.d', 'n.d', 'n.d', 'n.d'], ['PRKN', 'n.d', '✅', '✅', '✅', '❌', '✅']")
