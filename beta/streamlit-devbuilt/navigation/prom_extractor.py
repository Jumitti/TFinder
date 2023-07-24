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
import requests
import pandas as pd
import altair as alt
import math
import pickle
import time

def prom_extractor_page():
    # Reverse complement
    def reverse_complement(sequence):
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        reverse_sequence = sequence[::-1]
        complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
        return complement_sequence

    # Convert gene to ENTREZ_GENE_ID
    def convert_gene_to_entrez_id(gene, species):
        try:
            if gene.isdigit():
                return gene  # Already an ENTREZ_GENE_ID

            # Request for ENTREZ_GENE_ID
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene}[Gene%20Name]+AND+{species}[Organism]&retmode=json&rettype=xml"
            response = requests.get(url)

            if response.status_code == 200:
                response_data = response.json()

                if response_data['esearchresult']['count'] == '0':
                    raise Exception(f"No gene found for name: {gene}")

                else:
                    gene_id = response_data['esearchresult']['idlist'][0]
                    return gene_id

            else:
                raise Exception(f"Error during gene search: {response.status_code}")

        except Exception as e:
            raise Exception(f"Error: {str(e)}")

    # Get gene information
    def get_gene_info(gene_id, species):
        try:
            # Request gene information
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json&rettype=xml&species={species}"
            response = requests.get(url)

            if response.status_code == 200:
                response_data = response.json()
                gene_info = response_data['result'][str(gene_id)]
                return gene_info

            else:
                raise Exception(f"Error during extraction of gene information: {response.status_code}")

        except Exception as e:
            raise Exception(f"Error: {str(e)}")

    # Get DNA sequence
    def get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream):
        try:
            if prom_term == 'Promoter':
                if chrstop > chrstart:
                    start = chrstart - upstream
                    end = chrstart + downstream
                else:
                    start = chrstart + upstream
                    end = chrstart - downstream

                # Request for DNA sequence
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&from={start}&to={end}&rettype=fasta&retmode=text"
                response = requests.get(url)

                if response.status_code == 200:
                    # Extraction of DNA sequence
                    dna_sequence = response.text.split('\n', 1)[1].replace('\n', '')
                    if chrstop > chrstart:
                        sequence = dna_sequence
                    else:
                        sequence = reverse_complement(dna_sequence)

                    return sequence

                else:
                    raise Exception(f"An error occurred while retrieving the DNA sequence: {response.status_code}")
            elif prom_term == 'Terminator':
                if chrstop > chrstart:
                    start = chrstop - upstream
                    end = chrstop + downstream
                else:
                    start = chrstop + upstream
                    end = chrstop - downstream

                # Request for DNA sequence
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&from={start}&to={end}&rettype=fasta&retmode=text"
                response = requests.get(url)

                if response.status_code == 200:
                    # Extraction of DNA sequence
                    dna_sequence = response.text.split('\n', 1)[1].replace('\n', '')
                    if chrstop > chrstart:
                        sequence = dna_sequence
                    else:
                        sequence = reverse_complement(dna_sequence)

                    return sequence

                else:
                    raise Exception(f"An error occurred while retrieving the DNA sequence: {response.status_code}")
        except Exception as e:
            raise Exception(f"Error: {str(e)}")

    # Promoter Finder
    def find_promoters(gene_ids, species, upstream, downstream):
        try:
            
            for gene_id in gene_ids:
                if gene_id.isdigit():
                    entrez_id = gene_id
                else:
                    entrez_id = convert_gene_to_entrez_id(gene_id, species)

                gene_info = get_gene_info(entrez_id, species)
                gene_name = gene_info['name']
                chraccver = gene_info['genomicinfo'][0]['chraccver']
                chrstart = gene_info['genomicinfo'][0]['chrstart']
                chrstop = gene_info['genomicinfo'][0]['chrstop']
                species_API = gene_info['organism']['scientificname']

                dna_sequence = get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream)

                # Append the result to the result_promoter
                if prom_term == 'Promoter':
                    result_promoter.append(f">{gene_name} | {species_API} | {chraccver} | {prom_term} | TSS (on chromosome): {chrstart}\n{dna_sequence}\n")
                    st.session_state['result_promoter'] = result_promoter
                else:
                    result_promoter.append(f">{gene_name} | {species_API} | {chraccver} | {prom_term} | Gene end (on chromosome): {chrstop}\n{dna_sequence}\n")
                    st.session_state['result_promoter'] = result_promoter

            return result_promoter

        except Exception as e:
            raise Exception(f"Error retrieving gene information: {str(e)} for species {species}")

    # Promoter Finder
    st.subheader('🧬 Promoter and Terminator Extractor')
    result_promoter = []
    # Gene ID
    gene_id_entry = st.text_area("🔸 :red[**Step 1.1**] Gene ID:", value="PRKN\n5071")
    
    tab1, tab2 = st.tabs(['Default','Advance'])
    with tab1:
        # Species
        species = st.selectbox("🔸 :red[**Step 1.2**] Select species of gene names:", ["Human", "Mouse", "Rat", "Drosophila", "Zebrafish"], index=0)

        # Upstream/Downstream Promoter
        prom_term = st.radio(
            "🔸 :red[**Step 1.3**] Regulatory region:",
            ('Promoter', 'Terminator'))
        if prom_term == 'Promoter':
            updown_slide = st.slider("🔸 :red[**Step 1.4**] Upstream/downstream from the TSS (bp)", -10000, 10000, (-2000, 500), step=100)
            st.write("Upstream: ", min(updown_slide), " bp from TSS | Downstream: ", max(updown_slide), " bp from TSS")
            upstream_entry = -min(updown_slide)
            downstream_entry = max(updown_slide)
            st.session_state['upstream_entry'] = upstream_entry
        else:
            updown_slide = st.slider("🔸 :red[**Step 1.4**] Upstream/downstream from gene end (bp)", -10000, 10000, (-500, 2000), step=100)
            st.write("Upstream: ", min(updown_slide), " bp from gene end | Downstream: ", max(updown_slide), " bp from gene end")
            upstream_entry = -min(updown_slide)
            downstream_entry = max(updown_slide)
            st.session_state['upstream_entry'] = upstream_entry

        # Run Promoter Finder
        if prom_term == 'Promoter':
            if st.button("🧬 :red[**Step 1.5**] Extract promoter (~5sec/gene)"):
                with st.spinner("Finding promoters..."):
                    gene_ids = gene_id_entry.strip().split("\n")
                    upstream = int(upstream_entry)
                    downstream = int(downstream_entry)
                    try:
                        result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        st.success("Promoters extraction complete!")
                    except Exception as e:
                        st.error(f"Error finding promoters: {str(e)}")
        else:
            if st.button("🧬 :red[**Step 1.5**] Extract terminator (~5sec/gene)"):
                with st.spinner("Finding terminators..."):
                    gene_ids = gene_id_entry.strip().split("\n")
                    upstream = int(upstream_entry)
                    downstream = int(downstream_entry)
                    try:
                        result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        st.success("Terminators extraction complete!")
                    except Exception as e:
                        st.error(f"Error finding terminators: {str(e)}")

    with tab2:
        gene_list = gene_id_entry.strip().split('\n')
        
        with genetable1:
            for gene_input in gene_list:
                if not gene_input.isdigit():
                    species_list = ['human','mouse','rat','drosophila','zebrafish']
                    for species_test in species_list:
                        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_input}[Gene%20Name]+AND+{species_test}[Organism]&retmode=json&rettype=xml"
                        response = requests.get(url)

                        if response.status_code == 200:
                            response_data = response.json()

                            if response_data['esearchresult']['count'] == '0':
                                st.warning(f"No gene found for : {gene_input} for {species_test}", icon="⚠️")
                
            
        
        data_df = pd.DataFrame(
            {
                "Gene": gene_list,
                "human": [False] * len(gene_list),
                "mouse": [False] * len(gene_list),
                "rat": [False] * len(gene_list),
                "droso": [False] * len(gene_list),
                "zebra": [False] * len(gene_list),
                "prom": [False] * len(gene_list),
                "term": [False] * len(gene_list),
            }
        )
        
        genetable1, genetable2, genetable3 = st.columns([0.7,0.15,0.15])
        
        with genetable2:
        
            st.write('**Auto complete species**')
            
            all_human = st.checkbox("Human")
            if all_human:
                data_df["human"] = True
            else:
                data_df["human"] = False
                
            all_mouse = st.checkbox("Mouse")
            if all_mouse:
                data_df["mouse"] = True
            else:
                data_df["mouse"] = False
                
            all_rat = st.checkbox("Rat")
            if all_rat:
                data_df["rat"] = True
            else:
                data_df["rat"] = False
            
            all_droso = st.checkbox("Drosophila")
            if all_droso:
                data_df["droso"] = True
            else:
                data_df["droso"] = False
                
            all_zebra = st.checkbox("Zebrafish")
            if all_zebra:
                data_df["zebra"] = True
            else:
                data_df["zebra"] = False
        
        with genetable3:
        
            st.write('**Auto complete regions**')
            
            all_prom = st.checkbox("Promoter")
            if all_prom:
                data_df["prom"] = True
            else:
                data_df["prom"] = False
                
            all_term = st.checkbox("Terminator")
            if all_term:
                data_df["term"] = True
            else:
                data_df["term"] = False

        
        with genetable1:
            data_dff = st.data_editor(
                data_df,
                column_config={
                    "human": st.column_config.CheckboxColumn(
                        "Human",
                        default=False,
                    ),
                    "mouse": st.column_config.CheckboxColumn(
                        "Mouse",
                        default=False,
                    ),
                    "rat": st.column_config.CheckboxColumn(
                        "Rat",
                        default=False,
                    ),
                    "droso": st.column_config.CheckboxColumn(
                        "Drosophila",
                        default=False,
                    ),
                    "zebra": st.column_config.CheckboxColumn(
                        "Zebrafish",
                        default=False,
                    ),
                    "prom": st.column_config.CheckboxColumn(
                        "Promoter",
                        default=False,
                    ),
                    "term": st.column_config.CheckboxColumn(
                        "Terminator",
                        default=False,
                    )
                },
                hide_index=True,
            )
        
        promcol, termcol = st.columns(2)
        with promcol:
            updown_slide = st.slider("🔸 :red[**Step 1.4 Promoter**] Upstream/downstream from TSS (bp)", -10000, 10000, (-2000, 500), step=100)
            st.write("Upstream: ", min(updown_slide), " bp from TSS | Downstream: ", max(updown_slide), " bp from TSS")
            upstream_entry_prom = -min(updown_slide)
            downstream_entry_prom = max(updown_slide)
        with termcol:
            updown_slide = st.slider("🔸 :red[**Step 1.4 Terminator**] Upstream/downstream from gene end (bp)", -10000, 10000, (-500, 2000), step=100)
            st.write("Upstream: ", min(updown_slide), " bp from gene end | Downstream: ", max(updown_slide), " bp from gene end")
            upstream_entry_term = -min(updown_slide)
            downstream_entry_term = max(updown_slide)
        
        if st.button("🧬 :red[**Step 1.5**] Extract sequences (~5sec/seq)"):
            with st.spinner("Finding sequences..."):
                for i, gene_info in data_dff.iterrows():
                    gene_name = gene_info["Gene"]
                    human_checked = gene_info["human"]
                    mouse_checked = gene_info["mouse"]
                    rat_checked = gene_info["rat"]
                    droso_checked = gene_info["droso"]
                    zebra_checked = gene_info["zebra"]
                    prom_checked = gene_info["prom"]
                    term_checked = gene_info["term"]
                
                    if human_checked == True and prom_checked == True:
                        prom_term = 'Promoter'
                        gene_ids = gene_name.strip().split('\n')
                        upstream = int(upstream_entry_prom)
                        downstream = int(downstream_entry_prom)
                        species = 'human'
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        except Exception as e:
                            st.error(f"Error finding promoters: {str(e)}")
                    if mouse_checked == True and prom_checked == True:
                        prom_term = 'Promoter'
                        gene_ids = gene_name.strip().split('\n')
                        upstream = int(upstream_entry_prom)
                        downstream = int(downstream_entry_prom)
                        species = 'mouse'
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        except Exception as e:
                            st.error(f"Error finding promoters: {str(e)}")
                    if rat_checked == True and prom_checked == True:
                        prom_term = 'Promoter'
                        gene_ids = gene_name.strip().split('\n')
                        upstream = int(upstream_entry_prom)
                        downstream = int(downstream_entry_prom)
                        species = 'rat'
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        except Exception as e:
                            st.error(f"Error finding promoters: {str(e)}")
                    if droso_checked == True and prom_checked == True:
                        prom_term = 'Promoter'
                        gene_ids = gene_name.strip().split('\n')
                        upstream = int(upstream_entry_prom)
                        downstream = int(downstream_entry_prom)
                        species = 'drosophila'
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        except Exception as e:
                            st.error(f"Error finding promoters: {str(e)}")
                    if zebra_checked == True and prom_checked == True:
                        prom_term = 'Promoter'
                        gene_ids = gene_name.strip().split('\n')
                        upstream = int(upstream_entry_prom)
                        downstream = int(downstream_entry_prom)
                        species = 'zebrafish'
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        except Exception as e:
                            st.error(f"Error finding promoters: {str(e)}")
                    if human_checked == True and term_checked == True:
                        prom_term = 'Terminator'
                        gene_ids = gene_name.strip().split('\n')
                        upstream = int(upstream_entry_term)
                        downstream = int(downstream_entry_term)
                        species = 'human'
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        except Exception as e:
                            st.error(f"Error finding promoters: {str(e)}")
                    if mouse_checked == True and term_checked == True:
                        prom_term = 'Terminator'
                        gene_ids = gene_name.strip().split('\n')
                        upstream = int(upstream_entry_term)
                        downstream = int(downstream_entry_term)
                        species = 'mouse'
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        except Exception as e:
                            st.error(f"Error finding promoters: {str(e)}")
                    if rat_checked == True and term_checked == True:
                        prom_term = 'Terminator'
                        gene_ids = gene_name.strip().split('\n')
                        upstream = int(upstream_entry_term)
                        downstream = int(downstream_entry_term)
                        species = 'rat'
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        except Exception as e:
                            st.error(f"Error finding promoters: {str(e)}")
                    if droso_checked == True and term_checked == True:
                        prom_term = 'Terminator'
                        gene_ids = gene_name.strip().split('\n')
                        upstream = int(upstream_entry_term)
                        downstream = int(downstream_entry_term)
                        species = 'drosophila'
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        except Exception as e:
                            st.error(f"Error finding promoters: {str(e)}")
                    if zebra_checked == True and term_checked == True:
                        prom_term = 'Terminator'
                        gene_ids = gene_name.strip().split('\n')
                        upstream = int(upstream_entry_term)
                        downstream = int(downstream_entry_term)
                        species = 'zebrafish'
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                        except Exception as e:
                            st.error(f"Error finding promoters: {str(e)}")
        
        
    # Promoter output state
    if 'result_promoter' not in st.session_state:
        result_promoter = st.text_area("🔸 :red[**Step 1.6**] Sequences:", value="")
    else:
        result_promoter_text = "\n".join(st.session_state['result_promoter'])
        result_promoter = st.text_area("🔸 :red[**Step 1.6**] Sequences:", value=result_promoter_text, help='Copy: Click in sequence, CTRL+A, CTRL+C')  