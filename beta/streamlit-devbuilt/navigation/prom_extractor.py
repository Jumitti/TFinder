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
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json&rettype=xml"
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
            else:
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

    #Disposition
    st.subheader('🧬 Promoter and Terminator Extractor')
    colprom1, colprom2 = st.columns([0.8,1.2] , gap="small")
    
    # Promoter Finder
    with colprom1:        
        result_promoter = []

    # Gene ID
        st.markdown("🔸 :orange[**Step 1.1**] Gene ID:", help='NCBI gene name and NCBI gene ID allowed')
        gene_id_entry = st.text_area("🔸 :orange[**Step 1.1**] Gene ID:", value="PRKN\n351", label_visibility='collapsed')
        if st.button('🔎 Check genes avaibility', help='Sometimes genes do not have the same name in all species or do not exist.'):
            with st.spinner("Checking genes avaibility..."):
                gene_list = gene_id_entry.strip().split('\n')
                for gene_input in gene_list:
                    if not gene_input.isdigit():
                        species_list = ['human','mouse','rat','drosophila','zebrafish']
                        for species_test in species_list:
                            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_input}[Gene%20Name]+AND+{species_test}[Organism]&retmode=json&rettype=xml"
                            response = requests.get(url)

                            if response.status_code == 200:
                                response_data = response.json()

                                if response_data['esearchresult']['count'] == '0':
                                    st.warning(f"{gene_input} gene not found for {species_test}", icon="⚠️")
    
    with colprom2:
        tab1, tab2 = st.tabs(['Default','Advance'])
        
        with tab1:
            
        # Species
            st.markdown("🔸 :orange[**Step 1.2**] Select species of gene names:")
            species = st.selectbox("🔸 :orange[**Step 1.2**] Select species of gene names:", ["Human", "Mouse", "Rat", "Drosophila", "Zebrafish"], index=0, label_visibility='collapsed')

        # Upstream/Downstream Promoter
            st.markdown("🔸 :orange[**Step 1.3**] Regulatory region:")
            prom_term = st.radio("🔸 :orange[**Step 1.3**] Regulatory region:",('Promoter', 'Terminator'), label_visibility='collapsed')
            if prom_term == 'Promoter':
                st.markdown("🔸 :orange[**Step 1.4**] Upstream/downstream from the TSS (bp)")
                updown_slide = st.slider("🔸 :orange[**Step 1.4**] Upstream/downstream from the TSS (bp)", -10000, 10000, (-2000, 500), step=100, label_visibility='collapsed')
                st.write("Upstream: ", min(updown_slide), " bp from TSS | Downstream: ", max(updown_slide), " bp from TSS")
                upstream_entry = -min(updown_slide)
                downstream_entry = max(updown_slide)
                st.session_state['upstream_entry'] = upstream_entry
            else:
                st.markdown("🔸 :orange[**Step 1.4**] Upstream/downstream from gene end (bp)")
                updown_slide = st.slider("🔸 :orange[**Step 1.4**] Upstream/downstream from gene end (bp)", -10000, 10000, (-500, 2000), step=100, label_visibility='collapsed')
                st.write("Upstream: ", min(updown_slide), " bp from gene end | Downstream: ", max(updown_slide), " bp from gene end")
                upstream_entry = -min(updown_slide)
                downstream_entry = max(updown_slide)
                st.session_state['upstream_entry'] = upstream_entry
                
        # Run Promoter Finder
            if prom_term == 'Promoter':
                if st.button("🧬 :orange[**Step 1.5**] Extract promoter" ,help='(~5sec/gene)'):
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
                if st.button("🧬 :orange[**Step 1.5**] Extract terminator", help='(~5sec/gene)'):
                    with st.spinner("Finding terminators..."):
                        gene_ids = gene_id_entry.strip().split("\n")
                        upstream = int(upstream_entry)
                        downstream = int(downstream_entry)
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            st.success("Terminators extraction complete!")
                        except Exception as e:
                            st.error(f"Error finding terminators: {str(e)}")
    with colprom2:    
        with tab2:
            
            gene_list = gene_id_entry.strip().split('\n')
        
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
            
            st.markdown('**🔸 :orange[Step 1.2]** Select species for all genes:', help='Checking a box allows you to check all the corresponding boxes for each gene. Warning: if you have manually checked boxes in the table, they will be reset.')
            
            species1, species2, species3, species4, species5 = st.columns(5)
            
            with species1:
                all_human = st.checkbox("Human")
            with species2:    
                all_mouse = st.checkbox("Mouse")
            with species3:        
                all_rat = st.checkbox("Rat")
            with species4:    
                all_droso = st.checkbox("Drosophila")
            with species5:        
                all_zebra = st.checkbox("Zebrafish")
            
            st.markdown('**🔸 :orange[Step 1.2]** Select regions for all genes:', help='Checking a box allows you to check all the corresponding boxes for each gene. Warning: if you have manually checked boxes in the table, they will be reset.')
            
            region1, region2 = st.columns(2)
            
            with region1: 
                all_prom = st.checkbox("Promoter")
            with region2:         
                all_term = st.checkbox("Terminator")
                    
            if all_human:
                data_df["human"] = True
            if all_mouse:
                data_df["mouse"] = True
            if all_rat:
                data_df["rat"] = True
            if all_droso:
                data_df["droso"] = True
            if all_zebra:
                data_df["zebra"] = True
            if all_prom:
                data_df["prom"] = True
            if all_term:
                data_df["term"] = True
                
            st.markdown('**🔸 :orange[Step 1.2]** On demand genes table', help="Check the boxes for which you want to extract a sequence. Pay attention that the gene name is equivalent for each species. The choice of species is not available for gene IDs. Parameterize the table last, if you check the boxes above, it resets the whole table.")
                
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
                disabled=["Gene"],
                hide_index=True,
            ) 
            
            st.markdown("🔸 :red[**Step 1.3**] Upstream/downstream from TSS and gene end (bp)")
            updown_slide = st.slider("🔸 :red[**Step 1.3**] Upstream/downstream from TSS and gene end (bp)", -10000, 10000, (-2000, 2000), step=100, label_visibility='collapsed')
            st.write("Upstream: ", min(updown_slide), " bp from TSS and gene end | Downstream: ", max(updown_slide), " bp from TSS and gene end")
            upstream_entry = -min(updown_slide)
            downstream_entry = max(updown_slide)
            st.session_state['upstream_entry'] = upstream_entry
            
            if st.button("🧬 :red[**Step 1.4**] Extract sequences", help="(~5sec/seq)"):
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
                            upstream = int(upstream_entry)
                            downstream = int(downstream_entry)
                            species = 'human'
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            except Exception as e:
                                st.error(f"Error finding promoters: {str(e)}")
                        if mouse_checked == True and prom_checked == True:
                            prom_term = 'Promoter'
                            gene_ids = gene_name.strip().split('\n')
                            upstream = int(upstream_entry)
                            downstream = int(downstream_entry)
                            species = 'mouse'
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            except Exception as e:
                                st.error(f"Error finding promoters: {str(e)}")
                        if rat_checked == True and prom_checked == True:
                            prom_term = 'Promoter'
                            gene_ids = gene_name.strip().split('\n')
                            upstream = int(upstream_entry)
                            downstream = int(downstream_entry)
                            species = 'rat'
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            except Exception as e:
                                st.error(f"Error finding promoters: {str(e)}")
                        if droso_checked == True and prom_checked == True:
                            prom_term = 'Promoter'
                            gene_ids = gene_name.strip().split('\n')
                            upstream = int(upstream_entry)
                            downstream = int(downstream_entry)
                            species = 'drosophila'
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            except Exception as e:
                                st.error(f"Error finding promoters: {str(e)}")
                        if zebra_checked == True and prom_checked == True:
                            prom_term = 'Promoter'
                            gene_ids = gene_name.strip().split('\n')
                            upstream = int(upstream_entry)
                            downstream = int(downstream_entry)
                            species = 'zebrafish'
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            except Exception as e:
                                st.error(f"Error finding promoters: {str(e)}")
                        if human_checked == True and term_checked == True:
                            prom_term = 'Terminator'
                            gene_ids = gene_name.strip().split('\n')
                            upstream = int(upstream_entry)
                            downstream = int(downstream_entry)
                            species = 'human'
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            except Exception as e:
                                st.error(f"Error finding promoters: {str(e)}")
                        if mouse_checked == True and term_checked == True:
                            prom_term = 'Terminator'
                            gene_ids = gene_name.strip().split('\n')
                            upstream = int(upstream_entry)
                            downstream = int(downstream_entry)
                            species = 'mouse'
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            except Exception as e:
                                st.error(f"Error finding promoters: {str(e)}")
                        if rat_checked == True and term_checked == True:
                            prom_term = 'Terminator'
                            gene_ids = gene_name.strip().split('\n')
                            upstream = int(upstream_entry)
                            downstream = int(downstream_entry)
                            species = 'rat'
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            except Exception as e:
                                st.error(f"Error finding promoters: {str(e)}")
                        if droso_checked == True and term_checked == True:
                            prom_term = 'Terminator'
                            gene_ids = gene_name.strip().split('\n')
                            upstream = int(upstream_entry)
                            downstream = int(downstream_entry)
                            species = 'drosophila'
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            except Exception as e:
                                st.error(f"Error finding promoters: {str(e)}")
                        if zebra_checked == True and term_checked == True:
                            prom_term = 'Terminator'
                            gene_ids = gene_name.strip().split('\n')
                            upstream = int(upstream_entry)
                            downstream = int(downstream_entry)
                            species = 'zebrafish'
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            except Exception as e:
                                st.error(f"Error finding promoters: {str(e)}")

    # Promoter output state
    st.divider()
    if 'result_promoter' not in st.session_state:
        st.markdown("🔸 :orange[**Step 2.1**] Sequences:")
        result_promoter = st.text_area("🔸 :orange[**Step 2.1**] Sequences:", value="If Step 1 not used, paste sequences here (FASTA required for multiple sequences).", label_visibility='collapsed')
    else:
        st.markdown("🔸 :orange[**Step 2.1**] Sequences:", help='Copy: Click in sequence, CTRL+A, CTRL+C')
        result_promoter_text = "\n".join(st.session_state['result_promoter'])
        result_promoter = st.text_area("🔸 :orange[**Step 2.1**] Sequences:", value=result_promoter_text, label_visibility='collapsed') 