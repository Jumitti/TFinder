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
import numpy as np
import json
import logomaker
import random

def aio_page():
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
                    raise Exception(f"An error occurorange while retrieving the DNA sequence: {response.status_code}")
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
                    raise Exception(f"An error occurorange while retrieving the DNA sequence: {response.status_code}")
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
            raise Exception(f"Error retrieving gene information: {str(e)}")

    #Disposition

    col1, col2 = st.columns(2)

    # Promoter Finder
    with col1:
        st.header(':orange[Step 1] Promoter and Terminator Extractor')
        st.info("💡 If you have a FASTA sequence, go to :orange[**Step 2**]")
        
        result_promoter = []

    # Gene ID
        gene_id_entry = st.text_area("🔸 :orange[**Step 1.1**] Gene ID:", value="PRKN\n351")
        
        tab1, tab2 = st.tabs(['Default','Advance'])
        
        with tab1:
            
        # Species
            species = st.selectbox("🔸 :orange[**Step 1.2**] Select species of gene names:", ["Human", "Mouse", "Rat", "Drosophila", "Zebrafish"], index=0)

        # Upstream/Downstream Promoter
            prom_term = st.radio(
                "🔸 :orange[**Step 1.3**] Regulatory region:",
                ('Promoter', 'Terminator'))
            if prom_term == 'Promoter':
                updown_slide = st.slider("🔸 :orange[**Step 1.4**] Upstream/downstream from the TSS (bp)", -10000, 10000, (-2000, 500), step=100)
                st.write("Upstream: ", min(updown_slide), " bp from TSS | Downstream: ", max(updown_slide), " bp from TSS")
                upstream_entry = -min(updown_slide)
                downstream_entry = max(updown_slide)
                st.session_state['upstream_entry'] = upstream_entry
            else:
                updown_slide = st.slider("🔸 :orange[**Step 1.4**] Upstream/downstream from gene end (bp)", -10000, 10000, (-500, 2000), step=100)
                st.write("Upstream: ", min(updown_slide), " bp from gene end | Downstream: ", max(updown_slide), " bp from gene end")
                upstream_entry = -min(updown_slide)
                downstream_entry = max(updown_slide)
                st.session_state['upstream_entry'] = upstream_entry

        # Run Promoter Finder
            if prom_term == 'Promoter':
                if st.button("🧬 :orange[**Step 1.5**] Extract promoter (~5sec/gene)"):
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
                if st.button("🧬 :orange[**Step 1.5**] Extract terminator (~5sec/gene)"):
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
            
            st.write('**Auto complete species**')
            
            species1, species2, species3, species4, species5 = st.columns([1,1,0.9,1.1,1],gap="small")
            
            with species1:
                all_human = st.checkbox("Human")
                if all_human:
                    data_df["human"] = True
                else:
                    data_df["human"] = False
            with species2:    
                all_mouse = st.checkbox("Mouse")
                if all_mouse:
                    data_df["mouse"] = True
                else:
                    data_df["mouse"] = False
            with species3:        
                all_rat = st.checkbox("Rat")
                if all_rat:
                    data_df["rat"] = True
                else:
                    data_df["rat"] = False
            with species4:    
                all_droso = st.checkbox("Drosophila")
                if all_droso:
                    data_df["droso"] = True
                else:
                    data_df["droso"] = False
            with species5:        
                all_zebra = st.checkbox("Zebrafish")
                if all_zebra:
                    data_df["zebra"] = True
                else:
                    data_df["zebra"] = False
            
            st.write('**Auto complete regions**')
            
            region1, region2 = st.columns(2)
            
            with region1: 
                all_prom = st.checkbox("Promoter")
                if all_prom:
                    data_df["prom"] = True
                else:
                    data_df["prom"] = False
            with region2:         
                all_term = st.checkbox("Terminator")
                if all_term:
                    data_df["term"] = True
                else:
                    data_df["term"] = False
                    
            st.write('**On demand genes table**')      
                
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
            
            
            updown_slide = st.slider("🔸 :red[**Step 1.4**] Upstream/downstream from TSS and gene end (bp)", -10000, 10000, (-2000, 2000), step=100)
            st.write("Upstream: ", min(updown_slide), " bp from TSS and gene end | Downstream: ", max(updown_slide), " bp from TSS and gene end")
            upstream_entry = -min(updown_slide)
            downstream_entry = max(updown_slide)
            st.session_state['upstream_entry'] = upstream_entry
            
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
    with col2:
        st.header(':orange[Step 2] Binding Sites Finder')
        if 'result_promoter' not in st.session_state:
            result_promoter = st.text_area("🔸 :orange[**Step 2.1**] Sequences:", value="If Step 1 not used, paste sequences here (FASTA required for multiple sequences).")
        else:
            result_promoter_text = "\n".join(st.session_state['result_promoter'])
            result_promoter = st.text_area("🔸 :orange[**Step 2.1**] Sequences:", value=result_promoter_text, help='Copy: Click in sequence, CTRL+A, CTRL+C')

    # Responsive-Elements-Finder
        
    # Extract JASPAR matrix
    def matrix_extraction(sequence_consensus_input):
        jaspar_id = sequence_consensus_input
        url = f"https://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
        response = requests.get(url)
        if response.status_code == 200:
            response_data = response.json()
            matrix = response_data['pfm']
        else:
            messagebox.showerror("Erreur", f"Erreur lors de la récupération de la matrice de fréquence : {response.status_code}")
            return

        return transform_matrix(matrix)

    # Transform JASPAR matrix
    def transform_matrix(matrix):
        reversed_matrix = {base: list(reversed(scores)) for base, scores in matrix.items()}
        complement_matrix = {
            'A': matrix['T'],
            'C': matrix['G'],
            'G': matrix['C'],
            'T': matrix['A']
        }
        reversed_complement_matrix = {base: list(reversed(scores)) for base, scores in complement_matrix.items()}

        return {
            'Original': matrix,
            'Reversed': reversed_matrix,
            'Complement': complement_matrix,
            'Reversed Complement': reversed_complement_matrix
        }

    # Calculate score with JASPAR
    def calculate_score(sequence, matrix):
        score = 0
        for i, base in enumerate(sequence):
            if base in {'A', 'C', 'G', 'T'}:
                base_score = matrix[base]
                score += base_score[i]
        return score

    # Find with JASPAR and manual matrix
    def search_sequence(threshold, tis_value, result_promoter, matrices):
        global table2
        table2 = []
        
        # Promoter input type
        lines = result_promoter
        promoters = []

        first_line = lines
        if first_line.startswith(("A", "T", "C", "G")):
            shortened_promoter_name = "n.d."
            promoter_region = lines
            promoters.append((shortened_promoter_name, promoter_region))
        else:
            lines = result_promoter.split("\n")
            i = 0
            while i < len(lines):
                line = lines[i]
                if line.startswith(">"):
                    promoter_name = line[1:]
                    shortened_promoter_name = promoter_name[:15] if len(promoter_name) > 15 else promoter_name
                    if "Promoter" in promoter_name:
                        region = "Prom."
                    elif "Terminator" in promoter_name:
                        region = "Term."
                    else:
                        region = "n.d"
                    promoter_region = lines[i + 1]
                    promoters.append((shortened_promoter_name, promoter_region, region))
                    i += 2
                else:
                    i += 1
        
        if calc_pvalue:
            for matrix_name, matrix in matrices.items():
                seq_length = len(matrix['A'])
            
            for shortened_promoter_name, promoter_region, region in promoters:
                length_prom = len(promoter_region)
                        
                def generate_random_sequence(length, probabilities):
                    nucleotides = ['A', 'C', 'G', 'T']
                    sequence = random.choices(nucleotides, probabilities, k=length)
                    return ''.join(sequence)

                # Génération des séquences aléatoires une seule fois
                motif_length = seq_length  # Remplacer par la longueur de votre motif
                num_random_seqs = 500000

                count_a = promoter_region.count('A')
                count_t = promoter_region.count('T')
                count_g = promoter_region.count('G')
                count_c = promoter_region.count('C')

                length_prom = len(promoter_region)
                percentage_a = count_a / length_prom
                percentage_t = count_t / length_prom
                percentage_g = count_g / length_prom
                percentage_c = count_c / length_prom

                probabilities = [percentage_a, percentage_c, percentage_g, percentage_t]

                random_sequences = []
                for _ in range(num_random_seqs):
                    random_sequence = generate_random_sequence(motif_length, probabilities)
                    random_sequences.append(random_sequence)

                # Calcul des scores aléatoires à partir des différentes matrices
                random_scores = {}
        
        for matrix_name, matrix in matrices.items():
            seq_length = len(matrix['A'])

            # Max score per matrix
            max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
            min_score = sum(min(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))

            # REF
            for shortened_promoter_name, promoter_region, region in promoters:
                found_positions = []
                length_prom = len(promoter_region)

                if calc_pvalue :
                    
                    matrix_random_scores = []
                    for random_sequence in random_sequences:
                        random_score = calculate_score(random_sequence, matrix)
                        normalized_random_score = (random_score - min_score) / (max_score - min_score)
                        matrix_random_scores.append(normalized_random_score)

                    random_scores = np.array(matrix_random_scores)

                for i in range(len(promoter_region) - seq_length + 1):
                    seq = promoter_region[i:i + seq_length]
                    score = calculate_score(seq, matrix)
                    normalized_score = (score - min_score)/(max_score - min_score)
                    position = int(i)
                    
                    if calc_pvalue : 
                        p_value = (random_scores >= normalized_score).sum() / num_random_seqs

                        found_positions.append((position, seq, normalized_score, p_value))
                    else :
                        found_positions.append((position, seq, normalized_score))
                    
                # Sort positions in descending order of score percentage
                found_positions.sort(key=lambda x: x[1], reverse=True)

                # Creating a results table
                if len(found_positions) > 0:
                    if calc_pvalue :
                        for position, seq, normalized_score, p_value in found_positions:
                            start_position = max(0, position - 3)
                            end_position = min(len(promoter_region), position + len(seq) + 3)
                            sequence_with_context = promoter_region[start_position:end_position]

                            sequence_parts = []
                            for j in range(start_position, end_position):
                                if j < position or j >= position + len(seq):
                                    sequence_parts.append(sequence_with_context[j - start_position].lower())
                                else:
                                    sequence_parts.append(sequence_with_context[j - start_position].upper())

                            sequence_with_context = ''.join(sequence_parts)
                            tis_position = position - tis_value
                            
                            if normalized_score >= threshold:
                                row = [str(position).ljust(8),
                                       str(tis_position).ljust(15),
                                       sequence_with_context,
                                       "{:.6f}".format(normalized_score).ljust(12), "{:.3e}".format(p_value).ljust(12),
                                       shortened_promoter_name, region]
                                table2.append(row)
                    else:
                        for position, seq, normalized_score in found_positions:
                            start_position = max(0, position - 3)
                            end_position = min(len(promoter_region), position + len(seq) + 3)
                            sequence_with_context = promoter_region[start_position:end_position]

                            sequence_parts = []
                            for j in range(start_position, end_position):
                                if j < position or j >= position + len(seq):
                                    sequence_parts.append(sequence_with_context[j - start_position].lower())
                                else:
                                    sequence_parts.append(sequence_with_context[j - start_position].upper())

                            sequence_with_context = ''.join(sequence_parts)
                            tis_position = position - tis_value
                            
                            if normalized_score >= threshold:
                                row = [str(position).ljust(8),
                                       str(tis_position).ljust(15),
                                       sequence_with_context,
                                       "{:.6f}".format(normalized_score).ljust(12),
                                       shortened_promoter_name, region]
                                table2.append(row)

        if len(table2) > 0:
            if calc_pvalue :
                table2.sort(key=lambda x: float(x[3]), reverse=True)
                header = ["Position", "Relative position", "Sequence", "Relative Score", "p-value", "Gene", "Region"]
                table2.insert(0, header)
            else:
                table2.sort(key=lambda x: float(x[3]), reverse=True)
                header = ["Position", "Relative position", "Sequence", "Relative Score", "Gene", "Region"]
                table2.insert(0, header)
            
        else:
            no_consensus = "No consensus sequence found with the specified threshold."
            
        return table2
        
    # Responsive Elements Finder
    with col2:

    # RE entry
        jaspar = st.radio('🔸 :orange[**Step 2.2**] Responsive elements type:', ('Manual sequence','JASPAR_ID','Matrix'))
        if jaspar == 'JASPAR_ID':
            entry_sequence = st.text_input("🔸 :orange[**Step 2.3**] JASPAR ID:", value="MA0106.1")
            st.image(f"https://jaspar.genereg.net/static/logos/all/svg/{entry_sequence}.svg")
        elif jaspar == 'Matrix':
            matrix_type = st.radio('🔸 :orange[**Step 2.2bis**] Matrix:', ('With FASTA sequences','With PWM'))
            if matrix_type == 'With PWM':
                isUIPAC = True
                matrix_text = st.text_area("🔸 :orange[**Step 2.3**] Matrix:", value="A [ 20.0 0.0 0.0 0.0 0.0 0.0 0.0 100.0 0.0 60.0 20.0 ]\nT [ 60.0 20.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]\nG [ 0.0 20.0 100.0 0.0 0.0 100.0 100.0 0.0 100.0 40.0 0.0 ]\nC [ 20.0 60.0 0.0 100.0 100.0 0.0 0.0 0.0 0.0 0.0 80.0 ]", help="Only PWM generated with our tools are allowed")
            else:
                fasta_text = st.text_area("🔸 :orange[**Step 2.3**] Sequences:", value=">seq1\nCTGCCGGAGGA\n>seq2\nAGGCCGGAGGC\n>seq3\nTCGCCGGAGAC\n>seq4\nCCGCCGGAGCG\n>seq5\nAGGCCGGATCG", help='Put FASTA sequences. Same sequence length required ⚠️')
                isUIPAC = True
                def calculate_pwm(sequences):
                    num_sequences = len(sequences) 
                    sequence_length = len(sequences[0])
                    pwm = np.zeros((4, sequence_length))
                    for i in range(sequence_length):
                        counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
                        for sequence in sequences:
                            nucleotide = sequence[i]
                            if nucleotide in counts:
                                counts[nucleotide] += 1
                        pwm[0, i] = counts['A'] / num_sequences
                        pwm[1, i] = counts['T'] / num_sequences
                        pwm[2, i] = counts['G'] / num_sequences
                        pwm[3, i] = counts['C'] / num_sequences

                    return pwm

                def parse_fasta(fasta_text):
                    sequences = []
                    current_sequence = ""

                    for line in fasta_text.splitlines():
                        if line.startswith(">"):
                            if current_sequence:
                                sequences.append(current_sequence)
                            current_sequence = ""
                        else:
                            current_sequence += line

                    if current_sequence:
                        sequences.append(current_sequence)

                    return sequences
                    
                if fasta_text:
                    sequences = parse_fasta(fasta_text)
                    sequences = [seq.upper() for seq in sequences]

                    if len(sequences) > 0:
                        pwm = calculate_pwm(sequences)
                        bases = ['A', 'T', 'G', 'C']
                        pwm_text = ""
                        for i in range(len(pwm)):
                            base_name = bases[i]
                            base_values = pwm[i]

                            base_str = base_name + " ["
                            for value in base_values:
                                base_str += "\t" + format(value) + "\t" if np.isfinite(value) else "\t" + "NA" + "\t"

                            base_str += "]\n"
                            pwm_text += base_str

                        matrix_text = st.text_area("PWM:", value=pwm_text, help="Select and copy for later use. Don't modify.", key="non_editable_text")

                    else:
                        st.warning("You forget FASTA sequences :)")
                    
                    def create_web_logo(sequences):
                        matrix = logomaker.alignment_to_matrix(sequences)
                        logo = logomaker.Logo(matrix, color_scheme = 'classic')

                        return logo

                    sequences_text = fasta_text
                    sequences = []
                    current_sequence = ""
                    for line in sequences_text.splitlines():
                        line = line.strip()
                        if line.startswith(">"):
                            if current_sequence:
                                sequences.append(current_sequence)
                            current_sequence = ""
                        else:
                            current_sequence += line

                    if current_sequence:
                        sequences.append(current_sequence)

                    if sequences:
                        logo = create_web_logo(sequences)
                        st.pyplot(logo.fig)
              
        else:
            IUPAC = st.text_input("🔸 :orange[**Step 2.3**] Responsive element (IUPAC authorized):", value="ATGCN")
            
            IUPAC_code = ['A','T','G','C','R','Y','M','K','W','S','B','D','H','V','N']
            
            if all(char in IUPAC_code for char in IUPAC):
                isUIPAC = True
                # IUPAC code
                def generate_iupac_variants(sequence):
                    iupac_codes = {
                        "R": ["A", "G"],
                        "Y": ["C", "T"],
                        "M": ["A", "C"],
                        "K": ["G", "T"],
                        "W": ["A", "T"],
                        "S": ["C", "G"],
                        "B": ["C", "G", "T"],
                        "D": ["A", "G", "T"],
                        "H": ["A", "C", "T"],
                        "V": ["A", "C", "G"],
                        "N": ["A", "C", "G", "T"]
                    }

                    sequences = [sequence]
                    for i, base in enumerate(sequence):
                        if base.upper() in iupac_codes:
                            new_sequences = []
                            for seq in sequences:
                                for alternative in iupac_codes[base.upper()]:
                                    new_sequence = seq[:i] + alternative + seq[i + 1:]
                                    new_sequences.append(new_sequence)
                            sequences = new_sequences

                    return sequences
                   
                sequences = generate_iupac_variants(IUPAC)
                fasta_text = ""
                for i, seq in enumerate(sequences):
                    fasta_text += f">seq{i + 1}\n{seq}\n"
                    
                def calculate_pwm(sequences):
                    num_sequences = len(sequences)
                    sequence_length = len(sequences[0])
                    pwm = np.zeros((4, sequence_length))
                    for i in range(sequence_length):
                        counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
                        for sequence in sequences:
                            nucleotide = sequence[i]
                            if nucleotide in counts:
                                counts[nucleotide] += 1
                        pwm[0, i] = counts['A'] / num_sequences
                        pwm[1, i] = counts['T'] / num_sequences
                        pwm[2, i] = counts['G'] / num_sequences
                        pwm[3, i] = counts['C'] / num_sequences

                    return pwm

                def parse_fasta(fasta_text):
                    sequences = []
                    current_sequence = ""

                    for line in fasta_text.splitlines():
                        if line.startswith(">"):
                            if current_sequence:
                                sequences.append(current_sequence)
                            current_sequence = ""
                        else:
                            current_sequence += line

                    if current_sequence:
                        sequences.append(current_sequence)

                    return sequences
                    
                if fasta_text:
                    sequences = parse_fasta(fasta_text)
                    sequences = [seq.upper() for seq in sequences]

                    if len(sequences) > 0:
                        pwm = calculate_pwm(sequences)
                        bases = ['A', 'T', 'G', 'C']
                        pwm_text = ""
                        for i in range(len(pwm)):
                            base_name = bases[i]
                            base_values = pwm[i]

                            base_str = base_name + " ["
                            for value in base_values:
                                base_str += "\t" + format(value) + "\t" if np.isfinite(value) else "\t" + "NA" + "\t"

                            base_str += "]\n"
                            pwm_text += base_str

                        matrix_text = st.text_area("PWM:", value=pwm_text, help="Select and copy for later use. Dont't modify.", key="non_editable_text")

                    else:
                        st.warning("You forget FASTA sequences :)")
                    
                    def create_web_logo(sequences):
                        matrix = logomaker.alignment_to_matrix(sequences)
                        logo = logomaker.Logo(matrix, color_scheme = 'classic')

                        return logo

                    sequences_text = fasta_text
                    sequences = []
                    current_sequence = ""
                    for line in sequences_text.splitlines():
                        line = line.strip()
                        if line.startswith(">"):
                            if current_sequence:
                                sequences.append(current_sequence)
                            current_sequence = ""
                        else:
                            current_sequence += line

                    if current_sequence:
                        sequences.append(current_sequence)

                    if sequences:
                        logo = create_web_logo(sequences)
                        st.pyplot(logo.fig)
            else:
                isUIPAC = False

    # TSS entry
        if prom_term == 'Promoter':
            entry_tis = st.number_input("🔸 :orange[**Step 2.4**] Transcription Start Site (TSS) at (in bp):", -10000, 10000, st.session_state['upstream_entry'], help="Distance of TSS or gene end from begin of sequences. Do not modify if you use Step 1")
        elif prom_term == 'Terminator':
            entry_tis = st.number_input("🔸 :orange[**Step 2.4**] Gene end at (in bp):", -10000, 10000, st.session_state['upstream_entry'], help="Distance of TSS or gene end from begin of sequences. Do not modify if you use Step 1.")
        else:
            entry_tis = st.number_input("🔸 :orange[**Step 2.4**] Transcription Start Site (TSS) and gene end at (in bp):", -10000, 10000, st.session_state['upstream_entry'], help="Distance of TSS and gene end from begin of sequences. Do not modify if you use Step 1")

    # Threshold pvalue
        if jaspar:
            threshold_entry = st.slider("🔸 :orange[**Step 2.5**] Relative Score threshold", 0.0, 1.0 ,0.85, step= 0.05)
        
        calc_pvalue = st.checkbox('_p-value_ (Experimental, take more times)')

    # Run Responsive Elements finder
        if result_promoter.startswith(("A", "T", "G", "C", ">")):
            with st.spinner("Finding responsive elements..."):
                tis_value = int(entry_tis)
                threshold = float(threshold_entry)
                try:
                    if jaspar == 'JASPAR_ID':
                        sequence_consensus_input = entry_sequence
                        matrices = matrix_extraction(sequence_consensus_input)
                        table2 = search_sequence(threshold, tis_value, result_promoter, matrices)
                    else:
                        if isUIPAC == False:
                            st.error("Please use IUPAC code for Responsive Elements")
                        else:
                            matrix_lines = matrix_text.split('\n')
                            matrix = {}
                            for line in matrix_lines:
                                line = line.strip()
                                if line:
                                    key, values = line.split('[', 1)
                                    values = values.replace(']', '').split()
                                    values = [float(value) for value in values]
                                    matrix[key.strip()] = values
                            matrices = transform_matrix(matrix)
                            table2 = search_sequence(threshold, tis_value, result_promoter, matrices)            
                except Exception as e:
                    st.error(f"Error finding responsive elements: {str(e)}")
    
    # RE output
    if jaspar == 'JASPAR_ID':
        if 'table2' in locals():
            if len(table2) > 0:
                st.divider()
                jaspar_id = sequence_consensus_input
                url = f"https://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
                response = requests.get(url)
                response_data = response.json()
                TF_name = response_data['name']
                st.success(f"Finding responsive elements done for {TF_name}")
                outcol1, outcol2 = st.columns(2)
                with outcol1:
                    df = pd.DataFrame(table2[1:], columns=table2[0])
                    st.session_state['df'] = df
                    st.dataframe(df, hide_index=True)

                    csv = df.to_csv(index=False).encode('utf-8')

                    st.download_button("💾 Download (.csv)",csv,"file.csv","text/csv",key='download-csv')
                    
                with outcol2:
                    source = df
                    score_range = source['Relative Score'].astype(float)
                    ystart = score_range.min() - 0.02
                    ystop = score_range.max() + 0.02
                    scale = alt.Scale(scheme='category10')
                    color_scale = alt.Color("Gene:N", scale=scale)
                    
                    if calc_pvalue:
                        chart = alt.Chart(source).mark_circle().encode(
                            x=alt.X('Relative position:Q', axis=alt.Axis(title='Relative position (bp)'), sort='ascending'),
                            y=alt.Y('Relative Score:Q', axis=alt.Axis(title='Relative Score'), scale=alt.Scale(domain=[ystart, ystop])),
                            color=color_scale,
                            tooltip=['Relative position', 'Relative Score', 'p-value', 'Sequence', 'Gene', 'Region']
                        ).properties(width=600, height=400)
                                              
                        st.altair_chart(chart, use_container_width=True)
                    else:
                        chart = alt.Chart(source).mark_circle().encode(
                            x=alt.X('Relative position:Q', axis=alt.Axis(title='Relative position (bp)'), sort='ascending'),
                            y=alt.Y('Relative Score:Q', axis=alt.Axis(title='Relative Score'), scale=alt.Scale(domain=[ystart, ystop])),
                            color=color_scale,
                            tooltip=['Relative position', 'Relative Score', 'Sequence', 'Gene', 'Region']
                        ).properties(width=600, height=400)
                                              
                        st.altair_chart(chart, use_container_width=True)
            else: 
                jaspar_id = sequence_consensus_input
                url = f"https://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
                response = requests.get(url)
                response_data = response.json()
                TF_name = response_data['name']
                st.error(f"No consensus sequence found with the specified threshold for {TF_name}")
                
    else:
        if 'table2' in locals():
            if len(table2) > 0:
                st.divider()
                st.success(f"Finding responsive elements done")
                outcol1, outcol2 = st.columns(2)
                with outcol1:
                    df = pd.DataFrame(table2[1:], columns=table2[0])
                    st.session_state['df'] = df
                    st.dataframe(df)
                    
                    csv = df.to_csv(index=False).encode('utf-8')

                    st.download_button("💾 Download (.csv)",csv,"file.csv","text/csv",key='download-csv')
                    
                with outcol2:
                    source = df
                    score_range = source['Relative Score'].astype(float)
                    ystart = score_range.min() - 0.02
                    ystop = score_range.max() + 0.02
                    scale = alt.Scale(scheme='category10')
                    color_scale = alt.Color("Gene:N", scale=scale)
                    
                    if calc_pvalue:
                        chart = alt.Chart(source).mark_circle().encode(
                            x=alt.X('Relative position:Q', axis=alt.Axis(title='Relative position (bp)'), sort='ascending'),
                            y=alt.Y('Relative Score:Q', axis=alt.Axis(title='Relative Score'), scale=alt.Scale(domain=[ystart, ystop])),
                            color=color_scale,
                            tooltip=['Relative position', 'Relative Score', 'p-value', 'Sequence', 'Gene', 'Region']
                        ).properties(width=600, height=400)
                                              
                        st.altair_chart(chart, use_container_width=True)
                    else:
                        chart = alt.Chart(source).mark_circle().encode(
                            x=alt.X('Relative position:Q', axis=alt.Axis(title='Relative position (bp)'), sort='ascending'),
                            y=alt.Y('Relative Score:Q', axis=alt.Axis(title='Relative Score'), scale=alt.Scale(domain=[ystart, ystop])),
                            color=color_scale,
                            tooltip=['Relative position', 'Relative Score', 'Sequence', 'Gene', 'Region']
                        ).properties(width=600, height=400)
                                              
                        st.altair_chart(chart, use_container_width=True)
            else:
                st.error(f"No consensus sequence found with the specified threshold")