import streamlit as st
import requests
import pandas as pd
import altair as alt
import math
import pickle
import numpy as np

def REF_page():
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
            result_promoter = []
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

                dna_sequence = get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream)

                # Append the result to the result_promoter
                if prom_term == 'Promoter':
                    result_promoter.append(f">{gene_name} | {species} | {chraccver} | TSS (on chromosome): {chrstart}\n{dna_sequence}\n")
                    st.session_state['result_promoter'] = result_promoter
                else:
                    result_promoter.append(f">{gene_name} | {species} | {chraccver} | Gene end (on chromosome): {chrstop}\n{dna_sequence}\n")
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

    # Gene ID
        gene_id_entry = st.text_area("🔸 :orange[**Step 1.1**] Gene ID:", value="PRKN\n5071")

    # Species
        species_combobox = st.selectbox("🔸 :orange[**Step 1.2**] Species:", ["Human", "Mouse", "Rat", "Drosophila", "Zebrafish"], index=0)

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
                        result_promoter = find_promoters(gene_ids, species_combobox, upstream, downstream)
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
                        result_promoter = find_promoters(gene_ids, species_combobox, upstream, downstream)
                        st.success("Terminators extraction complete!")
                    except Exception as e:
                        st.error(f"Error finding terminators: {str(e)}")

    # Promoter output state
    with col2:
        st.header(':orange[Step 2] Binding Sites Finder')
        if prom_term == 'Promoter':
            if 'result_promoter' not in st.session_state:
                result_promoter = st.text_area("🔸 :orange[**Step 2.1**] Promoter:", value="If Step 1 not used, paste sequences here (FASTA requiorange for multiple sequences).")
            else:
                result_promoter_text = "\n".join(st.session_state['result_promoter'])
                result_promoter = st.text_area("🔸 :orange[**Step 2.1**] Promoter:", value=result_promoter_text, help='Copy: Click in sequence, CTRL+A, CTRL+C')
        else:
            if 'result_promoter' not in st.session_state:
                result_promoter = st.text_area("🔸 :orange[**Step 2.1**] Terminator:", value="If Step 1 not used, paste sequences here (FASTA requiorange for multiple sequences).")
            else:
                result_promoter_text = "\n".join(st.session_state['result_promoter'])
                result_promoter = st.text_area("🔸 :orange[**Step 2.1**] Terminator:", value=result_promoter_text, help='Copy: Click in sequence, CTRL+A, CTRL+C')

    # Responsive-Elements-Finder

    # Generation of all responsive elements
    def generate_variants(sequence):
        variants = []

        # Original sequence
        variants.append(sequence)

        # Reverse sequence
        variants.append(sequence[::-1])

        # Complementary sequence
        complement_sequence = "".join(reverse_complement(base) for base in sequence)
        variants.append(complement_sequence)
        complement_mirror_sequence = complement_sequence[::-1]
        variants.append(complement_mirror_sequence)

        return variants

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

    # Responsive Elements Finder (consensus sequence)
    def find_sequence_consensus(sequence_consensus_input, threshold, tis_value, result_promoter):
        global table
        table = []
        
        # Transform with IUPAC code
        sequence_consensus = generate_iupac_variants(sequence_consensus_input)

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
                    shortened_promoter_name = promoter_name[:10] if len(promoter_name) > 10 else promoter_name
                    promoter_region = lines[i + 1]
                    promoters.append((shortened_promoter_name, promoter_region))
                    i += 2
                else:
                    i += 1

        # REF
        for shortened_promoter_name, promoter_region in promoters:

            found_positions = []

            for consensus in sequence_consensus:
                variants = generate_variants(consensus)
                for variant in variants:

                    variant_length = len(variant)

                    for i in range(len(promoter_region) - variant_length + 1):
                        sequence = promoter_region[i:i + variant_length]

                        mismatches = sum(a != b for a, b in zip(sequence, variant))  # Mismatches

                        homology_percentage = (variant_length - mismatches) / variant_length * 100  # % Homology

                        # Find best homology sequence
                        better_homology = False
                        for position, _, _, _, best_homology_percentage in found_positions:
                            if abs(i - position) < 1 and homology_percentage <= best_homology_percentage:
                                better_homology = True
                                break

                        if not better_homology:

                            best_homology_percentage = (variant_length - mismatches) / variant_length * 100  # % Homology

                            found_positions.append((i, sequence, variant, mismatches, best_homology_percentage))

            # Sort positions in descending order of homology percentage
            found_positions.sort(key=lambda x: x[4], reverse=True)

            # Creating a results table
            if len(found_positions) > 0:
                for position, sequence, variant, mismatches, best_homology_percentage in found_positions:
                    start_position = max(0, position - 3)
                    end_position = min(len(promoter_region), position + len(sequence) + 3)
                    sequence_with_context = promoter_region[start_position:end_position]

                    sequence_parts = []
                    for j in range(start_position, end_position):
                        if j < position or j >= position + len(sequence):
                            sequence_parts.append(sequence_with_context[j - start_position].lower())
                        else:
                            sequence_parts.append(sequence_with_context[j - start_position].upper())

                    sequence_with_context = ''.join(sequence_parts)
                    tis_position = position - tis_value
                    
                    p_value = calculate_p_value(len(variant), mismatches, best_homology_percentage)

                    if best_homology_percentage >= threshold:
                        row = [str(position).ljust(8),
                               str(tis_position).ljust(15),
                               sequence_with_context,"{:.2e}".format(p_value),
                               "{:.1f}".format(best_homology_percentage).ljust(12),
                               variant,
                               shortened_promoter_name]
                        table.append(row)

        if len(table) > 0:
            if prom_term == 'Promoter':
                table.sort(key=lambda x: float(x[3]), reverse=False)
                header = ["Position", "Position (TSS)", "Sequence", "p-value","% Homology", "Ref seq", "Promoter"]
                table.insert(0, header)
            else:
                table.sort(key=lambda x: float(x[3]), reverse=False)
                header = ["Position", "Position (Gene end)", "Sequence", "p-value", "% Homology", "Ref seq", "Promoter"]
                table.insert(0, header)
        else:
            no_consensus = "No consensus sequence found with the specified threshold."
            
        return table
    
    #p-value calcul for manual sequence
    def calculate_p_value(sequence_length, mismatches, best_homology_percentage):
        p_value = 0.0

        for i in range(mismatches, sequence_length + 1):
            p = (math.factorial(sequence_length) / (math.factorial(i) * math.factorial(sequence_length - i))) \
                * (best_homology_percentage / 100) ** i * ((100 - best_homology_percentage) / 100) ** (sequence_length - i)
            p_value += p

        return 1 - p_value
        
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

    # Find with JASPAR
    def search_sequence(threshold, tis_value, result_promoter, matrices):
        global table2
        table2 = []
        
        for matrix_name, matrix in matrices.items():
            seq_length = len(matrix['A'])

            # Max score per matrix
            max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))

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
                        shortened_promoter_name = promoter_name[:10] if len(promoter_name) > 10 else promoter_name
                        promoter_region = lines[i + 1]
                        promoters.append((shortened_promoter_name, promoter_region))
                        i += 2
                    else:
                        i += 1

            # REF
            for shortened_promoter_name, promoter_region in promoters:
                found_positions = []
                total_promoter = len(promoters)

                for i in range(len(promoter_region) - seq_length + 1):
                    seq = promoter_region[i:i + seq_length]
                    score = calculate_score(seq, matrix)
                    normalized_score = (score / max_score) * 100
                    position = int(i)

                    found_positions.append((position, seq, normalized_score))

                # Sort positions in descending order of score percentage
                found_positions.sort(key=lambda x: x[1], reverse=True)

                # Creating a results table
                if len(found_positions) > 0:
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
                                   "{:.1f}".format(normalized_score).ljust(12),
                                   shortened_promoter_name]
                            table2.append(row)

        if len(table2) > 0:
            if prom_term == 'Promoter':
                table2.sort(key=lambda x: float(x[3]), reverse=False)
                header = ["Position", "Position (TSS)", "Sequence", "Score %", "Promoter"]
                table2.insert(0, header)
            else:
                table2.sort(key=lambda x: float(x[3]), reverse=False)
                header = ["Position", "Position (Gene end)", "Sequence", "Score %", "Promoter"]
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
        elif jaspar == 'Matrix':
            matrix_str = st.text_area("🔸 :orange[**Step 2.3**] Matrix:", value="A [ 20.0 0.0 0.0 0.0 0.0 0.0 0.0 100.0 0.0 60.0 20.0 ]\nT [ 60.0 20.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]\nG [ 0.0 20.0 100.0 0.0 0.0 100.0 100.0 0.0 100.0 40.0 0.0 ]\nC [ 20.0 60.0 0.0 100.0 100.0 0.0 0.0 0.0 0.0 0.0 80.0 ]")
        else:
            entry_sequence = st.text_input("🔸 :orange[**Step 2.3**] Responsive element (IUPAC authorized, take more time):", value="ATGCN")

    # TSS entry
        if prom_term == 'Promoter':
            entry_tis = st.number_input("🔸 :orange[**Step 2.4**] Transcription Start Site (TSS) at (in bp):", 0, 10000, st.session_state['upstream_entry'], help="Distance of TSS or gene end from begin of sequences. Do not modify if you use Step 1")
        else:
            entry_tis = st.number_input("🔸 :orange[**Step 2.4**] Gene end at (in bp):", 0, 10000, st.session_state['upstream_entry'], help="Distance of TSS or gene end from begin of sequences. Do not modify if you use Step 1.")

    # Threshold
        if jaspar == 'JASPAR_ID':
            threshold_entry = st.slider("🔸 :orange[**Step 2.5**] Score threshold (%)", 0, 100 ,90)
        else:
            threshold_entry = st.slider("🔸 :orange[**Step 2.5**] Homology threshold (%)", 0, 100 ,80)

    # Run Responsive Elements finder
        if st.button("🔎 :orange[**Step 2.6**] Find responsive elements"):
            with st.spinner("Finding responsive elements..."):
                tis_value = int(entry_tis)
                threshold = float(threshold_entry)
                try:
                    if jaspar == 'JASPAR_ID':
                        sequence_consensus_input = entry_sequence
                        matrices = matrix_extraction(sequence_consensus_input)
                        table2 = search_sequence(threshold, tis_value, result_promoter, matrices)
                    elif jaspar == 'Matrix':
                        matrix_lines = matrix_str.strip().split('\n')
                        matrix_data = {}
                        for line in matrix_lines:
                            parts = line.split('[')
                            base = parts[0].strip()
                            scores = [float(score.strip()) for score in parts[1].replace(']', '').split()]
                            matrix_data[base] = scores

                        matrix_array = np.array([matrix_data['A'], matrix_data['C'], matrix_data['G'], matrix_data['T']])

                        # Define the transformation function
                        def transform_matrix(matrix):
                            transformed_matrix = {
                                'A': [matrix[0][4], matrix[0][7], matrix[0][0], matrix[0][1], matrix[0][2], matrix[0][3], matrix[0][5], matrix[0][6], matrix[0][8]],
                                'C': [matrix[1][3], matrix[1][7], matrix[1][8], matrix[1][1], matrix[1][0], matrix[1][2], matrix[1][4], matrix[1][5], matrix[1][6]],
                                'G': [matrix[2][0], matrix[2][1], matrix[2][2], matrix[2][7], matrix[2][6], matrix[2][8], matrix[2][3], matrix[2][4], matrix[2][5]],
                                'T': [matrix[3][0], matrix[3][1], matrix[3][2], matrix[3][5], matrix[3][7], matrix[3][8], matrix[3][3], matrix[3][4], matrix[3][6]]
                            }
                            return transformed_matrix

                        matrices = transform_matrix(matrix_array)
                            
                        table2 = search_sequence_manual(threshold, tis_value, result_promoter, matrices)
                    else:
                        sequence_consensus_input = entry_sequence
                        table = find_sequence_consensus(sequence_consensus_input, threshold, tis_value, result_promoter)                
                except Exception as e:
                    st.error(f"Error finding responsive elements: {str(e)}")

    # RE output
    if jaspar == 'JASPAR_ID':
        if 'table2' in locals():
            if len(table2) > 0:
                jaspar_id = sequence_consensus_input
                url = f"https://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
                response = requests.get(url)
                response_data = response.json()
                TF_name = response_data['name']
                st.success(f"Finding responsive elements done for {TF_name}")
                st.image(f"https://jaspar.genereg.net/static/logos/all/svg/{jaspar_id}.svg")
                df = pd.DataFrame(table2[1:], columns=table2[0])
                st.session_state['df'] = df
                st.dataframe(df)
                st.info("⬆ Copy: select one cell, CTRL+A, CTRL+C, CTRL+V into spreadsheet software.")

                source = df
                score_range = source['Score %'].astype(float)
                ystart = math.floor(score_range.min() - 5)
                ystop = math.floor(score_range.max() + 5)
                scale = alt.Scale(scheme='category10')
                color_scale = alt.Color("Promoter:N", scale=scale)
                
                if prom_term == 'Promoter':
                    chart = alt.Chart(source).mark_circle().encode(
                        x=alt.X('Position (TSS):Q', axis=alt.Axis(title='Relative position to TSS (bp)'), sort='ascending'),
                        y=alt.Y('Score %:Q', axis=alt.Axis(title='Score %'), scale=alt.Scale(domain=[ystart, ystop])),
                        color=color_scale,
                        tooltip=['Position (TSS)', 'Score %', 'Sequence', 'Promoter']
                    ).properties(width=600, height=400)
                else:
                    chart = alt.Chart(source).mark_circle().encode(
                        x=alt.X('Position (Gene end):Q', axis=alt.Axis(title='Relative position to gene end (bp)'), sort='ascending'),
                        y=alt.Y('Score %:Q', axis=alt.Axis(title='Score %'), scale=alt.Scale(domain=[ystart, ystop])),
                        color=color_scale,
                        tooltip=['Position (Gene end)', 'Score %', 'Sequence', 'Promoter']
                    ).properties(width=600, height=400)
                                      
                st.altair_chart(chart, use_container_width=True)
            else: 
                jaspar_id = sequence_consensus_input
                url = f"https://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
                response = requests.get(url)
                response_data = response.json()
                TF_name = response_data['name']
                st.error(f"No consensus sequence found with the specified threshold for {TF_name}")
                st.image(f"https://jaspar.genereg.net/static/logos/all/svg/{jaspar_id}.svg")
        else:
            st.text("")
    if jaspar == 'Matrix':
        if 'table2' in locals():
            if len(table2) > 0:
                st.success(f"Finding responsive elements done")
                df = pd.DataFrame(table2[1:], columns=table2[0])
                st.session_state['df'] = df
                st.dataframe(df)
                st.info("⬆ Copy: select one cell, CTRL+A, CTRL+C, CTRL+V into spreadsheet software.")

                source = df
                score_range = source['Score %'].astype(float)
                ystart = math.floor(score_range.min() - 5)
                ystop = math.floor(score_range.max() + 5)
                scale = alt.Scale(scheme='category10')
                color_scale = alt.Color("Promoter:N", scale=scale)
                
                if prom_term == 'Promoter':
                    chart = alt.Chart(source).mark_circle().encode(
                        x=alt.X('Position (TSS):Q', axis=alt.Axis(title='Relative position to TSS (bp)'), sort='ascending'),
                        y=alt.Y('Score %:Q', axis=alt.Axis(title='Score %'), scale=alt.Scale(domain=[ystart, ystop])),
                        color=color_scale,
                        tooltip=['Position (TSS)', 'Score %', 'Sequence', 'Promoter']
                    ).properties(width=600, height=400)
                else:
                    chart = alt.Chart(source).mark_circle().encode(
                        x=alt.X('Position (Gene end):Q', axis=alt.Axis(title='Relative position to gene end (bp)'), sort='ascending'),
                        y=alt.Y('Score %:Q', axis=alt.Axis(title='Score %'), scale=alt.Scale(domain=[ystart, ystop])),
                        color=color_scale,
                        tooltip=['Position (Gene end)', 'Score %', 'Sequence', 'Promoter']
                    ).properties(width=600, height=400)
                                      
                st.altair_chart(chart, use_container_width=True)
            else:
                st.error(f"No consensus sequence found with the specified threshold")
        else:
            st.text("")        
    else:
        if 'table' in locals():
            if len(table) > 0 :
                st.success("Finding responsive elements done")
                df = pd.DataFrame(table[1:], columns=table[0])
                st.session_state['df'] = df
                st.dataframe(df)
                st.info("⬆ Copy: select one cell, CTRL+A, CTRL+C, CTRL+V into spreadsheet software.")

                source = df
                homology_range = source['% Homology'].astype(float)
                ystart = math.floor(homology_range.min() - 5)
                ystop = math.floor(homology_range.max() + 5)
                scale = alt.Scale(scheme='category10')
                color_scale = alt.Color("Promoter:N", scale=scale)
                
                if prom_term == 'Promoter':
                    chart = alt.Chart(source).mark_circle().encode(
                        x=alt.X('Position (TSS):Q', axis=alt.Axis(title='Relative position to TSS (bp)'), sort='ascending'),
                        y=alt.Y('% Homology:Q', axis=alt.Axis(title='Homology %'), scale=alt.Scale(domain=[ystart, ystop])),
                        color=color_scale,
                        tooltip=['Position (TSS)', '% Homology', 'Sequence', 'Ref seq', 'Promoter']
                    ).properties(width=600, height=400)
                else:
                    chart = alt.Chart(source).mark_circle().encode(
                        x=alt.X('Position (Gene end):Q', axis=alt.Axis(title='Relative position to gene end (bp)'), sort='ascending'),
                        y=alt.Y('% Homology:Q', axis=alt.Axis(title='Homology %'), scale=alt.Scale(domain=[ystart, ystop])),
                        color=color_scale,
                        tooltip=['Position (Gene end)', '% Homology', 'Sequence', 'Ref seq', 'Promoter']
                    ).properties(width=600, height=400)

                st.altair_chart(chart, use_container_width=True)
            else:
                st.error("No consensus sequence found with the specified threshold.")
        else:
            st.text("")