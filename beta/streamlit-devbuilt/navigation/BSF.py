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

def BSF_page():
    # Promoter output state

    st.subheader('🔎 Binding Sites Finder')
    result_promoter = st.text_area("🔸 :red[**Step 1.1**] Sequence:", value="Paste sequences here (FASTA required for multiple sequences).")
        
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

    # Find sequence in promoter
    def search_sequence(threshold, tis_value, result_promoter, matrices):
        global table2
        table2 = []
        
        for matrix_name, matrix in matrices.items():
            seq_length = len(matrix['A'])

            # Max score per matrix
            max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
            min_score = sum(min(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))

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

                def generate_random_sequence(length, probabilities):
                    nucleotides = ['A', 'C', 'G', 'T']
                    sequence = random.choices(nucleotides, probabilities, k=length)
                    return ''.join(sequence)

                def calculate_p_value(motif_score, motif_length, num_random_seqs, probabilities):
                    random_scores = []
                    for _ in range(num_random_seqs):
                        random_sequence = generate_random_sequence(motif_length, probabilities)
                        random_score = calculate_score(random_sequence, matrix)  # Remplacer cette fonction par votre calcul de score de motif
                        normalized__random_score = (random_score - min_score)/(max_score - min_score)
                        random_scores.append(random_score)

                    random_scores = np.array(normalized__random_score)
                    p_value = (random_scores >= motif_score).sum() / num_random_seqs

                    return p_value


                for i in range(len(promoter_region) - seq_length + 1):
                    seq = promoter_region[i:i + seq_length]
                    score = calculate_score(seq, matrix)
                    normalized_score = (score - min_score)/(max_score - min_score)
                    position = int(i)

                    found_positions.append((position, seq, normalized_score))
                    
                    motif_score = normalized_score  # Remplacer par votre score de motif
                    motif_length = seq_length  # Remplacer par la longueur de votre motif
                    num_random_seqs = promoter_region/seq_length  # Nombre de séquences aléatoires à générer
                    probabilities = [0.275, 0.225, 0.225, 0.275]  # Probabilités des nucléotides
                    p_value = calculate_p_value(motif_score, motif_length, num_random_seqs, probabilities)

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
                                   "{:.6f}".format(normalized_score).ljust(12), "{:.12e}".format(p_value).ljust(12),
                                   shortened_promoter_name]
                            table2.append(row)

        if len(table2) > 0:
            table2.sort(key=lambda x: float(x[3]), reverse=True)
            header = ["Position", "Relative position", "Sequence", "Score %", "p-value", "Promoter"]
            table2.insert(0, header)
            
        else:
            no_consensus = "No consensus sequence found with the specified threshold."
            
        return table2

    # RE entry
    jaspar = st.radio('🔸 :orange[**Step 2.2**] Responsive elements type:', ('Manual sequence','JASPAR_ID','Matrix'))
    if jaspar == 'JASPAR_ID':
        entry_sequence = st.text_input("🔸 :orange[**Step 2.3**] JASPAR ID:", value="MA0106.1")
        st.image(f"https://jaspar.genereg.net/static/logos/all/svg/{entry_sequence}.svg")
    elif jaspar == 'Matrix':
        matrix_type = st.radio('🔸 :orange[**Step 2.2bis**] Matrix:', ('With FASTA sequences','With PWM'))
        if matrix_type == 'With PWM':
            matrix_text = st.text_area("🔸 :orange[**Step 2.3**] Matrix:", value="A [ 20.0 0.0 0.0 0.0 0.0 0.0 0.0 100.0 0.0 60.0 20.0 ]\nT [ 60.0 20.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]\nG [ 0.0 20.0 100.0 0.0 0.0 100.0 100.0 0.0 100.0 40.0 0.0 ]\nC [ 20.0 60.0 0.0 100.0 100.0 0.0 0.0 0.0 0.0 0.0 80.0 ]", help="Only PWM generated with our tools are allowed")
        else:
            fasta_text = st.text_area("🔸 :orange[**Step 2.3**] Sequences:", value=">seq1\nCTGCCGGAGGA\n>seq2\nAGGCCGGAGGC\n>seq3\nTCGCCGGAGAC\n>seq4\nCCGCCGGAGCG\n>seq5\nAGGCCGGATCG", help='Put FASTA sequences. Same sequence length required ⚠️')
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
        IUPAC = st.text_input("🔸 :orange[**Step 2.3**] Responsive element (IUPAC authorized):", value="ATGCN")
        
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
        

# TSS entry
    entry_tis = st.number_input("🔸 :red[**Step 1.4**] Relative position to TSS or Gene End (in bp):", -10000, 10000, 0, help="Distance of TSS or gene end from begin of sequences. Same distance is required for multiple sequences. Leave '0' if you don't know")

# Threshold
    if jaspar:
            threshold_entry = st.slider("🔸 :orange[**Step 2.5**] Score threshold (%)", 0.0, 1.0 ,0.95, step= 0.05)

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
    else:
        if st.button("🔎 :orange[**Step 2.6**] Find responsive elements"):
            with st.spinner("Finding responsive elements..."):
                tis_value = int(entry_tis)
                threshold = float(threshold_entry)
                try:
                    if jaspar == 'JASPAR_ID':
                        sequence_consensus_input = entry_sequence
                        matrices = matrix_extraction(sequence_consensus_input)
                        table2 = search_sequence(threshold, tis_value, result_promoter, matrices)
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
                jaspar_id = sequence_consensus_input
                url = f"https://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
                response = requests.get(url)
                response_data = response.json()
                TF_name = response_data['name']
                st.success(f"Finding responsive elements done for {TF_name}")
                df = pd.DataFrame(table2[1:], columns=table2[0])
                st.session_state['df'] = df
                st.dataframe(df)
                st.info("⬆ Copy: select one cell, CTRL+A, CTRL+C, CTRL+V into spreadsheet software.")

                source = df
                score_range = source['Score %'].astype(float)
                ystart = score_range.min() - 0.05
                ystop = score_range.max() + 0.05
                scale = alt.Scale(scheme='category10')
                color_scale = alt.Color("Promoter:N", scale=scale)
                
                chart = alt.Chart(source).mark_circle().encode(
                    x=alt.X('Relative position:Q', axis=alt.Axis(title='Relative position (bp)'), sort='ascending'),
                    y=alt.Y('Score %:Q', axis=alt.Axis(title='Score %'), scale=alt.Scale(domain=[ystart, ystop])),
                    color=color_scale,
                    tooltip=['Relative position', 'Score %', 'Sequence', 'Promoter']
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
    else:
        if 'table2' in locals():
            if len(table2) > 0:
                st.success(f"Finding responsive elements done")
                df = pd.DataFrame(table2[1:], columns=table2[0])
                st.session_state['df'] = df
                st.dataframe(df)
                st.info("⬆ Copy: select one cell, CTRL+A, CTRL+C, CTRL+V into spreadsheet software.")

                source = df
                score_range = source['Score %'].astype(float)
                ystart = score_range.min() - 0.05
                ystop = score_range.max() + 0.05
                scale = alt.Scale(scheme='category10')
                color_scale = alt.Color("Promoter:N", scale=scale)
                
                chart = alt.Chart(source).mark_circle().encode(
                    x=alt.X('Relative position:Q', axis=alt.Axis(title='Relative position (bp)'), sort='ascending'),
                    y=alt.Y('Score %:Q', axis=alt.Axis(title='Score %'), scale=alt.Scale(domain=[ystart, ystop])),
                    color=color_scale,
                    tooltip=['Relative position', 'Score %', 'Sequence', 'Promoter']
                ).properties(width=600, height=400)
                                      
                st.altair_chart(chart, use_container_width=True)
            else:
                st.error(f"No consensus sequence found with the specified threshold")
        else:
            st.text("")