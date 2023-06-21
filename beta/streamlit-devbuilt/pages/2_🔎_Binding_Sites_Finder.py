import streamlit as st
import requests
import pandas as pd
import altair as alt
import math
import pickle

st.set_page_config(layout="wide")

# Credit Eastereggs
st.sidebar.markdown("[Github](https://github.com/Jumitti/Responsive-Elements-Finder) by Minniti Julien")

balloons = st.sidebar.checkbox('Balloons')

st.sidebar.title("Help")
with st.sidebar.expander("Video tutorials"):
    st.write("How to extract promoter and find responsive elements")
    st.video('https://www.youtube.com/watch?v=lknbKbZCXuo')
    st.write("How to use FASTA sequences and find responsive elements")
    st.video('https://www.youtube.com/watch?v=QelVLLuNJqs')
    st.write("How to use JASPAR option")
    st.video('https://www.youtube.com/watch?v=DH8PBVqa860')
    
with st.sidebar.expander("Promoter & Terminator Extractor"):
    st.subheader("Gene ID:")
    st.write("ENTREZ_GENE_ID of NCBI and gene names are allowed.")
    st.write("There is no limit to the number of gene names/ENTREZ_GENE_ID. Add them with a line break (like those displayed by default). You can mix ENTREZ_GENE_ID and gene names as long as they are of the same species.")
    st.subheader("Species:")
    st.write("Human, mouse, rat, drosophila and zebrafish are allowed.")
    st.write("If you use several ENTREZ_GENE_ID/gene names, make sure you select the correct species.")
    st.subheader("Upstream/Downstream:")
    st.write("Distance to Transcription Start Site (TSS) in bp.")
    st.image("https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/whatisagene.png")
    st.subheader("Promoter & Terminator:")
    st.write('Use "Find promoter/extractor" button or paste your sequences. FASTA format allowed and required for multiple sequences.')
    st.write('FASTA format: All sequences must have the TSS at the same distance, otherwise you assume the inconsistency of the positions of found sequences')
    
with st.sidebar.expander("Binding Sites Finder"):
    st.subheader("Responsive element:")
    st.write("To use the JASPAR option, check the box and use the JASPAR_ID of your transcription factor.")
    st.write('If you want to use your responsive element, do not check the JASPAR option.')
    st.write('IUPAC allowed')
    st.image("https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/IUPAC.png")
    st.subheader("Transcription Start Site (TSS):")
    st.write('Distance to Transcription Start Site (TSS) in bp')
    st.write('Note: If you use Step 1 , it will be defined automatically.')
    st.subheader("Threshold:")
    st.write('Eliminates responsive element with homology < threshold or score < threshold')
    st.write('Note for JASPAR option: Score is normalized to the maximum PWM score of the requested transcription factor. The result is displayed as a percentage')
    st.write('Note without JASPAR option: Homology is calculated between the responsive element in the promoter and the responsive element requested. The calculation uses the Hamming distance, counts the number of differences and gives a percentage score homology.')

# Streamlit app
st.title('Responsive Elements Finder 🧬🔎')

# Promoter output state

st.header('🔎 Binding Sites Finder')
result_promoter = st.text_area("🔸 :red[**Step 1.1**] Sequence:", value="Paste sequences here (FASTA required for multiple sequences).")

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

                if best_homology_percentage >= threshold:
                    row = [str(position).ljust(8),
                           str(tis_position).ljust(15),
                           sequence_with_context,
                           "{:.1f}".format(best_homology_percentage).ljust(12),
                           variant,
                           shortened_promoter_name]
                    table.append(row)

    if len(table) > 0:
        table.sort(key=lambda x: float(x[3]), reverse=True)
        header = ["Position", "Relative position", "Sequence", "% Homology", "Ref seq", "Promoter"]
        table.insert(0, header)
    else:
        no_consensus = "No consensus sequence found with the specified threshold."
        
    return table

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
def search_sequence(sequence_consensus_input, threshold, tis_value, result_promoter, matrices):
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
        table2.sort(key=lambda x: float(x[3]), reverse=True)
        header = ["Position", "Position (TSS)", "Sequence", "Score %", "Promoter"]
        table2.insert(0, header)
    else:
        no_consensus = "No consensus sequence found with the specified threshold."
        
    return table2
    
# Responsive Elements Finder

# RE entry
jaspar = st.radio('🔸 :red[**Step 1.2**] Respnsive elements type:', ('Manual sequence','JASPAR_ID'))
if jaspar == 'JASPAR_ID':
    entry_sequence = st.text_input("🔸 :red[**Step 1.3**] JASPAR ID:", value="MA0106.1")
else:
    entry_sequence = st.text_input("🔸 :red[**Step 1.3**] Responsive element (IUPAC authorized, take more time):", value="ATGCN")

# TSS entry
entry_tis = st.number_input("🔸 :red[**Step 1.4**] Transcription Start Site (TSS) at (in bp):", 0, 10000, 0, help="Distance of TSS or gene end from begin of sequences. Same distance is required for multiple sequences. Leave '0' if you don't know")

# Threshold
if jaspar == 'JASPAR_ID':
    threshold_entry = st.slider("🔸 :red[**Step 1.5**] Score threshold (%)", 0, 100 ,90)
else:
    threshold_entry = st.slider("🔸 :red[**Step 1.5**] Homology threshold (%)", 0, 100 ,80)

# Run Responsive Elements finder
if st.button("🔎 :red[**Step 1.6**] Find responsive elements"):
    with st.spinner("Finding responsive elements..."):
        sequence_consensus_input = entry_sequence
        tis_value = int(entry_tis)
        threshold = float(threshold_entry)
        try:
            if jaspar == 'JASPAR_ID':
                matrices = matrix_extraction(sequence_consensus_input)
                table2 = search_sequence(sequence_consensus_input, threshold, tis_value, result_promoter, matrices)
            else:
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
            
            chart = alt.Chart(source).mark_circle().encode(
                x=alt.X('Relative position:Q', axis=alt.Axis(title='Relative position (bp)'), sort='ascending'),
                y=alt.Y('% Homology:Q', axis=alt.Axis(title='Homology %'), scale=alt.Scale(domain=[ystart, ystop])),
                color=color_scale,
                tooltip=['Relative position', '% Homology', 'Sequence', 'Ref seq', 'Promoter']
            ).properties(width=600, height=400)

            st.altair_chart(chart, use_container_width=True)
        else:
            st.error("No consensus sequence found with the specified threshold.")
    else:
        st.text("")