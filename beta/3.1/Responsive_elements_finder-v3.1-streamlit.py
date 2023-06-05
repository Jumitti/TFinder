import streamlit as st
import requests
import pandas as pd

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
            result_promoter.append(f">{gene_name} | {species} | {chraccver} | TSS (on chromosome): {chrstart}\n{dna_sequence}\n")
            st.session_state['result_promoter'] = result_promoter
            st.session_state['upstream'] = upstream

        return result_promoter

    except Exception as e:
        raise Exception(f"Error retrieving gene information: {str(e)}")

# Streamlit app
st.title('Responsive Elements Finder')

# Promoter Finder
st.header('Promoter Finder')

# Gene ID
gene_id_entry = st.text_area("Gene ID:", value="PRKN\n5071")

# Species
species_combobox = st.selectbox("Species:", ["Human", "Mouse", "Rat"], index=0)

# Upstream
upstream_entry = st.text_input("Upstream:", value="2000")

# Downstream
downstream_entry = st.text_input("Downstream:", value="500")

# Run Promoter Finder
if st.button("Find promoter (~5sec/gene)"):
    with st.spinner("Finding promoters..."):
        gene_ids = gene_id_entry.strip().split("\n")
        upstream = int(upstream_entry)
        downstream = int(downstream_entry)
        try:
            result_promoter = find_promoters(gene_ids, species_combobox, upstream, downstream)
            st.success("Promoters extraction complete!")
        except Exception as e:
            st.error(f"Error finding promoters: {str(e)}")

# Promoter output state
if 'result_promoter' not in st.session_state:
    result_promoter = st.text_area("Promoter:", value="")
else:
    result_promoter_text = "\n".join(st.session_state['result_promoter'])
    result_promoter = st.text_area("Promoter:", value=result_promoter_text)
    st.text("Copy: CTRL+A CTRL+C")

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
        header = ["Position", "Position (TIS)", "Sequence", "% Homology", "Ref seq", "Prom."]
        table.insert(0, header)
    else:
        no_consensus = "No consensus sequence found with the specified threshold."
    return table


# Responsive Elements Finder
st.header('Responsive Elements Finder')

# RE entry
entry_sequence = st.text_input("Responsive element (IUPAC authorized):", value="RRRCWWGYYY")

# TSS entry
if 'upstream' not in locals():
    entry_tis = st.text_input("TSS:", value="0")
else:
    entry_tis_value = st.session_state['upstream']
    entry_tis = st.text_input("TSS:", value=entry_tis_value)

# Threshold
threshold_entry = st.text_input("Threshold (%)", value="80")

# Run Responsive Elements finder
if st.button("Find responsive elements"):
    with st.spinner("Finding responsive elements..."):
        sequence_consensus_input = entry_sequence
        tis_value = int(entry_tis_value)
        threshold = float(threshold_entry)
        try:
            table = find_sequence_consensus(sequence_consensus_input, threshold, tis_value, result_promoter)
            st.success("Finding responsive elements done")
        except Exception as e:
            st.error(f"Error finding responsive elements: {str(e)}")

# RE output
if 'table' in locals():
    for row in table:
        st.write("|".join(str(cell).ljust(15) for cell in row))
else:
    st.text("")

