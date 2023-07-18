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

# Streamlit app
st.title('Responsive Elements Finder 🧬🔎')

#Disposition

# Promoter Finder
st.header('🧬 Promoter and Terminator Extractor')

# Gene ID
gene_id_entry = st.text_area("🔸 :red[**Step 1.1**] Gene ID:", value="PRKN\n5071")

# Species
species_combobox = st.selectbox("🔸 :red[**Step 1.2**] Species:", ["Human", "Mouse", "Rat", "Drosophila", "Zebrafish"], index=0)

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
                result_promoter = find_promoters(gene_ids, species_combobox, upstream, downstream)
                st.success("Promoters extraction complete!")
                if balloons:
                    st.balloons()
            except Exception as e:
                st.error(f"Error finding promoters: {str(e)}")
else:
    if st.button("🧬 :red[**Step 1.5**] Extract terminator (~5sec/gene)"):
        with st.spinner("Finding terminators..."):
            gene_ids = gene_id_entry.strip().split("\n")
            upstream = int(upstream_entry)
            downstream = int(downstream_entry)
            try:
                result_promoter = find_promoters(gene_ids, species_combobox, upstream, downstream)
                st.success("Terminators extraction complete!")
                if balloons:
                    st.balloons()
            except Exception as e:
                st.error(f"Error finding terminators: {str(e)}")

# Promoter output state
if prom_term == 'Promoter':
    if 'result_promoter' not in st.session_state:
        result_promoter = st.text_area("🔸 :red[**Step 1.6**] Promoter:", value="")
    else:
        result_promoter_text = "\n".join(st.session_state['result_promoter'])
        result_promoter = st.text_area("🔸 :red[**Step 1.6**] Promoter:", value=result_promoter_text, help='Copy: Click in sequence, CTRL+A, CTRL+C')
else:
    if 'result_promoter' not in st.session_state:
        result_promoter = st.text_area("🔸 :red[**Step 1.6**] Terminator:", value="")
    else:
        result_promoter_text = "\n".join(st.session_state['result_promoter'])
        result_promoter = st.text_area("🔸 :red[**Step 1.6**] Terminator:", value=result_promoter_text, help='Copy: Click in sequence, CTRL+A, CTRL+C')