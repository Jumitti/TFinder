import streamlit as st
import requests
import itertools
import pandas as pd

#StreamLit

st.title('Responsive Elements Finder')

# Common features

# Reverse complement
def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
    return complement_sequence

#Promoter Finder

# Convert gene to ENTREZ_GENE_ID
def convert_gene_names_to_entrez_ids(gene_names):
    try:
        entrez_ids = []
        for gene_id in gene_names:
            species = species_combobox
            # Request for ENTREZ_GENE_ID
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_id}[Gene%20Name]+AND+{species}[Organism]&retmode=json&rettype=xml"
            response = requests.get(url)

            if response.status_code == 200:
                response_data = response.json()

                if response_data['esearchresult']['count'] == '0':
                    raise Exception(f"No gene found for name: {gene_name}")

                else:
                    gene_id = response_data['esearchresult']['idlist'][0]
                    entrez_ids.append(gene_id)

            else:
                raise Exception(f"Error during gene search: {response.status_code}")

        return entrez_ids

    except Exception as e:
        raise Exception(f"Error: {str(e)}")

#Gene informations
def get_gene_info(gene_id, species):
    try:
        # Request for gene informations
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json&rettype=xml&species={species}"
        response = requests.get(url)

        if response.status_code == 200:
            response_data = response.json()

            # Extraction of gene informations
            gene_info = response_data['result'][str(gene_id)]

            return gene_info

        else:
            raise Exception(f"Error during extraction of gene information : {response.status_code}")

    except Exception as e:
        raise Exception(f"Error : {str(e)}")

#Promoter finder
def get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream):
    try:
        # Calculating extraction coordinates on the chromosome
        if chrstop > chrstart:
            start = chrstart - upstream
            end = chrstart + downstream
            
            # Request for DNA sequence
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&from={start}&to={end}&rettype=fasta&retmode=text"
            response = requests.get(url)
            
            if response.status_code == 200:
                # Extraction of DNA sequence
                dna_sequence = response.text.split('\n', 1)[1].replace('\n', '')
                
                sequence = dna_sequence
            
            else:
                raise Exception(f"An error occurred while retrieving the DNA sequence : {response.status_code}")
        
        else:
            start = chrstart + upstream
            end = chrstart - downstream
            
            # Request for DNA sequence
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&from={end}&to={start}&rettype=fasta&retmode=text"
            response = requests.get(url)
            
            if response.status_code == 200:
                # Extraction of DNA sequence
                dna_sequence = response.text.split('\n', 1)[1].replace('\n', '')
                
                
                sequence = dna_sequence
                sequence = reverse_complement(sequence)
                
            else:
                raise Exception(f"An error occurred while retrieving the DNA sequence : {response.status_code}")

        return sequence

    except Exception as e:
        raise Exception(f"Error : {str(e)}")

# Display gene and promoter
def get_sequence(result_promoter):
    species = species_combobox
    gene_ids = gene_id_entry.strip().split("\n")
    total_gene_ids = len(gene_ids)
    upstream = int(upstream_entry)
    downstream = int(downstream_entry)
    
    for i, gene_id in enumerate(gene_ids, start=1):
        try:
            number_gene_id = i
            
            # gene name to ENTREZ_GENE_ID
            gene_names = []
            if not gene_id.isdigit():
                gene_id = gene_id.strip("'\'")
                gene_names.append(gene_id)
                gene_id = convert_gene_names_to_entrez_ids(gene_names)
                gene_entrez_id = gene_id.copy()
                
            else:
                gene_id = gene_id.strip("'\"")
                gene_entrez_id = [gene_id]
            
        except Exception as e:
            result_promoter += f"Error retrieving gene information for ID: {gene_id}\nError: {str(e)}\n"
                
        # Gene information retrieval
        for gene_id in gene_entrez_id:
            gene_info = get_gene_info(gene_id, species)
            gene_name = gene_info['name']
            chraccver = gene_info['genomicinfo'][0]['chraccver']
            chrstart = gene_info['genomicinfo'][0]['chrstart']
            chrstop = gene_info['genomicinfo'][0]['chrstop']

            # Promoter retrieval
            dna_sequence = get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream)

            # Append the result to the result_promoter
            result_promoter += f">{gene_name} | {species} | {chraccver} | TSS: {chrstart}\n{dna_sequence}\n\n"

# Promoter Finder STREAMLIT

st.header('Promoter Finder')

# Gene ID
gene_id_entry = st.text_area("Gene ID:", value="PRKN\n5071")

# Species
species_combobox = st.selectbox("Species:", ["Human", "Mouse", "Rat"], index=0)

#Upstream
upstream_entry = st.text_input("Upstream:", value="2000")

#Downstream
downstream_entry = st.text_input("Downstream:", value="500")

#Run Promoter Finder
if st.button("Find promoter (~30sec/gene)"):
    with st.spinner("Finding promoters..."):
        result_promoter = ""
        get_sequence(result_promoter)
    st.success("Promoters extraction complete!")

#Promoter
result_promoter = st.text_area("Promoter:")
st.text("Copy: CTRL+A CTRL+C")

# Responsive ELements Finder
st.header('Responsive Elements Finder')