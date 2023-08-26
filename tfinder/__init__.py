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

import time
import requests


def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
    return complement_sequence

# Get gene information
def get_gene_info(entrez_id):
    # Request gene information
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={entrez_id}&retmode=json&rettype=xml"
    response = requests.get(url)

    if response.status_code == 200:
        response_data = response.json()
        gene_info = response_data['result'][str(entrez_id)]
        return gene_info

class NCBI_dna:
    def __init__(self,
                 gene_id,
                 species=None,
                 upstream=None,
                 downstream=None,
                 prom_term=None):
        self.gene_id = gene_id
        self.upstream = int(upstream) if upstream is not None else None
        self.downstream = int(downstream) if downstream is not None else None
        self.prom_term = prom_term if prom_term is not None else None
        self.species = species if species is not None else None

    @staticmethod
    # Analyse if gene is available
    def analyse_gene(gene_id):
        disponibility_list = ['ID', 'Human', 'Mouse', 'Rat', 'Drosophila', 'Zebrafish']
        time.sleep(0.25)
        gene_analyse = [self.gene_id]
        for species_test in disponibility_list:
            if not self.gene_id.isdigit():
                if species_test == 'ID':
                    gene_analyse.append('n.d')
                else:
                    time.sleep(0.5)
                    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={self.gene_id}[Gene%20Name]+AND+{species_test}[Organism]&retmode=json&rettype=xml"
                    response = requests.get(url)

                    if response.status_code == 200:
                        response_data = response.json()

                        if response_data['esearchresult']['count'] != '0':
                            gene_analyse.append("✅")
                        else:
                            gene_analyse.append("❌")

            if self.gene_id.isdigit():
                if species_test != 'ID':
                    gene_analyse.append('n.d')
                else:
                    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={self.gene_id}&retmode=json&rettype=xml"
                    response = requests.get(url)

                    if response.status_code == 200:
                        response_data = response.json()

                        if 'chraccver' in str(response_data):
                            gene_analyse.append("✅")
                        else:
                            gene_analyse.append("❌")

        return gene_analyse

    # Sequence extractor
    def find_sequences(self):
        time.sleep(1)
        if self.gene_id.isdigit():
            entrez_id = self.gene_id
        else:
            entrez_id = self.convert_gene_to_entrez_id()
            if entrez_id != 'not_found':
                pass
            else:
                result_promoter = f'Please verify if {self.gene_id} exist for {self.species}'
                return result_promoter

        gene_info = get_gene_info(entrez_id)
        if 'chraccver' in str(gene_info):
            gene_name = gene_info['name']
            chraccver = gene_info['genomicinfo'][0]['chraccver']
            chrstart = int(gene_info['genomicinfo'][0]['chrstart'])
            chrstop = int(gene_info['genomicinfo'][0]['chrstop'])
            species_API = gene_info['organism']['scientificname']
        else:
            result_promoter = f'Please verify ID of {self.gene_id}'
            return result_promoter

        dna_sequence = self.get_dna_sequence(chraccver, chrstart, chrstop)

        if self.prom_term == 'Promoter':
            result_promoter = f">{gene_name} | {species_API} | {chraccver} | {self.prom_term} | TSS (on chromosome): {chrstart} | TSS (on sequence): {self.upstream}\n{dna_sequence}\n"
        else:
            result_promoter = f">{gene_name} | {species_API} | {chraccver} | {self.prom_term} | Gene end (on chromosome): {chrstop} | Gene end (on sequence): {self.upstream}\n{dna_sequence}\n"

        return result_promoter

    # Convert gene to ENTREZ_GENE_ID
    def convert_gene_to_entrez_id(self):
        if self.gene_id.isdigit():
            return gene  # Already an ENTREZ_GENE_ID

        # Request for ENTREZ_GENE_ID
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={self.gene_id}[Gene%20Name]+AND+{self.species}[Organism]&retmode=json&rettype=xml "
        response = requests.get(url)

        if response.status_code == 200:
            response_data = response.json()

            if response_data['esearchresult']['count'] == '0':
                gene_id = 'not_found'
                return gene_id

            else:
                gene_id = response_data['esearchresult']['idlist'][0]
                return gene_id

    # Get DNA sequence
    def get_dna_sequence(self, chraccver, chrstart, chrstop):
        # Determine sens of gene + coordinate for upstream and downstream
        if chrstop > chrstart:
            start = (chrstart if self.prom_term == 'Promoter' else chrstop) - self.upstream
            end = (chrstart if self.prom_term == 'Promoter' else chrstop) + self.downstream
        else:
            start = (chrstart if self.prom_term == 'Promoter' else chrstop) + self.upstream
            end = (chrstart if self.prom_term == 'Promoter' else chrstop) - self.downstream

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