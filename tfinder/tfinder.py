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

class NCBI_dna:
    species_list = ['Human', 'Mouse', 'Rat', 'Drosophila', 'Zebrafish']
    def __init__(self):
        self.gene_name = gene_name

    def reverse_complement(sequence):
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        reverse_sequence = sequence[::-1]
        complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
        return complement_sequence

    @classmethod
    def analyse_gene(gene_id):
        disponibility_list = ['ID', 'Human', 'Mouse', 'Rat', 'Drosophila', 'Zebrafish']
        time.sleep(0.25)
        gene_analyse = [gene_id]
        for species_test in disponibility_list:
            if not gene_id.isdigit():
                if species_test == 'ID':
                    gene_analyse.append('n.d')
                else:
                    time.sleep(0.5)
                    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_id}[Gene%20Name]+AND+{species_test}[Organism]&retmode=json&rettype=xml"
                    response = requests.get(url)

                    if response.status_code == 200:
                        response_data = response.json()

                        if response_data['esearchresult']['count'] != '0':
                            gene_analyse.append("✅")
                        else:
                            gene_analyse.append("❌")

            if gene_id.isdigit():
                if species_test != 'ID':
                    gene_analyse.append('n.d')
                else:
                    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json&rettype=xml"
                    response = requests.get(url)

                    if response.status_code == 200:
                        response_data = response.json()

                        if 'chraccver' in str(response_data):
                            gene_analyse.append("✅")
                        else:
                            gene_analyse.append("❌")

        return gene_analyse

    # Convert gene to ENTREZ_GENE_ID
    def convert_gene_to_entrez_id(self, gene, species):
        if gene.isdigit():
            return gene  # Already an ENTREZ_GENE_ID

        # Request for ENTREZ_GENE_ID
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene}[Gene%20Name]+AND+{species}[Organism]&retmode=json&rettype=xml "
        response = requests.get(url)

        if response.status_code == 200:
            response_data = response.json()

            if response_data['esearchresult']['count'] == '0':
                st.error(f"No gene found for: {gene} from {species}")
                gene_id = 'not_found'
                return gene_id

            else:
                gene_id = response_data['esearchresult']['idlist'][0]
                return gene_id

        else:
            raise Exception(f"Error during gene search: {response.status_code}")

    # Get gene information
    def get_gene_info(self, gene_id):
        try:
            # Request gene information
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json&rettype=xml"
            response = requests.get(url)

            if response.status_code == 200:
                response_data = response.json()
                gene_info = response_data['result'][str(gene_id)]
                return gene_info

        except Exception as e:
            raise Exception(f"Error: {str(e)}")

    # Get DNA sequence
    def get_dna_sequence(self, chraccver, chrstart, chrstop, upstream, downstream):
        try:
            # Determine sens of gene + coordinate for upstream and downstream
            if chrstop > chrstart:
                start = (chrstart if prom_term == 'Promoter' else chrstop) - upstream
                end = (chrstart if prom_term == 'Promoter' else chrstop) + downstream
            else:
                start = (chrstart if prom_term == 'Promoter' else chrstop) + upstream
                end = (chrstart if prom_term == 'Promoter' else chrstop) - downstream

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
    def find_promoters(self, gene_ids, species, upstream, downstream):
        try:
            for gene_id in gene_ids:
                time.sleep(1)
                if gene_id.isdigit():
                    entrez_id = gene_id
                else:
                    entrez_id = convert_gene_to_entrez_id(gene_id, species)
                    if entrez_id != 'not_found':
                        pass
                    else:
                        continue

                gene_info = get_gene_info(entrez_id)
                if 'chraccver' in str(gene_info):
                    gene_name = gene_info['name']
                    chraccver = gene_info['genomicinfo'][0]['chraccver']
                    chrstart = gene_info['genomicinfo'][0]['chrstart']
                    chrstop = gene_info['genomicinfo'][0]['chrstop']
                    species_API = gene_info['organism']['scientificname']
                else:
                    st.error(f'Please verify ID of {gene_id}')
                    continue

                dna_sequence = get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream)

                st.toast(f'{prom_term} **{gene_name}** from **{species_API}** extracted', icon='🧬')

                # Append the result to the result_promoter
                if prom_term == 'Promoter':
                    result_promoter.append(
                        f">{gene_name} | {species_API} | {chraccver} | {prom_term} | TSS (on chromosome): {chrstart} | TSS (on sequence): {upstream}\n{dna_sequence}\n")
                else:
                    result_promoter.append(
                        f">{gene_name} | {species_API} | {chraccver} | {prom_term} | Gene end (on chromosome): {chrstop} | Gene end (on sequence): {upstream}\n{dna_sequence}\n")

            return result_promoter

        except Exception as e:
            raise Exception(f"Error retrieving gene information: {entrez_id}")