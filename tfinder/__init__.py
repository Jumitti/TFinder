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

import logomaker
import numpy as np
import random
import requests
import time
from bs4 import BeautifulSoup


class NCBIdna:
    def __init__(self,
                 gene_id,
                 prom_term,
                 upstream,
                 downstream,
                 species=None):
        self.gene_id = gene_id
        self.prom_term = prom_term if prom_term is not None else None
        self.upstream = upstream if upstream is not None else None
        self.downstream = downstream if downstream is not None else None
        self.species = species if species is not None else None


    @staticmethod
    def XMNM_to_gene_ID(variant):
        url = f"https://www.ncbi.nlm.nih.gov/nuccore/{variant}"

        # Envoyer une requête GET pour récupérer le contenu de la page
        response = requests.get(url)

        # Vérifier si la requête a réussi (code 200)
        if response.status_code == 200:
            # Analyser le contenu HTML de la page
            soup = BeautifulSoup(response.text, 'html.parser')
            print(soup)

            # Trouver l'élément HTML contenant le titre de la page
            title_element = soup.find("title")

            # Extraire le texte du titre
            page_title = title_element.get_text() if title_element else "Titre non trouvé"

            # Imprimer le titre
            print("Titre de la page:", page_title)
        else:
            print("Erreur lors de la requête HTTP.")

    @staticmethod
    # Analyse if gene is available
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

    # Sequence extractor
    def find_sequences(self):
        time.sleep(1)
        if self.gene_id.isdigit():
            entrez_id = self.gene_id
        elif self.gene_id.startwith('XM_') or self.gene_id.startwith('NM_'):
            print('hello')
        else:
            entrez_id = NCBIdna.convert_gene_to_entrez_id(self.gene_id, self.species)
            if entrez_id != 'not_found':
                pass
            else:
                result_promoter = f'Please verify if {self.gene_id} exist for {self.species}'
                return result_promoter

        gene_info = NCBIdna.get_gene_info(entrez_id)
        if 'chraccver' in str(gene_info):
            gene_name = gene_info['name']
            chraccver = gene_info['genomicinfo'][0]['chraccver']
            chrstart = int(gene_info['genomicinfo'][0]['chrstart'])
            chrstop = int(gene_info['genomicinfo'][0]['chrstop'])
            species_API = gene_info['organism']['scientificname']
        else:
            result_promoter = f'Please verify ID of {self.gene_id}'
            return result_promoter

        prom_term = self.prom_term.lower()
        if prom_term not in ['promoter', 'terminator']:
            result_promoter = f"'{self.prom_term}' not valid. Please use 'Promoter' or 'Terminator'."
            return result_promoter

        if isinstance(self.upstream, int) and isinstance(self.downstream, int):
            upstream = int(self.upstream)
            downstream = int(self.downstream)
        else:
            result_window = f'Upstream {self.upstream} and Downstream {self.downstream} must be integer'
            return result_window

        dna_sequence = NCBIdna.get_dna_sequence(prom_term, upstream, downstream, chraccver, chrstart, chrstop)

        if prom_term == 'promoter':
            dna_sequence = f">{gene_name} | {species_API} | {chraccver} | {self.prom_term} | TSS (on chromosome): {chrstart} | TSS (on sequence): {self.upstream}\n{dna_sequence}"
        else:
            dna_sequence = f">{gene_name} | {species_API} | {chraccver} | {self.prom_term} | Gene end (on chromosome): {chrstop} | Gene end (on sequence): {self.upstream}\n{dna_sequence}"

        return dna_sequence

    @staticmethod
    # Convert gene to ENTREZ_GENE_ID
    def convert_gene_to_entrez_id(gene_name, species):
        if gene_name.isdigit():
            gene_name = 'Already gene ID'
            return gene_name  # Already an ENTREZ_GENE_ID

        # Request for ENTREZ_GENE_ID
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_name}[Gene%20Name]+AND+{species}[Organism]&retmode=json&rettype=xml "
        response = requests.get(url)

        if response.status_code == 200:
            response_data = response.json()

            if response_data['esearchresult']['count'] == '0':
                gene_name = 'not_found'
                return gene_name

            else:
                gene_name = response_data['esearchresult']['idlist'][0]
                return gene_name

    @staticmethod
    # Get gene information
    def get_gene_info(entrez_id):
        # Request gene information
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={entrez_id}&retmode=json&rettype=xml"
        response = requests.get(url)

        if response.status_code == 200:
            response_data = response.json()
            gene_info = response_data['result'][str(entrez_id)]
            if 'chraccver' in str(gene_info):
                return gene_info
            else:
                gene_info = int(str('0'))

                return gene_info

    @staticmethod
    # Get DNA sequence
    def get_dna_sequence(prom_term, upstream, downstream, chraccver, chrstart, chrstop):
        # Determine sens of gene + coordinate for upstream and downstream
        if chrstop > chrstart:
            start = (chrstart if prom_term.lower() == 'promoter' else chrstop) - upstream
            end = (chrstart if prom_term.lower() == 'promoter' else chrstop) + downstream
        else:
            start = (chrstart if prom_term.lower() == 'promoter' else chrstop) + upstream
            end = (chrstart if prom_term.lower() == 'promoter' else chrstop) - downstream

        # Request for DNA sequence
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&from={start}&to={end}&rettype=fasta&retmode=text"
        response = requests.get(url)

        if response.status_code == 200:
            # Extraction of DNA sequence
            dna_sequence = response.text.split('\n', 1)[1].replace('\n', '')
            if chrstop > chrstart:
                sequence = dna_sequence
            else:
                sequence = NCBIdna.reverse_complement(dna_sequence)

            return sequence

    @staticmethod
    def reverse_complement(dna_sequence):
        DNA_code = ["A", "T", "C", "G", "N", "a", "t", "c", "g", "n"]
        if not all(char in DNA_code for char in dna_sequence):
            isdna = 'Please use only A T G C'
            return isdna
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        reverse_sequence = dna_sequence[::-1].upper()
        complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
        return complement_sequence


class IMO:

    def __init__(self):
        pass

    @staticmethod
    # Extract JASPAR matrix
    def matrix_extraction(jaspar_id):
        url = f"https://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
        response = requests.get(url)
        if response.status_code == 200:
            response_data = response.json()
            TF_name = response_data['name']
            TF_species = response_data['species'][0]['name']
            matrix = response_data['pfm']
            weblogo = f"https://jaspar.genereg.net/static/logos/all/svg/{jaspar_id}.svg"
        else:
            TF_name = 'not found'
            TF_species = 'not found'
            matrix = 'not found'
            weblogo = 'not found'
        return TF_name, TF_species, matrix, weblogo

    @staticmethod
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

    @staticmethod
    # Generate random sequences for p_value
    def generate_ranseq(probabilities, seq_length, progress_bar, num_random_seqs):
        motif_length = seq_length
        random_sequences = []

        for _ in range(num_random_seqs):
            random_sequence = IMO.generate_random_sequence(motif_length, probabilities)
            random_sequences.append(random_sequence)
            progress_bar.update(1)

        return random_sequences

    @staticmethod
    # Calculate matrix score
    def calculate_score(sequence, matrix):
        score = 0
        for i, base in enumerate(sequence):
            if base in {'A', 'C', 'G', 'T'}:
                base_score = matrix[base]
                score += base_score[i]
        return score

    @staticmethod
    # Generate random sequences
    def generate_random_sequence(length, probabilities):
        nucleotides = ['A', 'C', 'G', 'T']
        sequence = random.choices(nucleotides, probabilities, k=length)
        return ''.join(sequence)

    @staticmethod
    # Analyse sequence for non-authorized characters
    def is_dna(dna_sequence):
        DNA_code = ["A", "T", "C", "G", "N", "a", "t", "c", "g", "n"]
        if not all(char in DNA_code for char in dna_sequence):
            isfasta = True
            return isfasta
        else:
            isfasta = False
            return isfasta

    @staticmethod
    # Find with JASPAR and manual matrix
    def individual_motif_finder(dna_sequences, threshold, matrix, progress_bar, calc_pvalue=None, tss_ge_distance=None):
        if calc_pvalue is not None:
            if calc_pvalue not in ["ATGCPreset", "ATGCProportion"]:
                raise ValueError("Use 'ATGCPreset' or 'ATGCProportion'")

        individual_motif_occurrences = []

        matrices = IMO.transform_matrix(matrix)

        num_random_seqs = 1000000

        seq_length = len(matrices['Original']['A'])

        if calc_pvalue == 'ATGCPreset':
            percentage_a = 0.275
            percentage_t = 0.275
            percentage_g = 0.225
            percentage_c = 0.225

            probabilities = [percentage_a, percentage_c, percentage_g, percentage_t]

            random_sequences = IMO.generate_ranseq(probabilities, seq_length, progress_bar, num_random_seqs)

        if calc_pvalue == 'ATGCPreset':
            random_scores = {}
            matrix_random_scores = []
            for matrix_name, matrix in matrices.items():
                max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                min_score = sum(min(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                for random_sequence in random_sequences:
                    sequence = random_sequence
                    random_score = IMO.calculate_score(sequence, matrix)
                    if max_score == min_score:
                        normalized_random_score = random_score / max_score
                    else:
                        normalized_random_score = (random_score - min_score) / (max_score - min_score)
                    matrix_random_scores.append(normalized_random_score)
                    progress_bar.update(1)

            random_scores = np.array(matrix_random_scores)

        for name, dna_sequence, species, region in dna_sequences:
            if calc_pvalue == 'ATGCProportion':
                count_a = dna_sequence.count('A')
                count_t = dna_sequence.count('T')
                count_g = dna_sequence.count('G')
                count_c = dna_sequence.count('C')

                length_prom = len(dna_sequence)
                percentage_a = count_a / length_prom
                percentage_t = count_t / length_prom
                percentage_g = count_g / length_prom
                percentage_c = count_c / length_prom

                probabilities = [percentage_a, percentage_c, percentage_g, percentage_t]

                random_sequences = IMO.generate_ranseq(probabilities, seq_length, progress_bar, num_random_seqs)

                if calc_pvalue == 'ATGCProportion':
                    random_scores = {}

            for matrix_name, matrix in matrices.items():
                found_positions = []

                # Max score per matrix
                max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                min_score = sum(min(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))

                if calc_pvalue == 'ATGCProportion':
                    matrix_random_scores = []
                    for random_sequence in random_sequences:
                        sequence = random_sequence
                        random_score = IMO.calculate_score(sequence, matrix)
                        if max_score == min_score:
                            normalized_random_score = random_score / max_score
                        else:
                            normalized_random_score = (random_score - min_score) / (max_score - min_score)
                        matrix_random_scores.append(normalized_random_score)
                        progress_bar.update(1)

                    random_scores = np.array(matrix_random_scores)

                for i in range(len(dna_sequence) - seq_length + 1):
                    seq = dna_sequence[i:i + seq_length]
                    score = IMO.calculate_score(seq, matrix)
                    if max_score == min_score:
                        normalized_score = score / max_score
                    else:
                        normalized_score = (score - min_score) / (max_score - min_score)
                    position = int(i) + 1

                    found_positions.append((position, seq, normalized_score))
                    progress_bar.update(1)

                # Sort positions in descending order of score percentage
                found_positions.sort(key=lambda x: x[1], reverse=True)

                # Creating a results table
                if len(found_positions) > 0:
                    if threshold < 0.5:
                        highest_normalized_score = max(
                            [normalized_score for _, _, normalized_score in found_positions])
                        if highest_normalized_score >= 0.6:
                            threshold = highest_normalized_score - 0.10
                        else:
                            threshold = 0.5

                    for position, seq, normalized_score in found_positions:
                        start_position = max(0, position - 4)
                        end_position = min(len(dna_sequence), position + len(seq) + 2)
                        sequence_with_context = dna_sequence[start_position:end_position]

                        sequence_parts = []
                        for j in range(start_position, end_position):
                            if j < position - 1 or j + 2 > position + len(seq):
                                sequence_parts.append(sequence_with_context[j - start_position].lower())
                            else:
                                sequence_parts.append(sequence_with_context[j - start_position].upper())

                        sequence_with_context = ''.join(sequence_parts)

                        if tss_ge_distance is not None:
                            tis_position = position - tss_ge_distance

                        if normalized_score >= threshold:
                            if calc_pvalue is not None:
                                p_value = (random_scores >= normalized_score).sum() / len(random_scores)

                            row = [str(position).ljust(8)]
                            if tss_ge_distance is not None:
                                row.append(str(tis_position).ljust(15))
                            row += [sequence_with_context,
                                    "{:.6f}".format(normalized_score).ljust(12)]
                            if calc_pvalue is not None:
                                row.append("{:.3e}".format(p_value).ljust(12))
                            row += [name, species, region]
                            individual_motif_occurrences.append(row)

        if len(individual_motif_occurrences) > 0:
            if tss_ge_distance is not None:
                individual_motif_occurrences.sort(key=lambda x: float(x[3]), reverse=True)
            else:
                individual_motif_occurrences.sort(key=lambda x: float(x[2]), reverse=True)
            header = ["Position"]
            if tss_ge_distance is not None:
                header.append("Rel Position")
            header += ["Sequence", "Rel Score"]
            if calc_pvalue is not None:
                header.append("p-value")
            header += ["Gene", "Species", "Region"]
            individual_motif_occurrences.insert(0, header)
        else:
            "No consensus sequence found with the specified threshold."

        return individual_motif_occurrences

    @staticmethod
    # IUPAC code
    def generate_iupac_variants(sequence, max_variant_allowed=None, progress_bar=None):
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
            "N": ["A", "C", "G", "T"],
            "-": ['-'],
            ".": ['-']
        }

        iupac_codes_score = {
            "R": 2,  # A or G
            "Y": 2,  # C or T
            "M": 2,  # A or C
            "K": 2,  # G or T
            "W": 2,  # A or T
            "S": 2,  # C or G
            "B": 3,  # C or G or T
            "D": 3,  # A or G or T
            "H": 3,  # A or C or T
            "V": 3,  # A or C or G
            "N": 4  # A or C or G or T
        }

        total_variants = 1
        for base in sequence:
            if base.upper() in iupac_codes_score:
                total_variants *= iupac_codes_score[base.upper()]
        if total_variants > max_variant_allowed:
            sequence = f'Too many variants. Use - or . instead N. Limit: {max_variant_allowed} | Total variants : {total_variants}'
            return sequence

        sequences = [sequence]
        for i, base in enumerate(sequence):
            if base.upper() in iupac_codes:
                new_sequences = []
                for seq in sequences:
                    for alternative in iupac_codes[base.upper()]:
                        new_sequence = seq[:i] + alternative + seq[i + 1:]
                        new_sequences.append(new_sequence)
                        if progress_bar is not None:
                            progress_bar.update(1)
                sequences = new_sequences

        return sequences

    @staticmethod
    # is PWM good ?
    def has_uniform_column_length(matrix_str):
        lines = matrix_str.strip().split('\n')
        values_lists = [line.split('[')[1].split(']')[0].split() for line in lines]
        column_lengths = set(len(values) for values in values_lists)
        if len(column_lengths) != 1:
            raise Exception('Invalid PWM length.')

    @staticmethod
    # Calculate PWM
    def calculate_pwm(sequences):
        num_sequences = len(sequences)
        sequence_length = len(sequences[0])
        pwm = np.zeros((4, sequence_length))
        for i in range(sequence_length):
            counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
            gap_found = False
            for sequence in sequences:
                nucleotide = sequence[i]
                if nucleotide == '_':
                    gap_found = True
                    break
                elif nucleotide in counts:
                    counts[nucleotide] += 1

            if gap_found:
                pwm[0, i] = 0
                pwm[1, i] = 0
                pwm[2, i] = 0
                pwm[3, i] = 0
            else:
                pwm[0, i] = counts['A'] / num_sequences
                pwm[1, i] = counts['T'] / num_sequences
                pwm[2, i] = counts['G'] / num_sequences
                pwm[3, i] = counts['C'] / num_sequences

        bases = ['A', 'T', 'G', 'C']
        pwm_text = ""
        for i in range(len(pwm)):
            base_name = bases[i]
            base_values = pwm[i]

            base_str = base_name + " ["
            for value in base_values:
                base_str += " " + format(value) + " " if np.isfinite(value) else " " + "NA" + " "

            base_str += "]\n"
            pwm_text += base_str

        matrix_lines = pwm_text.split('\n')
        matrix = {}
        for line in matrix_lines:
            line = line.strip()
            if line:
                key, values = line.split('[', 1)
                values = values.replace(']', '').split()
                values = [float(value) for value in values]
                matrix[key.strip()] = values
        return matrix

    @staticmethod
    # PWM with multiple FASTA
    def parse_fasta(individual_motif):
        sequences = []
        current_sequence = ""

        for line in individual_motif.splitlines():
            if line.startswith(">"):
                if current_sequence:
                    sequences.append(current_sequence)
                current_sequence = ""
            else:
                current_sequence += line

        if current_sequence:
            sequences.append(current_sequence)

        return sequences

    @staticmethod
    # generate Weblogo
    def create_web_logo(sequences):
        matrix = logomaker.alignment_to_matrix(sequences)
        logo = logomaker.Logo(matrix, color_scheme='classic')
        return logo

    @staticmethod
    # Individual motif PWM and weblogo
    def individual_motif_pwm(individual_motif):
        sequences = IMO.parse_fasta(individual_motif)
        sequences = [seq.upper() for seq in sequences]

        if len(sequences) > 0:
            sequence_length = len(sequences[0])
            inconsistent_lengths = False

            for sequence in sequences[1:]:
                if len(sequence) != sequence_length:
                    inconsistent_lengths = True
                    break

            if inconsistent_lengths:
                raise Exception(f"Sequence lengths are not consistent.")

            else:
                matrix = IMO.calculate_pwm(sequences)

                sequences_text = individual_motif
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

                sequences.append(current_sequence)

                weblogo = IMO.create_web_logo(sequences)

                return matrix, weblogo

        else:
            raise Exception(f"You forget FASTA sequences :)")
