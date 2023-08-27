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


    # Analyse if gene is available
    def analyse_gene(self):
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

        gene_info = self.get_gene_info(entrez_id)
        if 'chraccver' in str(gene_info):
            gene_name = gene_info['name']
            chraccver = gene_info['genomicinfo'][0]['chraccver']
            chrstart = int(gene_info['genomicinfo'][0]['chrstart'])
            chrstop = int(gene_info['genomicinfo'][0]['chrstop'])
            species_API = gene_info['organism']['scientificname']
        else:
            result_promoter = f'Please verify ID of {self.gene_id}'
            return result_promoter

        prom_term = self.prom_term
        upstream = self.upstream
        downstream = self.downstream

        dna_sequence = self.get_dna_sequence(prom_term, upstream, downstream, chraccver, chrstart, chrstop)

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
                sequence = NCBI_dna.reverse_complement(dna_sequence)

            return sequence
    @staticmethod
    def reverse_complement(dna_sequence):
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        reverse_sequence = dna_sequence[::-1]
        complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
        return complement_sequence

class IMO:
    def __init__(self):
        ok = ok

    # Extract JASPAR matrix
    def matrix_extraction(sequence_consensus_input):
        jaspar_id = sequence_consensus_input
        url = f"https://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
        response = requests.get(url)
        if response.status_code == 200:
            response_data = response.json()
            matrix = response_data['pfm']
        else:
            st.error(f"Error while retrieving PWM: {response.status_code}")
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

    # Generate random sequences for p_value
    def generate_ranseq(probabilities, seq_length, pbar, num_random_seqs):
        motif_length = seq_length
        random_sequences = []

        for _ in range(num_random_seqs):
            random_sequence = generate_random_sequence(motif_length, probabilities)
            random_sequences.append(random_sequence)
            pbar.update(1)

        return random_sequences

    # Calculate matrix score
    def calculate_score(sequence, matrix):
        score = 0
        for i, base in enumerate(sequence):
            if base in {'A', 'C', 'G', 'T'}:
                base_score = matrix[base]
                score += base_score[i]
        return score

    # Generate random sequences
    def generate_random_sequence(length, probabilities):
        nucleotides = ['A', 'C', 'G', 'T']
        sequence = random.choices(nucleotides, probabilities, k=length)
        return ''.join(sequence)

    # Analyse sequence for non authorized characters
    def isdna(promoter_region):
        DNA_code = ["A", "T", "C", "G", "N", "a", "t", "c", "g", "n"]
        if not all(char in DNA_code for char in promoter_region):
            isfasta = True
            return isfasta
        else:
            isfasta = False
            return isfasta

    # Find with JASPAR and manual matrix
    def search_sequence(threshold, tis_value, promoters, matrices, total_promoter_region_length, total_promoter):
        global table2
        table2 = []

        seq_length = len(matrices['Original']['A'])
        sequence_iteration = len(matrices.items()) * total_promoter_region_length

        num_random_seqs = 1000000
        if total_promoter <= 10:
            random_gen = len(promoters) * num_random_seqs
        else:
            random_gen = num_random_seqs
        random_score = random_gen * len(matrices.items())

        if calc_pvalue:
            total_iterations = sequence_iteration + random_gen + random_score
        else:
            total_iterations = sequence_iteration

        with stqdm(total=total_iterations,
                   desc='**:blue[Processing...] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**',
                   mininterval=0.1) as pbar:
            if calc_pvalue and total_promoter > 10:
                percentage_a = 0.275
                percentage_t = 0.275
                percentage_g = 0.225
                percentage_c = 0.225

                probabilities = [percentage_a, percentage_c, percentage_g, percentage_t]

                random_sequences = generate_ranseq(probabilities, seq_length, pbar, num_random_seqs)

            if calc_pvalue and total_promoter > 10:
                random_scores = {}
                matrix_random_scores = []
                for matrix_name, matrix in matrices.items():
                    max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                    min_score = sum(min(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                    for random_sequence in random_sequences:
                        sequence = random_sequence
                        random_score = calculate_score(sequence, matrix)
                        normalized_random_score = (random_score - min_score) / (max_score - min_score)
                        matrix_random_scores.append(normalized_random_score)
                        pbar.update(1)

                random_scores = np.array(matrix_random_scores)

            for shortened_promoter_name, promoter_region, found_species, region in promoters:
                if calc_pvalue and total_promoter <= 10:
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

                    random_sequences = generate_ranseq(probabilities, seq_length, pbar, num_random_seqs)

                    if calc_pvalue and total_promoter <= 10:
                        random_scores = {}

                for matrix_name, matrix in matrices.items():
                    found_positions = []

                    # Max score per matrix
                    max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                    min_score = sum(min(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))

                    if calc_pvalue and total_promoter <= 10:
                        matrix_random_scores = []
                        for random_sequence in random_sequences:
                            sequence = random_sequence
                            random_score = calculate_score(sequence, matrix)
                            normalized_random_score = (random_score - min_score) / (max_score - min_score)
                            matrix_random_scores.append(normalized_random_score)
                            pbar.update(1)

                        random_scores = np.array(matrix_random_scores)

                    for i in range(len(promoter_region) - seq_length + 1):
                        seq = promoter_region[i:i + seq_length]
                        score = calculate_score(seq, matrix)
                        normalized_score = (score - min_score) / (max_score - min_score)
                        position = int(i)

                        found_positions.append((position, seq, normalized_score))
                        pbar.update(1)

                    # Sort positions in descending order of score percentage
                    found_positions.sort(key=lambda x: x[1], reverse=True)

                    # Creating a results table
                    if len(found_positions) > 0:
                        if auto_thre:
                            highest_normalized_score = max(
                                [normalized_score for _, _, normalized_score in found_positions])
                            if highest_normalized_score >= 0.6:
                                threshold = highest_normalized_score - 0.10
                            else:
                                threshold = 0.5

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
                                if calc_pvalue:
                                    p_value = (random_scores >= normalized_score).sum() / len(random_scores)

                                row = [str(position).ljust(8),
                                       str(tis_position).ljust(15),
                                       sequence_with_context,
                                       "{:.6f}".format(normalized_score).ljust(12)]
                                if calc_pvalue:
                                    row.append("{:.3e}".format(p_value).ljust(12))
                                row += [shortened_promoter_name, found_species, region]
                                table2.append(row)

        if len(table2) > 0:
            table2.sort(key=lambda x: float(x[3]), reverse=True)
            header = ["Position", "Rel Position", "Sequence", "Rel Score"]
            if calc_pvalue:
                header.append("p-value")
            header += ["Gene", "Species", "Region"]
            table2.insert(0, header)
        else:
            "No consensus sequence found with the specified threshold."

        return table2

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

    # is PWM good ?
    def has_uniform_column_length(pwm):
        column_lengths = set(len(column) for column in pwm)
        if len(column_lengths) != 1:
            raise Exception('Invalid PWM lenght.')

    # Calculate PWM
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

    # PWM with multiple FASTA
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

    # generate Weblogo
    def create_web_logo(sequences):
        matrix = logomaker.alignment_to_matrix(sequences)
        logo = logomaker.Logo(matrix, color_scheme='classic')
        return logo

    # Individual motif PWM and weblogo
    def im(fasta_text):
        sequences = parse_fasta(fasta_text)
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
                pwm = calculate_pwm(sequences)
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

                with REcol2:
                    st.markdown("PWM", help="Modification not allowed. Still select and copy for later use.")
                    matrix_text = st.text_area("PWM:", value=pwm_text,
                                               label_visibility='collapsed',
                                               height=125,
                                               disabled=True)

                with REcol2:
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

                    sequences.append(current_sequence)

                    logo = create_web_logo(sequences)
                    st.pyplot(logo.fig)
                    buffer = io.BytesIO()
                    logo.fig.savefig(buffer, format='png')
                    buffer.seek(0)

                    st.session_state['buffer'] = buffer

                    return matrix_text, buffer

        else:
            raise Exception(f"You forget FASTA sequences :)")