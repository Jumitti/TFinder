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
import io
from openpyxl import Workbook
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.image import MIMEImage
from email import encoders
import base64
import datetime
import matplotlib.pyplot as plt
from PIL import Image
import time
from tqdm import tqdm
from stqdm import stqdm


def aio_page():
    # Reverse complement
    def reverse_complement(sequence):
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        reverse_sequence = sequence[::-1]
        complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
        return complement_sequence

    # Convert gene to ENTREZ_GENE_ID
    def convert_gene_to_entrez_id(gene, species):
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
    def get_gene_info(gene_id):
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
    def get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream):
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
    def find_promoters(gene_ids, species, upstream, downstream):
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
    def search_sequence(threshold, tis_value, promoters, matrices, total_promoter_region_length):
        global table2
        table2 = []

        for matrix_name, matrix in matrices.items():
            seq_length = len(matrix['A'])

        sequence_iteration = len(matrices.items()) * total_promoter_region_length
        random_gen = len(promoters) * 1000000
        random_score = random_gen * len(matrices.items())

        if calc_pvalue:
            total_iterations = sequence_iteration + random_gen + random_score
        else:
            total_iterations = sequence_iteration

        with stqdm(total=total_iterations, desc='Calculating scores', mininterval=0.1) as pbar:

            if calc_pvalue:
                for shortened_promoter_name, promoter_region, found_species, region in promoters:

                    # Generate random sequences
                    motif_length = seq_length
                    num_random_seqs = 1000000

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

                    random_sequences = []
                    for _ in range(num_random_seqs):
                        random_sequence = generate_random_sequence(motif_length, probabilities)
                        random_sequences.append(random_sequence)
                        pbar.update(1)

                    # Calculation of random scores from the different matrices
                    random_scores = {}

            for matrix_name, matrix in matrices.items():
                seq_length = len(matrix['A'])

                # Max score per matrix
                max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                min_score = sum(min(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))

                # REF
                for shortened_promoter_name, promoter_region, found_species, region in promoters:
                    found_positions = []

                    if calc_pvalue:

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

                        if calc_pvalue:
                            p_value = (random_scores >= normalized_score).sum() / num_random_seqs
                        else:
                            p_value = 0

                        found_positions.append((position, seq, normalized_score, p_value))
                        pbar.update(1)

                    # Sort positions in descending order of score percentage
                    found_positions.sort(key=lambda x: x[1], reverse=True)

                    # Creating a results table
                    if len(found_positions) > 0:
                        if auto_thre:
                            highest_normalized_score = max(
                                [normalized_score for _, _, normalized_score, _ in found_positions])
                            if highest_normalized_score >= 0.6:
                                threshold = highest_normalized_score - 0.10
                            else:
                                threshold = 0.5

                        for position, seq, normalized_score, p_value in found_positions:
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
                    st.markdown("PWM",  help="Modification not allowed. Still select and copy for later use.")
                    matrix_text = st.text_area("PWM:", value=pwm_text,
                                               label_visibility = 'collapsed',
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

    def email(excel_file, txt_output, email_receiver, body):
        try:
            current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            subject = f'Results TFinder - {current_date_time}'
            email_sender = st.secrets['sender']
            password = st.secrets['password']

            msg = MIMEMultipart()
            msg['From'] = email_sender
            msg['To'] = email_receiver
            msg['Subject'] = subject

            msg.attach(MIMEText(body, 'plain'))

            attachment_excel = MIMEBase('application', 'octet-stream')
            attachment_excel.set_payload(excel_file.getvalue())
            encoders.encode_base64(attachment_excel)
            attachment_excel.add_header('Content-Disposition', 'attachment',
                                        filename=f'Results_TFinder_{current_date_time}.xlsx')
            msg.attach(attachment_excel)

            if jaspar == 'PWM':
                if matrix_type == 'With FASTA sequences':
                    image = MIMEImage(st.session_state['buffer'].read(), name=f'LOGOMAKER_{current_date_time}.jpg')
                    msg.attach(image)
            elif jaspar == 'Manual sequence':
                image = MIMEImage(st.session_state['buffer'].read(), name=f'LOGOMAKER_{current_date_time}.jpg')
                msg.attach(image)

            attachment_text = MIMEText(txt_output, 'plain', 'utf-8')
            attachment_text.add_header('Content-Disposition', 'attachment',
                                       filename=f'Sequences_{current_date_time}.fasta')
            msg.attach(attachment_text)

            server = smtplib.SMTP('smtp.gmail.com', 587)
            server.starttls()
            server.login(email_sender, password)
            server.sendmail(email_sender, email_receiver, msg.as_string())
            server.quit()
            st.toast('Email sent successfully !', icon='🚀')

        except smtplib.SMTPAuthenticationError:
            st.toast("Failed to authenticate. Please check your email and password.")
        except smtplib.SMTPServerDisconnected:
            st.toast("Failed to connect to the SMTP server. Please check your internet connection.")
        except smtplib.SMTPRecipientsRefused:
            st.toast(f"Error sending email to {email_receiver}")
        except smtplib.SMTPException as e:
            st.toast(f"Error sending email: {e}")
        except Exception as e:
            st.toast(f"Unknown error occurred: {e}")

    def result_table_output(df):
        source = df
        score_range = source['Rel Score'].astype(float)
        ystart = score_range.min() - 0.02
        ystop = score_range.max() + 0.02
        source['Gene_Region'] = source['Gene'] + " " + source['Species'] + " " + source['Region']
        source['Beginning of sequences'] = source['Position']
        source['From TSS/gene end'] = source['Rel Position']
        scale = alt.Scale(scheme='category10')
        color_scale = alt.Color("Gene_Region:N", scale=scale)
        gene_region_selection = alt.selection_point(fields=['Gene_Region'], on='click', bind='legend')

        dropdown = alt.binding_select(options=['Beginning of sequences', 'From TSS/gene end'], name='(X-axis) Position from:')
        xcol_param = alt.param(value='Beginning of sequences', bind=dropdown)

        if 'p-value' in source:
            ispvalue = True

        chart = alt.Chart(source).mark_circle().encode(
            x=alt.X('x:Q').title('Position (bp)'),
            y=alt.Y('Rel Score:Q', axis=alt.Axis(title='Relative Score'),
                    scale=alt.Scale(domain=[ystart, ystop])),
            color=alt.condition(gene_region_selection, color_scale, alt.value('lightgray')),
            tooltip=['Position', 'Rel Position', 'Rel Score'] + (
                ['p-value'] if 'p-value' in source else []) + ['Sequence', 'Gene', 'Species', 'Region'],
            opacity=alt.condition(gene_region_selection, alt.value(0.8), alt.value(0.2))
        ).transform_calculate(x=f'datum[{xcol_param.name}]').properties(width=600, height=400).interactive().add_params(gene_region_selection, xcol_param)
        st.altair_chart(chart, theme=None, use_container_width=True)

    # Disposition
    st.subheader(':blue[Step 1] Promoter and Terminator Extractor')
    colprom1, colprom2 = st.columns([0.8, 1.2], gap="small")

    # Promoter Finder
    with colprom1:
        st.info("💡 If you have a FASTA sequence, go to :blue[**Step 2**]")

        result_promoter = []
        upstream_entry = []

        # Gene ID
        st.markdown("🔹 :blue[**Step 1.1**] Gene ID:", help='NCBI gene name and NCBI gene ID allowed')
        gene_id_entry = st.text_area("🔹 :blue[**Step 1.1**] Gene ID:", value="PRKN\n351",
                                     label_visibility='collapsed')
        gene_list = gene_id_entry.strip().split('\n')
        gene_ids = gene_id_entry.strip().split("\n")

        # Verify if gene is available for all species
        if st.button('🔎 Check genes avaibility',
                     help='Sometimes genes do not have the same name in all species or do not exist.'):
            species_list = ['Human', 'Mouse', 'Rat', 'Drosophila', 'Zebrafish']
            results_gene_list = []
            data = []
            for gene_input in stqdm(gene_list, desc="Analyse genes...", mininterval=0.1):
                time.sleep(0.25)
                if not gene_input.isdigit():
                    row = [gene_input]
                    for species_test in species_list:
                        time.sleep(0.5)
                        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_input}[Gene%20Name]+AND+{species_test}[Organism]&retmode=json&rettype=xml"
                        response = requests.get(url)

                        if response.status_code == 200:
                            response_data = response.json()

                            if response_data['esearchresult']['count'] != '0':
                                row.append("✅")
                            else:
                                row.append("❌")

                    data.append(row)

                if gene_input.isdigit():
                    gene_id = gene_input
                    gene_info = get_gene_info(gene_id)
                    if not 'chraccver' in str(gene_info):
                        st.error(f'Please verify ID of {gene_id}')

            species_columns = ['Gene'] + species_list
            dfgene = pd.DataFrame(data, columns=species_columns)
            st.session_state['dfgene'] = dfgene
        if 'dfgene' in st.session_state:
            st.dataframe(st.session_state['dfgene'], hide_index=True)

    with colprom2:
        tab1, tab2 = st.tabs(['Default', 'Advance'])

        with tab1:
            # Species
            st.markdown("🔹 :blue[**Step 1.2**] Select species of gene names:")
            species = st.selectbox("🔹 :blue[**Step 1.2**] Select species of gene names:",
                                   ["Human", "Mouse", "Rat", "Drosophila", "Zebrafish"], index=0,
                                   label_visibility='collapsed')

            # Upstream/Downstream Promoter
            st.markdown("🔹 :blue[**Step 1.3**] Regulatory region:")
            prom_term = st.radio("🔹 :blue[**Step 1.3**] Regulatory region:", ('Promoter', 'Terminator'),
                                 label_visibility='collapsed')
            if prom_term == 'Promoter':
                st.markdown("🔹 :blue[**Step 1.4**] Upstream/downstream from the TSS (bp)")
            else:
                st.markdown("🔹 :blue[**Step 1.4**] Upstream/downstream from gene end (bp)")

            updown_slide = st.slider("🔹 :blue[**Step 1.4**] Upstream/downstream", -10000, 10000,
                                     (-2000, 2000), step=100, label_visibility='collapsed')
            if prom_term == 'Promoter':
                st.write("Upstream: ", min(updown_slide), " bp from TSS | Downstream: ", max(updown_slide),
                         " bp from TSS")
            else:
                st.write("Upstream: ", min(updown_slide), " bp from gene end | Downstream: ", max(updown_slide),
                         " bp from gene end")

            upstream_entry = -min(updown_slide)
            downstream_entry = max(updown_slide)

            upstream = int(upstream_entry)
            st.session_state['upstream'] = upstream
            downstream = int(downstream_entry)

            # Run Promoter Finder
            if st.button(f"🧬 :blue[**Step 1.5**] Extract {prom_term}", help='(~5sec/gene)'):
                with colprom1:
                    with st.spinner("Finding promoters..."):
                        if 'result_promoter' in st.session_state:
                            del st.session_state['result_promoter']
                        try:
                            result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                            st.session_state['result_promoter'] = result_promoter
                            st.success(f"{prom_term} extraction complete !")
                            st.toast(f"{prom_term} extraction complete !", icon='😊')
                        except Exception as e:
                            st.error(f"Error finding {prom_term}: {str(e)}")

        with tab2:
            # Advance mode extraction
            data_df = pd.DataFrame(
                {
                    "Gene": gene_list,
                    "human": [False] * len(gene_list),
                    "mouse": [False] * len(gene_list),
                    "rat": [False] * len(gene_list),
                    "drosophila": [False] * len(gene_list),
                    "zebrafish": [False] * len(gene_list),
                    "promoter": [False] * len(gene_list),
                    "terminator": [False] * len(gene_list),
                }
            )

            species_list = ['human', 'mouse', 'rat', 'drosophila', 'zebrafish']
            search_types = ['promoter', 'terminator']

            st.markdown('**🔹 :blue[Step 1.2]** Select species for all genes:',
                        help='Checking a box allows you to check all the corresponding boxes for each gene. Warning: if you have manually checked boxes in the table, they will be reset.')

            species1, species2, species3, species4, species5 = st.columns(5)

            with species1:
                all_human = st.checkbox("Human")
            with species2:
                all_mouse = st.checkbox("Mouse")
            with species3:
                all_rat = st.checkbox("Rat")
            with species4:
                all_droso = st.checkbox("Drosophila")
            with species5:
                all_zebra = st.checkbox("Zebrafish")

            st.markdown('**🔹 :blue[Step 1.2]** Select regions for all genes:',
                        help='Checking a box allows you to check all the corresponding boxes for each gene. Warning: if you have manually checked boxes in the table, they will be reset.')

            region1, region2 = st.columns(2)

            with region1:
                all_prom = st.checkbox("Promoter")
            with region2:
                all_term = st.checkbox("Terminator")

            if all_human:
                data_df["human"] = True
            if all_mouse:
                data_df["mouse"] = True
            if all_rat:
                data_df["rat"] = True
            if all_droso:
                data_df["drosophila"] = True
            if all_zebra:
                data_df["zebrafish"] = True
            if all_prom:
                data_df["promoter"] = True
            if all_term:
                data_df["terminator"] = True

            st.markdown('**🔹 :blue[Step 1.2]** On demand genes table',
                        help="Check the boxes for which you want to extract a sequence. Pay attention that the gene name is equivalent for each species. The choice of species is not available for gene IDs. Parameterize the table last, if you check the boxes above, it resets the whole table.")

            data_dff = st.data_editor(
                data_df,
                column_config={
                    "human": st.column_config.CheckboxColumn(
                        "Human",
                        default=False,
                    ),
                    "mouse": st.column_config.CheckboxColumn(
                        "Mouse",
                        default=False,
                    ),
                    "rat": st.column_config.CheckboxColumn(
                        "Rat",
                        default=False,
                    ),
                    "drosophila": st.column_config.CheckboxColumn(
                        "Drosophila",
                        default=False,
                    ),
                    "zebrafish": st.column_config.CheckboxColumn(
                        "Zebrafish",
                        default=False,
                    ),
                    "promoter": st.column_config.CheckboxColumn(
                        "Promoter",
                        default=False,
                    ),
                    "terminator": st.column_config.CheckboxColumn(
                        "Terminator",
                        default=False,
                    )
                },
                disabled=["Gene"],
                hide_index=True,
            )

            updown_slide = st.slider("🔹 :blue[**Step 1.3**] Upstream/downstream from TSS and gene end (bp)",
                                     -10000,
                                     10000, (-2000, 2000), step=100, label_visibility='collapsed')
            st.write("Upstream: ", min(updown_slide), " bp from TSS and gene end | Downstream: ",
                     max(updown_slide),
                     " bp from TSS and gene end")
            upstream_entry = -min(updown_slide)
            downstream_entry = max(updown_slide)

            if st.button("🧬 :blue[**Step 1.4**] Extract sequences", help="(~5sec/seq)"):
                with colprom1:
                    with st.spinner("Finding sequences..."):
                        st.session_state['upstream'] = upstream_entry
                        upstream = int(upstream_entry)
                        downstream = int(downstream_entry)
                        for gene_info in data_dff.itertuples(index=False):
                            gene_name = gene_info.Gene
                            gene_ids = gene_name.strip().split('\n')
                            if gene_name.isdigit():
                                for search_type in search_types:
                                    if getattr(gene_info, f'{search_type}'):
                                        prom_term = search_type.capitalize()
                                        species = 'human'  # This is just a remnant of the past
                                        try:
                                            result_promoter = find_promoters(gene_ids, species, upstream,
                                                                             downstream)
                                            st.session_state['result_promoter'] = result_promoter
                                        except Exception as e:
                                            st.error(f"Error finding {gene_ids}: {str(e)}")
                            else:
                                for species in species_list:
                                    for search_type in search_types:
                                        if getattr(gene_info, f'{species}') and getattr(gene_info,
                                                                                        f'{search_type}'):
                                            prom_term = search_type.capitalize()
                                            try:
                                                result_promoter = find_promoters(gene_ids, species, upstream,
                                                                                 downstream)
                                                st.session_state['result_promoter'] = result_promoter
                                            except Exception as e:
                                                st.error(f"Error finding {gene_ids}: {str(e)}")

                        st.success(f"{prom_term} extraction complete !")
                        st.toast(f"{prom_term} extraction complete !", icon='😊')

    # Promoter output state
    st.divider()
    st.subheader(':blue[Step 2] Binding Sites Finder')
    promcol1, promcol2 = st.columns([0.9, 0.1], gap='small')
    with promcol1:
        st.markdown("🔹 :blue[**Step 2.1**] Sequences:", help='Copy: Click in sequence, CTRL+A, CTRL+C')
        if 'result_promoter' in st.session_state:
            result_promoter_text = "\n".join(st.session_state['result_promoter'])
        result_promoter = st.text_area("🔹 :blue[**Step 2.1**] Sequences:", value=result_promoter_text if 'result_promoter' in st.session_state else '', placeholder='If Step 1 not used, paste sequences here (FASTA required for multiple sequences).',
                                       key="result_promoter_key", label_visibility='collapsed')
    with promcol2:
        st.markdown('')
        st.markdown('')
        st.markdown('')
        current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        txt_output = f"{result_promoter}"
        st.download_button(label="💾 Download (.fasta)", data=txt_output,
                           file_name=f"Sequences_{current_date_time}.fasta", mime="text/plain")

    # Promoter detection information
    lines = result_promoter
    promoters = []
    if lines.startswith(("A", "T", "C", "G", "N", "a", "t", "c", "g", "n")):
        promoter_region = lines.upper()
        isfasta = isdna(promoter_region)
        shortened_promoter_name = "n.d."
        found_species = "n.d"
        region = "n.d"
        promoters.append((shortened_promoter_name, promoter_region, found_species, region))
    elif lines.startswith(">"):
        lines = result_promoter.split("\n")
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.startswith(">"):
                species_prom = ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Drosophila melanogaster',
                                'Danio rerio']
                promoter_name = line[1:]
                words = promoter_name.lstrip('>').split()
                shortened_promoter_name = words[0]
                for species in species_prom:
                    if species.lower() in promoter_name.lower():
                        found_species = species
                        break
                    else:
                        found_species = "n.d"
                regions_prom = ['Promoter', 'Terminator']
                for regions in regions_prom:
                    if regions.lower() in promoter_name.lower():
                        region = regions[:4] + "."
                        break
                    else:
                        region = "n.d"
                promoter_region = lines[i + 1].upper()
                isfasta = isdna(promoter_region)
                promoters.append((shortened_promoter_name, promoter_region, found_species, region))
                i += 1
            else:
                i += 1
    elif not lines.startswith(("A", "T", "C", "G", "N", "a", "t", "c", "g", "n", "I", "i", "")):
        isfasta = True
    else:
        isfasta = False

    total_promoter_region_length = sum(len(promoter_region) for _, promoter_region, _, _ in promoters)
    total_promoter = len(promoters)

    # RE entry
    REcol1, REcol2 = st.columns([0.30, 0.70])
    with REcol1:
        st.markdown('🔹 :blue[**Step 2.2**] Responsive elements type:')
        jaspar = st.radio('🔹 :blue[**Step 2.2**] Responsive elements type:', ('Manual sequence', 'JASPAR_ID', 'PWM'),
                          label_visibility='collapsed')
    if jaspar == 'JASPAR_ID':
        with REcol1:
            st.markdown("🔹 :blue[**Step 2.3**] JASPAR ID:")
            entry_sequence = st.text_input("🔹 :blue[**Step 2.3**] JASPAR ID:", value="MA0106.1",
                                           label_visibility='collapsed')
            url = f"https://jaspar.genereg.net/api/v1/matrix/{entry_sequence}/"
            response = requests.get(url)
            if response.status_code == 200:
                response_data = response.json()
                TF_name = response_data['name']
                TF_species = response_data['species'][0]['name']
                st.success(f"{TF_species} transcription factor {TF_name}")
                matrix = response_data['pfm']
                with REcol2:
                    st.image(f"https://jaspar.genereg.net/static/logos/all/svg/{entry_sequence}.svg")
                button = False
                error_input_im = True
            else:
                button = True
                error_input_im = False
                st.error('Wrong JASPAR_ID')

    elif jaspar == 'PWM':
        with REcol1:
            st.markdown('🔹 :blue[**Step 2.2bis**] Matrix:')
            matrix_type = st.radio('🔹 :blue[**Step 2.2bis**] Matrix:', ('With FASTA sequences', 'With PWM'),
                                   label_visibility='collapsed')
        if matrix_type == 'With PWM':
            isUIPAC = True
            with REcol2:
                st.markdown("🔹 :blue[**Step 2.3**] Matrix:", help="Only PWM generated with our tools are allowed")
                matrix_text = st.text_area("🔹 :blue[**Step 2.3**] Matrix:",
                                           value="A [ 20.0 0.0 0.0 0.0 0.0 0.0 0.0 100.0 0.0 60.0 20.0 ]\nT [ 60.0 20.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]\nG [ 0.0 20.0 100.0 0.0 0.0 100.0 100.0 0.0 100.0 40.0 0.0 ]\nC [ 20.0 60.0 0.0 100.0 100.0 0.0 0.0 0.0 0.0 0.0 80.0 ]",
                                           label_visibility='collapsed')

                pwm_rows = kmatrix_text.strip().split('\n')
                pwm = [list(map(str, row.split())) for row in pwm_rows]

                try:
                    has_uniform_column_length(pwm)
                    error_input_im = True
                except Exception as e:
                    error_input_im = False
                    st.error(e)
        else:
            with REcol1:
                st.markdown("🔹 :blue[**Step 2.3**] Sequences:",
                            help='Put FASTA sequences. Same sequence length required ⚠')
                fasta_text = st.text_area("🔹 :blue[**Step 2.3**] Sequences:",
                                          value=">seq1\nCTGCCGGAGGA\n>seq2\nAGGCCGGAGGC\n>seq3\nTCGCCGGAGAC\n>seq4\nCCGCCGGAGCG\n>seq5\nAGGCCGGATCG",
                                          label_visibility='collapsed')
                fasta_text = fasta_text.upper()
            isUIPAC = True

            try:
                matrix_text, buffer = im(fasta_text)
                error_input_im = True
            except Exception as e:
                error_input_im = False
                st.error(e)

    else:
        with REcol1:
            st.markdown("🔹 :blue[**Step 2.3**] Responsive element:", help="IUPAC authorized")
            IUPAC = st.text_input("🔹 :blue[**Step 2.3**] Responsive element (IUPAC authorized):", value="GGGRNYYYCC",
                                  label_visibility='collapsed')
            IUPAC = IUPAC.upper()

        IUPAC_code = ['A', 'T', 'G', 'C', 'R', 'Y', 'M', 'K', 'W', 'S', 'B', 'D', 'H', 'V', 'N']

        if all(char in IUPAC_code for char in IUPAC):
            isUIPAC = True

            sequences = generate_iupac_variants(IUPAC)
            fasta_text = ""
            for i, seq in enumerate(sequences):
                fasta_text += f">seq{i + 1}\n{seq}\n"

            try:
                matrix_text, buffer = im(fasta_text)
                error_input_im = True
            except Exception as e:
                error_input_im = False
                st.error(e)

        else:
            isUIPAC = False

    # TSS entry
    BSFcol1, BSFcol2, BSFcol3 = st.columns([2, 2, 2], gap="medium")
    with BSFcol1:
        st.markdown("🔹 :blue[**Step 2.4**] Transcription Start Site (TSS)/gene end at (in bp):",
                    help="Distance of TSS and gene end from begin of sequences. If you use Step 1, it is positive value of upstream")
        if 'upstream' not in st.session_state:
            entry_tis = st.number_input("🔹 :blue[**Step 2.4**] Transcription Start Site (TSS)/gene end at (in bp):",
                                        -10000, 10000, 0, label_visibility="collapsed")
        else:
            entry_tis = st.number_input("🔹 :blue[**Step 2.4**] Transcription Start Site (TSS)/gene end at (in bp):",
                                        -10000, 10000, st.session_state['upstream'], label_visibility="collapsed")

    # Threshold pvalue

    with BSFcol2:
        st.markdown("🔹 :blue[**Step 2.5**] Relative Score threshold")
        auto_thre = st.checkbox("Automatic threshold", value=True)
        if auto_thre:
            threshold_entry = 0
        else:
            threshold_entry = st.slider("🔹 :blue[**Step 2.5**] Relative Score threshold", 0.5, 1.0, 0.85, step=0.05,
                                    label_visibility="collapsed")
    with BSFcol3:
        st.markdown("🔹 :blue[**_Experimental_**] Calcul _p-value_", help='Experimental, take more times. 10 sequences max.')
        if total_promoter > 10:
            calc_pvalue_stop = True
            st.warning('⚠️_p-value_ not allowed. 10 sequences max. Insufficient server resource.')
        else:
            calc_pvalue_stop = False
        calc_pvalue = st.checkbox('_p-value_', disabled=calc_pvalue_stop)

    # Run Responsive Elements finder
    tis_value = int(entry_tis)
    threshold = float(threshold_entry)
    if jaspar == 'JASPAR_ID':
        sequence_consensus_input = entry_sequence
    else:
        if not isUIPAC:
            st.error("Please use IUPAC code for Responsive Elements")
            button = True
        elif not error_input_im:
            button = True
        elif isfasta:
            st.error("Please use only A, T, G, C, N in your sequence")
            button = True
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
            button = False
    st.markdown("")
    with st.form("runBSF"):
        if st.form_submit_button("🔹 :blue[**Step 2.6**] Click here to find motif in your sequences 🔎 🧬", use_container_width=True, disabled=button):
            matrices = transform_matrix(matrix)
            table2 = search_sequence(threshold, tis_value, promoters, matrices, total_promoter_region_length)
            st.session_state['table2'] = table2

    st.divider()
    if 'table2' in st.session_state:
        if len(st.session_state['table2']) > 1:
            current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            st.subheader(':blue[Results]')

            df = pd.DataFrame(st.session_state['table2'][1:], columns=st.session_state['table2'][0])
            st.session_state['df'] = df

            st.markdown('**Table**')
            tablecol1, tablecol2 = st.columns([0.75, 0.25])
            with tablecol1:
                st.dataframe(df, hide_index=True)
                excel_file = io.BytesIO()
                df.to_excel(excel_file, index=False, sheet_name='Sheet1')
                excel_file.seek(0)

            with tablecol2:
                st.success(f"Finding responsive elements done !")

            st.markdown("")
            st.markdown('**Graph**',
                        help='Zoom +/- with the mouse wheel. Drag while pressing the mouse to move the graph. Selection of a group by clicking on a point of the graph (double click de-selection). Double-click on a point to reset the zoom and the moving of graph.')

            result_table_output(df)

            with tablecol2:
                st.download_button("💾 Download table (.xlsx)", excel_file,
                                   file_name=f'Results_TFinder_{current_date_time}.xlsx',
                                   mime="application/vnd.ms-excel", key='download-excel')
                email_receiver = st.text_input('Send results by email ✉',
                                               value='Send results by email ✉',
                                               label_visibility="collapsed")
                if st.button("Send ✉"):
                    if jaspar == 'PWM':
                        if matrix_type == 'With PWM':
                            body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nPosition Weight Matrix:\n{matrix_text}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                        if matrix_type == 'With FASTA sequences':
                            body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nResponsive Elements:\n{fasta_text}\n\nPosition Weight Matrix:\n{matrix_text}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                    elif jaspar == 'JASPAR_ID':
                        body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nJASPAR_ID: {sequence_consensus_input} | Transcription Factor name: {TF_name}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                    else:
                        body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nResponsive Elements:\n{IUPAC}\n\nPosition Weight Matrix:\n{matrix_text}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                    email(excel_file, txt_output, email_receiver, body)
        else:
            st.error(f"No consensus sequence found with the specified threshold")