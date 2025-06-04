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

import random
import re
import time
import xml.etree.ElementTree as ET

import altair as alt
import math
import logomaker
import numpy as np
import pandas as pd
import requests
import streamlit as st
from bs4 import BeautifulSoup
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
from Bio import motifs
from Bio.motifs.jaspar import calculate_pseudocounts


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
}


class NCBIdna:
    def __init__(self, gene_id, species=None, seq_type="mrna", upstream=2000, downstream=2000, genome_version="Current",
                 all_slice_forms=None, max_attempts=25):
        self.gene_id = gene_id
        self.species = species if species is not None else "human"
        self.seq_type = seq_type if seq_type is not None else "mrna"
        self.upstream = upstream if upstream is not None and seq_type in ["promoter", "terminator"] else None
        self.downstream = downstream if downstream is not None and seq_type in ["promoter", "terminator"] else None
        self.genome_version = genome_version if genome_version is not None else "current"
        self.all_slice_forms = True if all_slice_forms is True else False
        self.max_attempts = max_attempts

    @staticmethod
    def XMNM_to_gene_ID(variant, max_attempts=25):
        global headers

        attempt = 0

        while attempt < max_attempts:
            attempt += 1

            uids = f"https://www.ncbi.nlm.nih.gov/nuccore/{variant}"

            response = requests.get(uids, headers=headers)

            if response.status_code == 200:
                soup = BeautifulSoup(response.text, 'html.parser')

                pattern = r"list_uids=(\d+)"
                matches = re.search(pattern, str(soup))

                if matches:
                    entrez_id = matches.group(1)
                else:
                    entrez_id = "UIDs not founds"

                return entrez_id
            else:
                print(bcolors.FAIL + "Error during process of retrieving UIDs" + bcolors.ENDC)

        return "UIDs not founds"

    # Sequence extractor
    def find_sequences(self):
        time.sleep(1)
        if self.gene_id.upper().startswith(('XM_', 'NM_', 'XR_', 'NR_', 'YP_')):
            if '.' in self.gene_id:
                self.gene_id = self.gene_id.split('.')[0]
            entrez_id = NCBIdna.XMNM_to_gene_ID(self.gene_id, self.max_attempts)
            if entrez_id == 'UIDs not founds':
                result_promoter = f'Please verify {self.gene_id} variant'
                return result_promoter, result_promoter
        else:
            if self.gene_id.isdigit():
                entrez_id = self.gene_id

            else:
                entrez_id, message = NCBIdna.convert_gene_to_entrez_id(self.gene_id, self.species, self.max_attempts)
                if entrez_id == "Error":
                    return entrez_id, message

        all_variants, message = NCBIdna.all_variant(entrez_id, self.genome_version, self.all_slice_forms,
                                                    self.gene_id if
                                                    self.gene_id.upper().startswith(('XM_', 'NM_', 'XR_', 'NR_', 'YP_')) else None,
                                                    self.max_attempts)
        if "Error" not in all_variants:
            for nm_id, data in all_variants.items():
                exon_coords = data.get('exon_coords')
                data['upstream'] = self.upstream
                data['seq_type'] = self.seq_type
                sequence, message_seq = NCBIdna.get_dna_sequence(data.get("entrez_id"), data.get("chraccver"), exon_coords[0][0],
                                                    exon_coords[-1][1], self.seq_type, self.upstream, self.downstream,
                                                    self.max_attempts)
                if sequence != "Error":
                    if self.seq_type in ['mrna']:
                        if self.seq_type == 'mrna':
                            data['sequence'] = "".join(
                                sequence[start:end + 1] for start, end in data["normalized_exon_coords"])
                    else:
                        data['sequence'] = sequence

            return all_variants, "OK"
        else:
            return all_variants, f"Error: {message}"

    @staticmethod
    # Convert gene to ENTREZ_GENE_ID
    def convert_gene_to_entrez_id(gene_name, species, max_attempts=25):
        global headers

        attempt = 0

        while attempt < max_attempts:
            attempt += 1
            # Request for ENTREZ_GENE_ID
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term="{gene_name}"[Gene%20Name]+AND+{species}[Organism]&retmode=json&rettype=xml'
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                response_data = response.json()

                if 'count' in response_data.get('esearchresult', {}):
                    if response_data['esearchresult']['count'] == '0':
                        print(bcolors.WARNING + f"Please verify if {gene_name} exist for {species}" + bcolors.ENDC)
                        return "Error", f"Please verify if {gene_name} exist for {species}"

                    else:
                        gene_id = response_data['esearchresult']['idlist'][0]
                        print(
                            bcolors.OKGREEN + f"Response 200: ID found for {species} {gene_name}: {gene_id}" + bcolors.ENDC)
                        return gene_id, bcolors.OKGREEN + f"Response 200: ID found for {species} {gene_name}: {gene_id}" + bcolors.ENDC
                else:
                    print(
                        bcolors.FAIL + f"Error 200: Issues for {species} {gene_name}, try again: {response.text}" + bcolors.ENDC)
                    time.sleep(random.uniform(0.25, 0.5))

            elif response.status_code == 429:
                print(
                    bcolors.FAIL + f"Error 429: API rate limit exceeded during get ID of {species} {gene_name}, try again: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))
            else:
                print(bcolors.FAIL + f"Error {response.status_code}, try again: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))

        return "Error", f"{response.status_code} for {gene_name} {species}"

    @staticmethod
    # Get gene information
    def all_variant(entrez_id, genome_version="current", all_slice_forms=False, specific_transcript=None, max_attempts=25):
        global headers

        attempt = 0
        while attempt < max_attempts:
            attempt += 1

            url2 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={entrez_id}&retmode=json&rettype=xml"
            response = requests.get(url2, headers=headers)

            if response.status_code == 200:
                response_data = response.json()
                try:
                    gene_info = response_data['result'][str(entrez_id)]
                    species_API = gene_info['organism']['scientificname']
                    gene_name = gene_info['name']
                    title, chrloc, chraccver, coords = NCBIdna.extract_genomic_info(entrez_id, response_data,
                                                                                    genome_version, species_API, max_attempts)
                    chrstart = coords[0]
                    chrstop = coords[1]
                    print(
                        bcolors.OKGREEN + f"Response 200: Chromosome {chrloc} {chraccver} found for {entrez_id}: {response.text}" + bcolors.ENDC)
                    break

                except Exception as e:
                    print(
                        bcolors.WARNING + f"Response 200: Chromosome not found for {entrez_id}: {response.text} {e} {traceback.print_exc()}" + bcolors.ENDC)
                    print(
                        bcolors.WARNING + f"Response 200: Transcript not found(s) for {entrez_id}." + bcolors.ENDC)
                    return [("Error", None, None, None, None, None)], f"Transcript not found(s) for {entrez_id}."

            elif response.status_code == 429:
                print(
                    bcolors.ENDC + f"Error {response.status_code}: API rate limit exceeded during get chromosome of {entrez_id}: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))
            else:
                print(
                    bcolors.ENDC + f"Error {response.status_code}: Error during get chromosome of {entrez_id}: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))

        if attempt == max_attempts:
            return [("Error", None, None, None, None, None)], f"{response.status_code} {entrez_id}."

        attempt = 0
        while attempt < max_attempts:
            attempt += 1

            all_variants = {}

            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={entrez_id}&retmode=xml"
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                root = ET.fromstring(response.text)

                tv = []
                variants = []
                chromosome = ""

                for elem in root.iter():
                    if elem.tag == "Gene-commentary_label":
                        if elem.text.startswith('transcript variant'):
                            if elem.text not in tv:
                                tv.append(elem.text)
                    if elem.tag == "Gene-commentary_accession":
                        if elem.text.upper().startswith(('XM_', 'NM_', 'XR_', 'NR_', 'YP_')):
                            if elem.text not in variants:
                                variants.append(elem.text)

                    if elem.tag == "Gene-commentary_type":
                        if elem.attrib.get("value") in ["tRNA", "rRNA", "d-segment", "v-segment", "j-segment"]:
                            all_variants[entrez_id] = {
                                'entrez_id': entrez_id,
                                'gene_name': gene_name,
                                'genomic_info': title,
                                'chraccver': chraccver,
                                'strand': "plus" if chrstart < chrstop else "minus",
                                'exon_coords': [(chrstart if chrstart < chrstop else chrstop,
                                                 chrstop if chrstart < chrstop else chrstart)],
                                'normalized_exon_coords': [(0, abs(chrstart - chrstop))],
                                'species': species_API
                            }
                            return all_variants, f"Transcript(s) found(s) for {entrez_id}: {list(all_variants.keys())}"

                for elem in root.iter('Gene-commentary_accession'):
                    if elem.text.startswith(("NC_", "NT_")):
                        chromosome = elem.text
                        break

                def calc_exon(root, variants):

                    for variant in variants:
                        exon_coords = []
                        found_variant = False
                        k_found = False
                        orientation = ""
                        for elem in root.iter():
                            if elem.tag == "Gene-commentary_accession" and elem.text != variant:
                                if elem.text == chromosome:
                                    k_found = True
                                elif len(exon_coords) == 0:
                                    continue
                                else:
                                    break

                            if k_found and elem.tag == "Gene-commentary_accession":
                                found_variant = True if elem.text == variant else False
                            elif k_found and found_variant and elem.tag == "Seq-interval_from":
                                start = int(elem.text)
                            elif k_found and found_variant and elem.tag == "Seq-interval_to":
                                end = int(elem.text)
                                exon_coords.append((start, end))
                            elif k_found and found_variant and elem.tag == "Na-strand" and orientation == "":
                                orientation += elem.attrib.get("value")

                            elif elem.tag == "Org-ref_taxname":
                                species_API = elem.text

                        if exon_coords:
                            if orientation == "minus":
                                exon_coords = [(end, start) for start, end in exon_coords]

                            first_exon_start = exon_coords[0][0]

                            normalized_exon_coords = [
                                (abs(start - first_exon_start), abs(end - first_exon_start)) for start, end in
                                exon_coords]

                            all_variants[variant] = {
                                'entrez_id': entrez_id,
                                'gene_name': gene_name,
                                'genomic_info': title,
                                'chraccver': chraccver,
                                'strand': orientation,
                                'exon_coords': exon_coords,
                                'normalized_exon_coords': normalized_exon_coords,
                                'species': species_API
                            }

                    if len(all_variants) > 0:
                        print(
                            bcolors.OKGREEN + f"Response 200: Transcript(s) found(s) for {entrez_id}: {all_variants}" + bcolors.ENDC)
                        return all_variants, f"Transcript(s) found(s) for {entrez_id}: {list(all_variants.keys())}"
                    else:
                        all_variants["Error"] = {
                            "entrez_id": f"Transcript not found for {entrez_id}.",
                            "gene_name": None,
                            "chraccver": None,
                            "exon_coords": None,
                            "normalized_exon_coords": None,
                            "species": None
                        }
                        print(
                            bcolors.WARNING + f"Error 200: Transcript not found(s) for {entrez_id}." + bcolors.ENDC)
                        return all_variants, f"Error 200: Transcript not found(s) for {entrez_id}."

                if all_slice_forms is True:
                    all_variants, message = calc_exon(root, variants)
                    return all_variants, message

                elif all_slice_forms is False:
                    variant = None

                    if specific_transcript and re.search(r"\.\d+$", specific_transcript):
                        specific_transcript = re.sub(r"\.\d+$", "", specific_transcript)

                    if len(tv) > 0:
                        if specific_transcript and specific_transcript in variants:
                            variant = specific_transcript
                        elif "transcript variant 1" in tv:
                            associations = dict(zip(tv, variants))
                            variant = associations["transcript variant 1"]
                        else:
                            variant = variants[0]

                    elif len(tv) == 0 and len(variants) > 0:
                        variant = variants[0]

                    if variant is not None:
                        all_variants, message = calc_exon(root, [variant])
                        return all_variants, message

            elif response.status_code == 429:
                print(
                    bcolors.FAIL + f"Error {response.status_code}: API rate limit exceeded while searching for {entrez_id} transcripts: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))
            else:
                print(
                    bcolors.FAIL + f"Error {response.status_code}: Error while searching for {entrez_id} transcripts: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))

        return [("Error", None, None, None, None, None)], f"{response.status_code} {entrez_id}."

    @staticmethod
    # Get DNA sequence
    def get_dna_sequence(gene_name, chraccver, chrstart, chrstop, seq_type, upstream=2000, downstream=2000, max_attempts=25):
        global headers

        if seq_type in ['mrna', 'rna']:
            if chrstop > chrstart:
                start = chrstart + 1
                end = chrstop + 1
            else:
                start = chrstop + 1
                end = chrstart + 1
        elif seq_type in ['promoter', 'terminator']:
            if chrstop > chrstart:
                start = (chrstart if seq_type == 'promoter' else chrstop) - upstream + 1
                end = (chrstart if seq_type == 'promoter' else chrstop) + downstream
            else:
                start = (chrstart if seq_type == 'promoter' else chrstop) + upstream + 1
                end = (chrstart if seq_type == 'promoter' else chrstop) - downstream + 2

        attempt = 0
        while attempt < max_attempts:
            attempt += 1

            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&from={start}&to={end}&rettype=fasta&retmode=text"
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                sequence = response.text.split('\n', 1)[1].replace('\n',
                                                                       '')

                if chrstop < chrstart:
                    sequence = NCBIdna.reverse_complement(sequence)

                print(
                    bcolors.OKGREEN + f"Response 200: DNA sequence for {gene_name} extracted: {sequence}" + bcolors.ENDC)
                return sequence, "OK"

            elif response.status_code == 429:
                print(
                    bcolors.FAIL + f"Error 429: API rate limit exceeded for DNA extraction of {gene_name}, try again." + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))
            else:
                print(bcolors.OKGREEN + f"Error {response.status_code}: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))

        return "Error", f"{response.status_code}: during extraction of DNA sequence"

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

    @staticmethod
    def extract_genomic_info(gene_id, gene_info, genome_version, species=None, max_attempts=25):
        if gene_info and 'result' in gene_info and gene_id in gene_info['result']:
            accession_dict = {}
            gene_details = gene_info['result'][gene_id]

            time.sleep(1)

            location_hist = gene_details.get('locationhist', [])
            if len(location_hist) == 0:
                location_hist = gene_details.get('genomicinfo', [])

            for loc in location_hist:

                nc_loc = loc.get('chrloc')
                nc_accver = loc.get('chraccver')
                chrstart = loc.get('chrstart')
                chrstop = loc.get('chrstop')

                if nc_accver:
                    base_accession = nc_accver
                    if base_accession not in accession_dict:
                        accession_dict[base_accession] = (chrstart, chrstop)
                    else:
                        existing_start, existing_stop = accession_dict[base_accession]
                        accession_dict[base_accession] = (min(existing_start, chrstart), max(existing_stop, chrstop))

            nc_dict = accession_dict

            nc_dict = {base_accver: (chrstart, chrstop) for base_accver, (chrstart, chrstop) in nc_dict.items() if
                       base_accver.startswith(("NC_", "NT_"))}

            if species == "Rattus norvegicus" and 'locationhist' in gene_details:
                location_hist = gene_details['locationhist']
                rs_annotations = [
                    loc for loc in location_hist if loc.get('annotationrelease', "").startswith("RS_")
                ]
                rs_annotations.sort(key=lambda x: int(x['annotationrelease'].split('_')[1]))
                if rs_annotations:
                    selected_annotation = rs_annotations[-1] if genome_version == "current" else rs_annotations[0]
                    selected_nc_accver = selected_annotation.get('chraccver')
                    if selected_nc_accver:
                        nc_dict = {
                            selected_nc_accver: (selected_annotation['chrstart'], selected_annotation['chrstop'])
                        }

            if nc_dict:
                first_base = next(iter(nc_dict)).split('.')[0]
                nc_dict = {base_accver: (chrstart, chrstop) for base_accver, (chrstart, chrstop) in nc_dict.items() if
                           base_accver.split('.')[0] == first_base}

            max_version = -1
            max_accver = None
            max_coords = None
            min_version = float('inf')
            min_accver = None
            min_coords = None

            for base_accver in nc_dict.keys():
                version = int(base_accver.split('.')[1])

                if version > max_version:
                    max_version = version
                    max_accver = base_accver
                    max_coords = nc_dict[base_accver]

                if version < min_version:
                    min_version = version
                    min_accver = base_accver
                    min_coords = nc_dict[base_accver]

            if genome_version != "current":
                title = NCBIdna.fetch_nc_info(min_accver, max_attempts)
                return title, nc_loc, min_accver, min_coords
            else:
                title = NCBIdna.fetch_nc_info(max_accver, max_attempts)
                return title, nc_loc, max_accver, max_coords

    @staticmethod
    def fetch_nc_info(nc_accver, max_attempts=25):
        global headers

        attempt = 0
        while attempt < max_attempts:
            attempt += 1

            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={nc_accver}&retmode=json"
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                nc_info = response.json()
                try:
                    uid = nc_info['result']['uids'][0]
                    title = nc_info['result'][uid]['title']
                    return title
                except Exception as e:
                    time.sleep(random.uniform(0.25, 0.5))
            else:
                time.sleep(random.uniform(0.25, 0.5))

        return f"Error {response.status_code} fetch NC info"


class IMO:

    def __init__(self):
        pass

    @staticmethod
    # Extract JASPAR matrix
    def matrix_extraction(jaspar_id):
        url = f"https://jaspar.elixir.no/api/v1/matrix/{jaspar_id}/"
        response = requests.get(url)
        if response.status_code == 200:
            response_data = response.json()
            TF_name = response_data['name']
            TF_species = response_data['species'][0]['name']
            matrix = response_data['pfm']
            weblogo = f"https://jaspar.elixir.no/static/logos/all/svg/{jaspar_id}.svg"
        else:
            TF_name = 'not found'
            TF_species = 'not found'
            matrix = 'not found'
            weblogo = 'not found'
        return TF_name, TF_species, matrix, weblogo

    @staticmethod
    # Transform JASPAR matrix
    def transform_matrix(matrix, alldirection):
        if alldirection is True:
            reversed_matrix = {base: list(reversed(scores)) for base, scores in matrix.items()}
        complement_matrix = {
            'A': matrix['T'],
            'C': matrix['G'],
            'G': matrix['C'],
            'T': matrix['A']
        }
        reversed_complement_matrix = {base: list(reversed(scores)) for base, scores in complement_matrix.items()}

        if alldirection is True:
            return {
                '+ f': {
                    'A': matrix['A'],
                    'C': matrix['C'],
                    'G': matrix['G'],
                    'T': matrix['T']
                },
                '+ r': reversed_matrix,
                '- f': complement_matrix,
                '- r': reversed_complement_matrix
            }
        else:
            return {
                '+ f': {
                    'A': matrix['A'],
                    'C': matrix['C'],
                    'G': matrix['G'],
                    'T': matrix['T']
                },
                '- r': reversed_complement_matrix
            }

    @staticmethod
    # Generate random sequences for p_value
    def generate_ranseq(probabilities, seq_length, progress_bar=None):
        motif_length = seq_length
        random_sequences = []

        for _ in range(1000000):
            random_sequence = IMO.generate_random_sequence(motif_length, probabilities)
            random_sequences.append(random_sequence)
            if progress_bar:
                progress_bar.update(1)

        return random_sequences

    @staticmethod
    # Calculate matrix score
    def calculate_score(sequence, matrix):
        score = sum(
            min(matrix[base][position] for base in matrix) if base not in matrix else matrix[base][position]
            for position, base in enumerate(sequence)
        )
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
        if all(char in DNA_code for char in dna_sequence):
            isfasta = True
        else:
            isfasta = False

        return isfasta

    @staticmethod
    def transform_PWM(matrix, pseudocount=None, bg=None):
        motif = motifs.Motif(counts=matrix)
        pseudocount_auto = calculate_pseudocounts(motif)
        pwm = motif.counts.normalize(pseudocounts=pseudocount_auto if pseudocount is None else pseudocount)
        background_preset = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
        log_odds_matrix = pwm.log_odds(background=bg if bg is not None else background_preset)

        return pseudocount_auto, log_odds_matrix

    @staticmethod
    # Find with JASPAR and manual matrix
    def individual_motif_finder(dna_sequences, threshold, matrix_input, progress_bar=None, calc_pvalue=None,
                                tss_ge_distance=None, alldirection=None, pseudocount=None, bgnf=None):
        if calc_pvalue is not None:
            if calc_pvalue not in ["ATGCPreset", "ATGCProportion"]:
                raise ValueError("Use 'ATGCPreset' or 'ATGCProportion'")

        individual_motif_occurrences = []
        _, log_odds_matrix = IMO.transform_PWM(matrix_input, pseudocount)
        matrices = IMO.transform_matrix(log_odds_matrix, alldirection)

        seq_length = len(matrices['+ f']['A'])

        if calc_pvalue == 'ATGCPreset':
            percentage_a = 0.275
            percentage_t = 0.275
            percentage_g = 0.225
            percentage_c = 0.225

            probabilities = [percentage_a, percentage_c, percentage_g, percentage_t]

            random_sequences = IMO.generate_ranseq(probabilities, seq_length, progress_bar)

        random_scores = {}
        matrix_random_scores = []
        for matrix_name, matrix in matrices.items():
            if calc_pvalue == 'ATGCPreset':
                max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                min_score = sum(min(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                for random_sequence in random_sequences:
                    sequence = random_sequence
                    random_score = IMO.calculate_score(sequence, matrix)
                    normalized_random_score = (random_score - min_score) / (max_score - min_score)
                    matrix_random_scores.append(normalized_random_score)
                    if progress_bar:
                        progress_bar.update(1)
                random_scores = np.array(matrix_random_scores)

        for name, dna_sequence, species, region, strand_seq, tss_ch in dna_sequences:
            count_a = dna_sequence.count('A')
            count_t = dna_sequence.count('T')
            count_g = dna_sequence.count('G')
            count_c = dna_sequence.count('C')

            length_prom = len(dna_sequence)
            percentage_a = count_a / length_prom
            percentage_t = count_t / length_prom
            percentage_g = count_g / length_prom
            percentage_c = count_c / length_prom

            background_freq = {'A': percentage_a, 'C': percentage_c,
                               'G': percentage_g, 'T': percentage_t} if bgnf is None else bgnf

            _, pwm = IMO.transform_PWM(matrix_input, pseudocount, background_freq)
            pwm_weight = IMO.transform_matrix(pwm, alldirection)

            if calc_pvalue == 'ATGCProportion':
                random_sequences = IMO.generate_ranseq([percentage_a, percentage_c, percentage_g, percentage_t],
                                                       seq_length, progress_bar)

                if calc_pvalue == 'ATGCProportion':
                    random_scores = {}

            for matrix_name, matrix in matrices.items():
                if "+" in matrix_name:
                    strand = "+"
                else:
                    strand = "-"
                if "f" in matrix_name:
                    direction = "→"
                else:
                    direction = "←"
                found_positions = []

                # Max score per matrix
                max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                min_score = sum(min(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
                max_score_pwm = sum(
                    max(pwm_weight[matrix_name][base][i] for base in pwm_weight[matrix_name].keys()) for i in
                    range(seq_length))
                min_score_pwm = sum(
                    min(pwm_weight[matrix_name][base][i] for base in pwm_weight[matrix_name].keys()) for i in
                    range(seq_length))

                if calc_pvalue == 'ATGCProportion':
                    matrix_random_scores = []
                    for random_sequence in random_sequences:
                        sequence = random_sequence
                        random_score = IMO.calculate_score(sequence, matrix)
                        normalized_random_score = (random_score - min_score) / (max_score - min_score)
                        matrix_random_scores.append(normalized_random_score)
                        if progress_bar:
                            progress_bar.update(1)

                    random_scores = np.array(matrix_random_scores)

                for i in range(len(dna_sequence) - seq_length + 1):
                    seq = dna_sequence[i:i + seq_length]
                    score = IMO.calculate_score(seq, matrix)
                    normalized_score = (score - min_score) / (max_score - min_score)
                    position = int(i) + 1

                    found_positions.append((position, seq, normalized_score, score))

                    if progress_bar:
                        pb = 1 if region in ["prom.", "term.", "n.d"] else 2
                        progress_bar.update(pb)

                # Sort positions in descending order of score percentage
                found_positions.sort(key=lambda x: x[1], reverse=True)

                # Creating a results table
                if len(found_positions) > 0:
                    if threshold < 0.5:
                        highest_normalized_score = max(
                            [normalized_score for _, _, normalized_score, _, in found_positions])
                        if highest_normalized_score >= 0.6:
                            threshold = highest_normalized_score - 0.10
                        else:
                            threshold = 0.5

                    for position, seq, normalized_score, score in found_positions:
                        if normalized_score >= threshold:

                            weight = IMO.calculate_score(seq, pwm_weight[matrix_name])
                            relative_weight = (weight - min_score_pwm) / (max_score_pwm - min_score_pwm)

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

                            tis_position = None
                            ch_pos = None
                            if region in ["prom.", "term."]:
                                if tss_ge_distance is not None:
                                    tis_position = position - tss_ge_distance - 1

                                    if strand_seq in ["plus", "minus"] and tss_ch != 0:
                                        ch_pos = int(tss_ch) - tis_position if strand_seq == "minus" else int(
                                            tss_ch) + tis_position

                            elif region in ["rna.", "mrna."]:
                                if strand_seq in ["plus", "minus"] and tss_ch != 0:
                                    ch_pos = int(tss_ch) - position if strand_seq == "minus" else int(
                                        tss_ch) + position

                            if calc_pvalue is not None:
                                p_value = (random_scores >= normalized_score).sum() / len(random_scores)

                            row = [position]
                            if tis_position is not None:
                                row.append(tis_position)
                            if ch_pos is not None:
                                row.append(ch_pos)
                            row += [sequence_with_context,
                                    "{:.6f}".format(normalized_score).ljust(12),
                                    "{:.6f}".format(score).ljust(12),
                                    "{:.6f}".format(relative_weight).ljust(12),
                                    "{:.6f}".format(weight).ljust(12)]
                            if calc_pvalue is not None:
                                row.append("{:.6e}".format(p_value))
                            row += [strand, direction, name, species, region]
                            individual_motif_occurrences.append(row)

                if region in ["rna.", "mrna."]:
                    break

        if len(individual_motif_occurrences) > 0:
            header = ["Position"]
            if tis_position is not None:
                header.append("Rel Position")
            if ch_pos is not None:
                header.append("Ch Position")
            header += ["Sequence", "Rel Score", 'Score', 'Rel Score Adj', 'Score Adj']
            if calc_pvalue is not None:
                header.append("p-value")
            header += ["Strand", "Direction", "Gene", "Species", "Region"]
            individual_motif_occurrences.insert(0, header)

            df = pd.DataFrame(individual_motif_occurrences[1:], columns=individual_motif_occurrences[0])
            df = df.sort_values(by="Rel Score", ascending=False)

            return df, True
        else:
            return "No consensus sequence found with the specified threshold.", False

    @staticmethod
    def generate_iupac_variants(sequence):
        iupac_codes_score = {
            "R": 2,
            "Y": 2,
            "M": 2,
            "K": 2,
            "W": 2,
            "S": 2,
            "B": 3,
            "D": 3,
            "H": 3,
            "V": 3,
            "N": 4
        }

        total_variants = 1
        for base in sequence:
            if base.upper() in iupac_codes_score:
                total_variants *= iupac_codes_score[base.upper()]

        iupac_frequencies = {
            "A": {"A": 1.0, "C": 0.0, "G": 0.0, "T": 0.0},
            "C": {"A": 0.0, "C": 1.0, "G": 0.0, "T": 0.0},
            "G": {"A": 0.0, "C": 0.0, "G": 1.0, "T": 0.0},
            "T": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 1.0},
            "R": {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0},
            "Y": {"A": 0.0, "C": 0.5, "G": 0.0, "T": 0.5},
            "M": {"A": 0.5, "C": 0.5, "G": 0.0, "T": 0.0},
            "K": {"A": 0.0, "C": 0.0, "G": 0.5, "T": 0.5},
            "W": {"A": 0.5, "C": 0.0, "G": 0.0, "T": 0.5},
            "S": {"A": 0.0, "C": 0.5, "G": 0.5, "T": 0.0},
            "B": {"A": 0.0, "C": 0.333, "G": 0.333, "T": 0.333},
            "D": {"A": 0.333, "C": 0.0, "G": 0.333, "T": 0.333},
            "H": {"A": 0.333, "C": 0.333, "G": 0.0, "T": 0.333},
            "V": {"A": 0.333, "C": 0.333, "G": 0.333, "T": 0.0},
            "N": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            "-": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0},
            ".": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0}}

        pwm_dict = {'A': [], 'C': [], 'G': [], 'T': []}
        for char in sequence:
            freq = iupac_frequencies.get(char, {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0})

            pwm_dict['A'].append(math.ceil(freq['A'] * total_variants))
            pwm_dict['C'].append(math.ceil(freq['C'] * total_variants))
            pwm_dict['G'].append(math.ceil(freq['G'] * total_variants))
            pwm_dict['T'].append(math.ceil(freq['T'] * total_variants))

        return pwm_dict

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
                pwm[0, i] = counts['A']
                pwm[1, i] = counts['T']
                pwm[2, i] = counts['G']
                pwm[3, i] = counts['C']

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
    def create_web_logo(sequences, type_input):
        if type_input == "PWM":
            matrix = logomaker.alignment_to_matrix(sequences)
        elif type_input == "UIPAC":
            matrix = pd.DataFrame(sequences)
        logo = logomaker.Logo(matrix, color_scheme='classic')
        logo.style_spines(visible=True)
        logo.style_xticks(fmt='%d')
        logo.ax.set_ylabel("Bits")
        return logo

    @staticmethod
    def individual_motif_pwm(individual_motif, type_input="PWM"):
        if type_input == "PWM":
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

                    weblogo = IMO.create_web_logo(sequences, type_input)

                    return matrix, weblogo

            else:
                raise Exception(f"You forget to input something :)")

        elif type_input == "UIPAC":
            weblogo = IMO.create_web_logo(individual_motif, type_input)

            return None, weblogo

    @staticmethod
    def PWM_to_weblogo(pwm_str):
        lines = pwm_str.strip().split('\n')
        pwm_dict = {}
        for line in lines:
            parts = line.strip().split()
            key = parts[0]
            values = [float(x) for x in parts[1:] if x != '[' and x != ']']
            pwm_dict[key] = values

        pwm_df = pd.DataFrame(pwm_dict)

        # Créer le weblogo
        weblogo = logomaker.Logo(pwm_df, color_scheme='classic')

        # Style du weblogo
        weblogo.style_spines(visible=True)
        weblogo.style_xticks(fmt='%d')
        weblogo.ax.set_ylabel("Bits")

        return weblogo

    @staticmethod
    def normalize_matrix(matrix):
        normalized_matrix = {}
        num_rows = len(next(iter(matrix.values())))
        for key in matrix.keys():
            normalized_matrix[key] = [matrix[key][i] / sum(matrix[k][i] for k in matrix) for i in range(num_rows)]
        return normalized_matrix

    @staticmethod
    def generate_sequences(matrix):
        normalized_matrix = IMO.normalize_matrix(matrix)
        keys = list(normalized_matrix.keys())
        sequence_length = len(next(iter(normalized_matrix.values())))

        generated_sequences = []

        def generate_sequence_helper(current_sequence, position):
            if position == sequence_length:
                generated_sequences.append(current_sequence)
                return

            for key in keys:
                if normalized_matrix[key][position] > 0:
                    generate_sequence_helper(current_sequence + key, position + 1)

        generate_sequence_helper("", 0)

        return generated_sequences
