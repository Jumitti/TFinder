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

headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
}


class NCBIdna:
    def __init__(self,
                 gene_id,
                 prom_term,
                 upstream,
                 downstream,
                 species=None, gr="Current", all_slice_forms=None):
        self.gene_id = gene_id
        self.prom_term = prom_term if prom_term is not None else None
        self.upstream = upstream if upstream is not None else None
        self.downstream = downstream if downstream is not None else None
        self.species = species if species is not None else None
        self.gr = gr if gr is not None else "Current"
        self.all_slice_forms = True if all_slice_forms is True else False

    @staticmethod
    def XMNM_to_gene_ID(variant):
        uids = f"https://www.ncbi.nlm.nih.gov/nuccore/{variant}"

        response = requests.get(uids)

        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')

            pattern = r"list_uids=(\d+)"
            matches = re.search(pattern, str(soup))

            if matches:
                entrez_id = matches.group(1)
            else:
                entrez_id = "UIDs not founds"
        else:
            entrez_id = "Error during process of retrieving UIDs"

        return entrez_id

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
        if self.gene_id.startswith('XM_') or self.gene_id.startswith('NM_') or self.gene_id.startswith(
                'XR_') or self.gene_id.startswith('NR_'):
            entrez_id = NCBIdna.XMNM_to_gene_ID(self.gene_id)
            if entrez_id == 'UIDs not founds' or entrez_id == 'Error during process of retrieving UIDs':
                result_promoter = f'Please verify {self.gene_id} variant'
                return result_promoter
            else:
                variant, gene_name, title, chraccver, chrstart, chrstop, species_API = NCBIdna.get_variant_info(
                    entrez_id,
                    self.gene_id)
        else:
            if self.gene_id.isdigit():
                entrez_id = self.gene_id

            else:
                entrez_id, message = NCBIdna.convert_gene_to_entrez_id(self.gene_id, self.species)
                if entrez_id == "Error 200":
                    return entrez_id, message

            if not self.all_slice_forms:
                variant, gene_name, title, chraccver, chrstart, chrstop, species_API, message = NCBIdna.get_gene_info(
                    entrez_id, self.gr, gene_name_error=self.gene_id)
                if variant == "Error 200":
                    return variant, message

            elif self.all_slice_forms:
                all_variants, message = NCBIdna.all_variant(entrez_id)
                if all_variants == "Error 200":
                    return all_variants, message

        prom_term = self.prom_term.lower()
        if prom_term not in ['promoter', 'terminator']:
            result_promoter = f"'{self.prom_term}' not valid. Please use 'Promoter' or 'Terminator'."
            return result_promoter, "OK"

        if isinstance(self.upstream, int) and isinstance(self.downstream, int):
            upstream = int(self.upstream)
            downstream = int(self.downstream)
        else:
            result_window = f'Upstream {self.upstream} and Downstream {self.downstream} must be integer'
            return result_window, "OK"

        if not self.all_slice_forms or self.all_slice_forms and self.gene_id.startswith(
                'XM_') or self.gene_id.startswith('NM_') or self.gene_id.startswith('XR_') or self.gene_id.startswith('NR_'):
            dna_sequence = NCBIdna.get_dna_sequence(prom_term, upstream, downstream, chraccver, chrstart, chrstop)

            if prom_term == 'promoter':
                dna_sequence = f">{variant} {gene_name} | {title} {chraccver} | {self.prom_term} | TSS (on chromosome): {chrstart + 1} | TSS (on sequence): {self.upstream}\n{dna_sequence}\n"
            else:
                dna_sequence = f">{variant} {gene_name} | {title} {chraccver} | {self.prom_term} | Gene end (on chromosome): {chrstop} | Gene end (on sequence): {self.upstream}\n{dna_sequence}\n"

            return dna_sequence, "OK"

        elif self.all_slice_forms:
            result_compil = []
            for variant, gene_name, chraccver, chrstart, chrstop, species_API in all_variants:
                dna_sequence = NCBIdna.get_dna_sequence(prom_term, upstream, downstream, chraccver, chrstart, chrstop)
                if prom_term == 'promoter':
                    results = f">{variant} {gene_name} | {species_API} | {chraccver} | {self.prom_term} | TSS (on chromosome): {chrstart} | TSS (on sequence): {self.upstream}\n{dna_sequence}\n"
                else:
                    results = f">{variant} {gene_name} | {species_API} | {chraccver} | {self.prom_term} | Gene end (on chromosome): {chrstop} | Gene end (on sequence): {self.upstream}\n{dna_sequence}\n"

                result_compil.append(results)

            result_output = "\n".join(result_compil)
            return result_output, "OK"

    @staticmethod
    # Convert gene to ENTREZ_GENE_ID
    def convert_gene_to_entrez_id(gene_name, species):
        if gene_name.isdigit():
            gene_name = 'Already gene ID'
            return gene_name  # Already an ENTREZ_GENE_ID

        while True:
            # Request for ENTREZ_GENE_ID
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term="{gene_name}"[Gene%20Name]+AND+{species}[Organism]&retmode=json&rettype=xml'
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                response_data = response.json()

                if response_data['esearchresult']['count'] == '0':
                    return "Error 200", f"Please verify if {gene_name} exist for {species}"

                else:
                    gene_id = response_data['esearchresult']['idlist'][0]
                    return gene_id, f"Response 200: ID found for {species} {gene_name}: {gene_id}"

            elif response.status_code == 429:
                time.sleep(random.uniform(0.25, 0.5))
            else:
                time.sleep(random.uniform(0.25, 0.5))

    @staticmethod
    # Get gene information
    def get_gene_info(entrez_id, gr="Current", from_id=True, gene_name_error=None):

        while True:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={entrez_id}&retmode=json&rettype=xml"
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                response_data = response.json()
                try:
                    gene_info = response_data['result'][str(entrez_id)]
                    gene_name = gene_info['name']
                    title, chraccver, chrstart, chrstop = NCBIdna.extract_genomic_info(entrez_id, response_data, gr)
                    species_API = gene_info['organism']['scientificname']

                    if from_id:
                        variant = NCBIdna.all_variant(entrez_id, from_id=True)

                    return variant, gene_name, title, chraccver, chrstart, chrstop, species_API, (
                        f"Response 200: Info for {entrez_id} {variant} {gene_name} retrieved."
                        f"Entrez_ID: {entrez_id} | Gene name: {gene_name} | Genome Assembly: {title} | ChrAccVer: {chraccver}"
                        f"ChrStart/ChrStop: {chrstart}/{chrstop} | Species: {species_API}")
                except Exception as e:
                    return "Error 200", None, None, None, None, None, None, f"Info for {gene_name_error} {entrez_id} not found."

            elif response.status_code == 429:
                time.sleep(random.uniform(0.25, 0.5))
            else:
                time.sleep(random.uniform(0.25, 0.5))

    @staticmethod
    def extract_genomic_info(gene_id, gene_info, gr):
        if gene_info and 'result' in gene_info and gene_id in gene_info['result']:
            accession_dict = {}
            gene_details = gene_info['result'][gene_id]

            time.sleep(1)

            # Extraction depuis la section locationhist
            location_hist = gene_details.get('locationhist', [])
            if len(location_hist) == 0:
                location_hist = gene_details.get('genomicinfo', [])
            for loc in location_hist:
                nc_accver = loc.get('chraccver')  # Extrait le NC_XXXXX.YY
                chrstart = loc.get('chrstart')
                chrstop = loc.get('chrstop')

                if nc_accver:
                    base_accession = nc_accver
                    if base_accession not in accession_dict:
                        accession_dict[base_accession] = (chrstart, chrstop)
                    else:
                        # Conserver la première occurrence des coordonnées
                        existing_start, existing_stop = accession_dict[base_accession]
                        accession_dict[base_accession] = (min(existing_start, chrstart), max(existing_stop, chrstop))

            nc_dict = accession_dict

            # Affichage des NC_ et leurs informations de position
            # for base_accver, (chrstart, chrstop) in nc_dict.items():
            #     print(f"{base_accver}: chrstart={chrstart}, chrstop={chrstop}")

            # Filtrer pour garder uniquement les entrées qui commencent par "NC_"
            nc_dict = {base_accver: (chrstart, chrstop) for base_accver, (chrstart, chrstop) in nc_dict.items() if
                       base_accver.startswith("NC_")}

            # Si le dictionnaire n'est pas vide, récupérez la base avant le point du premier élément
            if nc_dict:
                first_base = next(iter(nc_dict)).split('.')[0]  # Récupérer la base avant le point de la première clé

                # Filtrer le dictionnaire pour ne garder que les éléments avec la même base
                nc_dict = {base_accver: (chrstart, chrstop) for base_accver, (chrstart, chrstop) in nc_dict.items() if
                           base_accver.split('.')[0] == first_base}

            # Votre code existant pour afficher le dictionnaire filtré
            # print("\nDictionnaire filtré :")
            # for base_accver, (chrstart, chrstop) in nc_dict.items():
            #     print(f"{base_accver}: chrstart={chrstart}, chrstop={chrstop}")

            # Trouver l'entrée avec le chiffre le plus grand et le plus petit après le point
            max_version = -1
            max_accver = None
            max_coords = None  # Pour stocker les coordonnées associées à la version maximale
            min_version = float('inf')  # Initialiser à l'infini pour trouver la plus petite version
            min_accver = None
            min_coords = None  # Pour stocker les coordonnées associées à la version minimale

            for base_accver in nc_dict.keys():
                version = int(base_accver.split('.')[1])  # Extraire la version après le point

                # Mise à jour pour la version maximale
                if version > max_version:
                    max_version = version
                    max_accver = base_accver
                    max_coords = nc_dict[base_accver]  # Stocker les coordonnées

                # Mise à jour pour la version minimale
                if version < min_version:
                    min_version = version
                    min_accver = base_accver
                    min_coords = nc_dict[base_accver]  # Stocker les coordonnées

            # Afficher les résultats
            # if max_accver:
            #     print(f"\nL'accès avec la version la plus élevée est : {max_accver} avec la version {max_version}.")
            #     print(f"Coordonnées : chrstart={max_coords[0]}, chrstop={max_coords[1]}")
            # else:
            #     print("\nAucun accès trouvé.")
            #
            # if min_accver:
            #     print(f"L'accès avec la version la plus basse est : {min_accver} avec la version {min_version}.")
            #     print(f"Coordonnées : chrstart={min_coords[0]}, chrstop={min_coords[1]}")
            # else:
            #     print("Aucun accès trouvé.")

            if gr != "Current":
                title = NCBIdna.fetch_nc_info(min_accver)
                return title, min_accver, min_coords[0], min_coords[1]
            else:
                title = NCBIdna.fetch_nc_info(max_accver)
                return title, max_accver, max_coords[0], max_coords[1]

    @staticmethod
    def fetch_nc_info(nc_accver):
        while True:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={nc_accver}&retmode=json"
            response = requests.get(url)
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

    @staticmethod
    # Get gene information
    def get_variant_info(entrez_id, variant):

        variant = variant.split(".")[0]

        while True:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={entrez_id}&retmode=xml"
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                root = ET.fromstring(response.text)

                chromosome = ""
                found_variant = False
                k_found = False
                start_coords = []
                end_coords = []
                orientation = ""
                break

            elif response.status_code == 429:
                time.sleep(random.uniform(0.25, 0.5))
            else:
                time.sleep(random.uniform(0.25, 0.5))

        if response.status_code == 200:
            for elem in root.iter('Gene-commentary_accession'):
                if elem.text.startswith('NC_'):
                    chromosome = elem.text
                    break

            for elem in root.iter():
                if elem.tag == "Gene-commentary_accession" and elem.text != variant:
                    if elem.text == chromosome:
                        k_found = True
                    elif len(start_coords) < 2 and len(end_coords) < 2:
                        continue
                    else:
                        break

                if k_found and elem.tag == "Gene-commentary_accession":
                    found_variant = True if elem.text == variant else False
                elif k_found and found_variant and elem.tag == "Seq-interval_from":
                    start_coords.append(elem.text)
                elif k_found and found_variant and elem.tag == "Seq-interval_to":
                    end_coords.append(elem.text)
                elif k_found and found_variant and elem.tag == "Na-strand" and orientation == "":
                    orientation += elem.attrib.get("value")

                elif elem.tag == "Org-ref_taxname":
                    species_API = elem.text

                elif elem.tag == 'Gene-ref_locus':
                    gene_name = elem.text

            if orientation != "minus":
                chrstart = int(start_coords[0]) + 1
                chrstop = int(end_coords[-1]) + 1
            else:
                chrstart = int(end_coords[0]) + 1
                chrstop = int(start_coords[-1]) + 1

            while True:
                url2 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={entrez_id}&retmode=json&rettype=xml"
                response = requests.get(url2, headers=headers)
                if response.status_code == 200:
                    response_data = response.json()
                    if 'result' in response_data and str(entrez_id) in response_data['result']:
                        gene_info = response_data['result'][str(entrez_id)]
                        if 'chraccver' in str(gene_info):
                            chraccver = gene_info['genomicinfo'][0]['chraccver']
                            title = NCBIdna.fetch_nc_info(chraccver)
                    else:
                        gene_name = 'Bad ID'

                    return variant, gene_name, title, chraccver, chrstart, chrstop, species_API

                elif response.status_code == 429:
                    time.sleep(random.uniform(0.25, 0.5))
                else:
                    time.sleep(random.uniform(0.25, 0.5))

    @staticmethod
    # Get gene information
    def all_variant(entrez_id, from_id=False):

        if not from_id:
            while True:
                url2 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={entrez_id}&retmode=json&rettype=xml"
                response = requests.get(url2, headers=headers)

                if response.status_code == 200:
                    response_data = response.json()
                    try:
                        gene_info = response_data['result'][str(entrez_id)]
                        chraccver = gene_info['genomicinfo'][0]['chraccver']
                        break

                    except Exception as e:
                        all_variants = [("Error 200", None, None, None, None, None)]
                        return "Error 200", f"Transcript not found(s) for {entrez_id}."

                elif response.status_code == 429:
                    time.sleep(random.uniform(0.25, 0.5))
                else:
                    time.sleep(random.uniform(0.25, 0.5))

        while True:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={entrez_id}&retmode=xml"
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                root = ET.fromstring(response.text)

                tv = []
                variants = []
                gene_name = []
                species_API = []
                chromosome = ""

                all_variants = []

                for elem in root.iter():
                    if elem.tag == "Gene-commentary_label":
                        if elem.text.startswith('transcript variant'):
                            if elem.text not in tv:
                                tv.append(elem.text)
                    if elem.tag == "Gene-commentary_accession":
                        if elem.text.startswith('NM_') or elem.text.startswith('XM_') or elem.text.startswith(
                                'NR_') or elem.text.startswith('XR_'):
                            if elem.text not in variants:
                                variants.append(elem.text)

                    elif elem.tag == "Org-ref_taxname":
                        species_API = elem.text

                    elif elem.tag == 'Gene-ref_locus':
                        gene_name = elem.text

                if not from_id:
                    for elem in root.iter('Gene-commentary_accession'):
                        if elem.text.startswith('NC_'):
                            chromosome = elem.text
                            break

                    for variant in variants:
                        start_coords = []
                        end_coords = []
                        found_variant = False
                        k_found = False
                        orientation = ""
                        for elem in root.iter():
                            if elem.tag == "Gene-commentary_accession" and elem.text != variant:
                                if elem.text == chromosome:
                                    k_found = True
                                elif len(start_coords) < 2 and len(end_coords) < 2:
                                    continue
                                else:
                                    break

                            if k_found and elem.tag == "Gene-commentary_accession":
                                found_variant = True if elem.text == variant else False
                            elif k_found and found_variant and elem.tag == "Seq-interval_from":
                                start_coords.append(elem.text)
                            elif k_found and found_variant and elem.tag == "Seq-interval_to":
                                end_coords.append(elem.text)
                            elif k_found and found_variant and elem.tag == "Na-strand" and orientation == "":
                                orientation += elem.attrib.get("value")

                            elif elem.tag == "Org-ref_taxname":
                                species_API = elem.text

                            elif elem.tag == 'Gene-ref_locus':
                                gene_name = elem.text

                        if orientation != "minus":
                            chrstart = int(start_coords[0]) + 1
                            chrstop = int(end_coords[-1]) + 1
                        else:
                            chrstart = int(end_coords[0]) + 1
                            chrstop = int(start_coords[-1]) + 1

                        all_variants.append((variant, gene_name, chraccver, chrstart, chrstop, species_API))

                    if len(all_variants) > 0:
                        return all_variants, f"Transcript(s) found(s) for {entrez_id}: {all_variants}"
                    else:
                        gene_name, title, chraccver, chrstart, chrstop, species_API, message = NCBIdna.get_gene_info(
                            entrez_id, from_id=False)
                        if gene_name == "Error 200":
                            all_variants.append(("Error 200", None, None, None, None, None))
                            return "Error 200", f"Transcript not found(s) for {entrez_id}."
                        else:
                            all_variants.append((gene_name, gene_name, chraccver, chrstart, chrstop, species_API))
                            return all_variants, f"Transcript(s) found(s) for {entrez_id}: {all_variants}"

                elif from_id:
                    if len(tv) > 0:
                        associations = dict(zip(tv, variants))
                        return associations["transcript variant 1"]
                    else:
                        return variants[0]

            elif response.status_code == 429:
                time.sleep(random.uniform(0.25, 0.5))
            else:
                time.sleep(random.uniform(0.25, 0.5))

    @staticmethod
    # Get DNA sequence
    def get_dna_sequence(prom_term, upstream, downstream, chraccver, chrstart, chrstop):
        # Determine sens of gene + coordinate for upstream and downstream
        if chrstop > chrstart:
            start = (chrstart if prom_term.lower() == 'promoter' else chrstop) - upstream + 1
            end = (chrstart if prom_term.lower() == 'promoter' else chrstop) + downstream
        else:
            start = (chrstart if prom_term.lower() == 'promoter' else chrstop) + upstream + 1
            end = (chrstart if prom_term.lower() == 'promoter' else chrstop) - downstream + 2

        # Request for DNA sequence
        while True:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&from={start}&to={end}&rettype=fasta&retmode=text"
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                # Extraction of DNA sequence
                dna_sequence = response.text.split('\n', 1)[1].replace('\n', '')
                if chrstop > chrstart:
                    sequence = dna_sequence
                else:
                    sequence = NCBIdna.reverse_complement(dna_sequence)

                return sequence

            elif response.status_code == 429:
                time.sleep(random.uniform(0.25, 0.5))
            else:
                time.sleep(random.uniform(0.25, 0.5))

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
                '+ f': matrix,
                '+ r': reversed_matrix,
                '- f': complement_matrix,
                '- r': reversed_complement_matrix
            }
        else:
            return {
                '+ f': matrix,
                '- r': reversed_complement_matrix
            }

    @staticmethod
    # Generate random sequences for p_value
    def generate_ranseq(probabilities, seq_length, progress_bar):
        motif_length = seq_length
        random_sequences = []

        for _ in range(1000000):
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
                if i < len(base_score):
                    score += base_score[i]
                else:
                    score += 0
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
    def individual_motif_finder(dna_sequences, threshold, matrix, progress_bar=None, calc_pvalue=None,
                                tss_ge_distance=None, alldirection=None):
        if calc_pvalue is not None:
            if calc_pvalue not in ["ATGCPreset", "ATGCProportion"]:
                raise ValueError("Use 'ATGCPreset' or 'ATGCProportion'")

        individual_motif_occurrences = []

        matrices = IMO.transform_matrix(matrix, alldirection)

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

                random_sequences = IMO.generate_ranseq(probabilities, seq_length, progress_bar)

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
                            tis_position = position - tss_ge_distance - 1

                        if normalized_score >= threshold:

                            if calc_pvalue is not None:
                                p_value = (random_scores >= normalized_score).sum() / len(random_scores)

                            row = [position]
                            if tss_ge_distance is not None:
                                row.append(tis_position)
                            row += [sequence_with_context,
                                    "{:.6f}".format(normalized_score).ljust(12)]
                            if calc_pvalue is not None:
                                row.append("{:.6e}".format(p_value))
                            row += [strand, direction, name, species, region]
                            individual_motif_occurrences.append(row)

        if len(individual_motif_occurrences) > 0:
            if tss_ge_distance is not None and calc_pvalue is not None:
                individual_motif_occurrences.sort(key=lambda x: (-float(x[3]), -float(x[4])))
            elif calc_pvalue is not None:
                individual_motif_occurrences.sort(key=lambda x: (-float(x[2]), -float(x[3])))
            elif tss_ge_distance is not None:
                individual_motif_occurrences.sort(key=lambda x: (-float(x[3])))
            else:
                individual_motif_occurrences.sort(key=lambda x: (-float(x[2])))
            header = ["Position"]
            if tss_ge_distance is not None:
                header.append("Rel Position")
            header += ["Sequence", "Rel Score"]
            if calc_pvalue is not None:
                header.append("p-value")
            header += ["Strand", "Direction", "Gene", "Species", "Region"]
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
        if max_variant_allowed is not None:
            if total_variants > max_variant_allowed:
                sequence = f'Too many variants. Use - or . instead N. Limit: {max_variant_allowed} | Total variants : {total_variants}'
                return sequence

        if progress_bar is True:
            pbar = tqdm(total=total_variants, desc='Generate variants from IUPAC code...', mininterval=0.1)
        sequences = [sequence]
        for i, base in enumerate(sequence):
            if base.upper() in iupac_codes:
                new_sequences = []
                for seq in sequences:
                    for alternative in iupac_codes[base.upper()]:
                        new_sequence = seq[:i] + alternative + seq[i + 1:]
                        new_sequences.append(new_sequence)
                        if progress_bar is True:
                            pbar.update(1)
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
        logo.style_spines(visible=True)
        logo.style_xticks(fmt='%d')
        logo.ax.set_ylabel("Bits")
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
            raise Exception(f"You forget to input something :)")

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
