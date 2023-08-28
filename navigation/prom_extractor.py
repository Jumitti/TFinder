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

import datetime
import io
import random
import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

import altair as alt
import logomaker
import numpy as np
import pandas as pd
import requests
import streamlit as st
from stqdm import stqdm

from tfinder import NCBIdna


def prom_extractor_page():
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
                    st.session_state['result_promoter'] = result_promoter
                else:
                    result_promoter.append(
                        f">{gene_name} | {species_API} | {chraccver} | {prom_term} | Gene end (on chromosome): {chrstop} | Gene end (on sequence): {upstream}\n{dna_sequence}\n")
                    st.session_state['result_promoter'] = result_promoter

            return result_promoter

        except Exception as e:
            raise Exception(f"Error retrieving gene information: {entrez_id}")

    st.subheader('🧬 Gene Region Extractor')
    colprom1, colprom2 = st.columns([0.8, 1.2], gap="small")
    # Promoter Finder
    with colprom1:

        result_promoter = []
        upstream_entry = []

        # Gene ID
        st.markdown("🔹 :blue[**Step 1.1**] Gene ID:", help='NCBI gene name and NCBI gene ID allowed')
        gene_id_entry = st.text_area("🔹 :blue[**Step 1.1**] Gene ID:", value="PRKN\n351",
                                     label_visibility='collapsed')
        gene_list = gene_id_entry.strip().split('\n')

        # Verify if gene is available for all species
        if st.button('🔎 Check genes avaibility',
                     help='Sometimes genes do not have the same name in all species or do not exist.'):
            with st.spinner("Checking genes avaibility..."):
                species_list = ['Human', 'Mouse', 'Rat', 'Drosophila', 'Zebrafish']
                results_gene_list = []
                data = []
                for gene_input in gene_list:
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
                df = pd.DataFrame(data, columns=species_columns)
                st.dataframe(df, hide_index=True)

    with colprom2:
        tab1, tab2 = st.tabs(['Default', 'Advance'])

        with tab1:
            with st.form("Default"):

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

                # Run Promoter Finder
                if st.form_submit_button(f"🧬 :blue[**Step 1.5**] Extract {prom_term}", help='(~5sec/gene)'):
                    with colprom1:
                        with st.spinner("Finding promoters..."):
                            gene_ids = gene_id_entry.strip().split("\n")
                            upstream = int(upstream_entry)
                            st.session_state['upstream'] = upstream
                            downstream = int(downstream_entry)
                            try:
                                result_promoter = find_promoters(gene_ids, species, upstream, downstream)
                                st.success(f"{prom_term} extraction complete !")
                                st.toast(f"{prom_term} extraction complete !", icon='😊')
                            except Exception as e:
                                st.error(f"Error finding {prom_term}: {str(e)}")

        with tab2:
            with st.form("Advance"):
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

                if st.form_submit_button("🧬 :blue[**Step 1.4**] Extract sequences", help="(~5sec/seq)"):
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
                                                except Exception as e:
                                                    st.error(f"Error finding {gene_ids}: {str(e)}")

    # Promoter output state
    st.divider()
    st.subheader(':blue[Step 2] Binding Sites Finder')
    promcol1, promcol2 = st.columns([0.9, 0.1], gap='small')
    with promcol1:
        if 'result_promoter' not in st.session_state:
            st.markdown("🔹 :blue[**Step 2.1**] Sequences:")
            result_promoter = st.text_area("🔹 :blue[**Step 2.1**] Sequences:",
                                           value="If Step 1 not used, paste sequences here (FASTA required for multiple sequences).",
                                           label_visibility='collapsed')
        else:
            st.markdown("🔹 :blue[**Step 2.1**] Sequences:", help='Copy: Click in sequence, CTRL+A, CTRL+C')
            result_promoter_text = "\n".join(st.session_state['result_promoter'])
            result_promoter = st.text_area("🔹 :blue[**Step 2.1**] Sequences:", value=result_promoter_text,
                                           label_visibility='collapsed')
    with promcol2:
        st.markdown('')
        st.markdown('')
        st.markdown('')
        current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        txt_output = f"{result_promoter}"
        st.download_button(label="💾 Download (.fasta)", data=txt_output,
                           file_name=f"Sequences_{current_date_time}.fasta", mime="text/plain")
