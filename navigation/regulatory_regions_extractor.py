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

import pandas as pd
import requests
import streamlit as st

from tfinder import NCBIdna


def prom_extractor_page():
    result = extract()
    dna_sequence = fasta(result)
    return dna_sequence


def extract():
    st.subheader(':blue[Step 1] Promoter and Terminator Extractor')
    colprom1, colprom2 = st.columns([0.8, 1.2], gap="small")

    # Extraction of DNA sequence
    with colprom1:
        st.info("💡 If you have a FASTA sequence, go to :blue[**Step 2**]")

        result_promoter = []
        upstream_entry = []

        # Gene ID
        st.markdown("🔹 :blue[**Step 1.1**] Gene ID:", help='NCBI gene name and NCBI gene ID allowed')
        gene_id_entry = st.text_area("🔹 :blue[**Step 1.1**] Gene ID:", value="PRKN\n351\nNM_003130.4",
                                     label_visibility='collapsed')
        gene_ids = gene_id_entry.strip().split("\n")

    with colprom2:
        tab1, tab2 = st.tabs(['Default', 'Advance'])

        with tab1:
            # Species
            st.markdown("🔹 :blue[**Step 1.2**] Species of gene names and sliced variants:")
            col1, col2, col3 = st.columns(3)
            species = col2.selectbox("Species:", ["Human", "Mouse", "Rat", "Drosophila", "Zebrafish"], index=0)
            genomes = {
                "Human": {"Current": "GrCh38", "Previous": "GrCh37"},
                "Mouse": {"Current": "GRCm39", "Previous": "GRCm38"},
                "Rat": {"Current": "GRCr8", "Previous": "mRatBN7.2"},
                "Drosophila": {"Current": "BDGP6", "Previous": "BDGP5"},
                "Zebrafish": {"Current": "GRCz11", "Previous": "GRCz10"}
            }
            genome_choices = genomes.get(species, {})
            gr = col1.selectbox(
                "Genome Version:",
                options=list(genome_choices.keys()),
                format_func=lambda x: f"{x} ({genome_choices[x]})",
                index=0,
                help=f'Selection of genome versions for {species}'
            )

            col3.markdown("")
            col3.markdown("")
            all_slice_forms = col3.toggle(label='All variants')

            st.markdown("🔹 :blue[**Step 1.3**] Sequence type:")
            seq_type_display = st.radio(
                "🔹 :blue[**Step 1.3**] Sequence type:",
                ('Promoter', 'Terminator', "RNA (exon + intron)", "mRNA (exon)"),
                horizontal=True,
                label_visibility='collapsed'
            )
            seq_type_mapping = {
                "Promoter": "Promoter",
                "Terminator": "Terminator",
                "RNA (exon + intron)": "RNA",
                "mRNA (exon)": "mRNA"
            }

            seq_type = seq_type_mapping[seq_type_display]

            if seq_type in ['Promoter', 'Terminator']:
                if seq_type == 'Promoter':
                    st.markdown("🔹 :blue[**Step 1.4**] Upstream/downstream from the TSS (bp)")
                elif seq_type == "Terminator":
                    st.markdown("🔹 :blue[**Step 1.4**] Upstream/downstream from gene end (bp)")

                updown_slide = st.slider("🔹 :blue[**Step 1.4**] Upstream/downstream", -10000, 10000,
                                         (-2000, 2000), step=100, label_visibility='collapsed', disabled=True if seq_type in ['mRNA', "RNA"] else False)
                if seq_type == 'Promoter':
                    st.write("Upstream: ", min(updown_slide), " bp from TSS | Downstream: ", max(updown_slide),
                             " bp from TSS")
                elif seq_type == "Terminator":
                    st.write("Upstream: ", min(updown_slide), " bp from gene end | Downstream: ", max(updown_slide),
                             " bp from gene end")
            else:
                updown_slide = [0, 0]

            upstream_entry = -min(updown_slide)
            downstream_entry = max(updown_slide)

            upstream = int(upstream_entry)
            st.session_state['upstream'] = upstream
            downstream = int(downstream_entry)

            # Run Promoter Finder
            if st.button(f"🧬 :blue[**Step 1.5**] Extract {seq_type}", help='(~5sec/gene)'):
                response = requests.get(
                    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=nos2[Gene%20Name]+AND+human[Organism]&retmode=json&rettype=xml')

                ncbi_status = True if response.status_code == 200 else False

                if ncbi_status is True:
                    st.session_state['result_promoter_text'] = ""
                    with st.spinner('Please wait...'):
                        with colprom1:

                            pbar = st.progress(0,
                                               text='**:blue[Extract sequence...] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                            for i, gene_id in enumerate(gene_ids):
                                pbar.progress(i / len(gene_ids),
                                              text=f'**:blue[Extract sequence... {gene_id}] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                                all_variants, message = NCBIdna(gene_id, species, seq_type.lower(), upstream, downstream, gr.lower(),
                                                                          all_slice_forms=True if all_slice_forms else False).find_sequences()
                                if message == "OK":
                                    for nm_id, data in all_variants.items():
                                        exon_coords = data.get('exon_coords', [])
                                        upstream_seq = data.get('upstream')

                                        dna_sequence = (
                                            f">{nm_id} {data.get('gene_name', 'Unknown')} | "
                                            f"{data.get('genomic_info', 'No Title')} {data.get('chraccver', 'No chraccver')} | "
                                            f"Strand: {data.get('strand', 'Unknown')} | "
                                            f"Type: {data.get('seq_type', 'Unknown')} | "
                                        )

                                        if data.get('seq_type') == 'promoter':
                                            coord = (
                                                f"TSS (on chromosome): {exon_coords[0][0] + 1} "
                                                f"{f'| TSS (on sequence): {upstream_seq}' if upstream_seq is not None else ''}\n"
                                            )
                                        elif data.get('seq_type') == 'terminator':
                                            coord = (
                                                f"Gene end (on chromosome): {exon_coords[-1][1]} "
                                                f"{f'| Gene end (on sequence): {upstream_seq}' if upstream_seq is not None else ''}\n"
                                            )
                                        else:
                                            coord = (
                                                f"TSS (on chromosome): {exon_coords[0][0] + 1} "
                                                f"| Gene end (on chromosome): {exon_coords[-1][1]}\n"
                                            )

                                        sequence = f"{data.get('sequence', 'No sequence available')}\n"

                                        result_promoter.append(dna_sequence + coord + sequence)

                                        pbar.progress(
                                            (i + 1) / len(gene_ids),
                                            text=(
                                                f"**:blue[Extract sequence... {gene_id}] "
                                                f"⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**"
                                            )
                                        )

                                        st.toast(
                                            f"{data.get('seq_type', 'Unknown')} **{gene_id}** from **{species}** extracted",
                                            icon='🧬'
                                        )

                                else:
                                    st.error(message)
                                    continue

                            result_promoter_text = "\n".join(result_promoter)

                            st.session_state['result_promoter_text'] = result_promoter_text

                            st.success(f"{seq_type} extraction complete !")
                            st.toast(f"{seq_type} extraction complete !", icon='😊')

                elif ncbi_status is False:
                    st.warning("⚠ NCBI servers are under maintenance or have an error")

        with tab2:
            # Advance mode extraction
            data_df = pd.DataFrame(
                {
                    "Gene": gene_ids,
                    "human": [False] * len(gene_ids),
                    "mouse": [False] * len(gene_ids),
                    "rat": [False] * len(gene_ids),
                    "drosophila": [False] * len(gene_ids),
                    "zebrafish": [False] * len(gene_ids),
                    "promoter": [False] * len(gene_ids),
                    "terminator": [False] * len(gene_ids),
                    "rna": [False] * len(gene_ids),
                    "mrna": [False] * len(gene_ids),
                    "current": [False] * len(gene_ids),
                    "previous": [False] * len(gene_ids),
                }
            )

            species_list = ['human', 'mouse', 'rat', 'drosophila', 'zebrafish']
            search_types = ['promoter', 'terminator', "rna", "mrna"]
            genome_type = ["current", "previous"]

            st.markdown('**🔹 :blue[Step 1.2]** Select species for all genes:',
                        help='Checking a box allows you to check all the corresponding boxes for each gene. Warning: if you have manually checked boxes in the table, they will be reset.')

            species1, species2, species3, species4, species5 = st.columns(5)

            with species1:
                all_human = st.toggle("Human")
            with species2:
                all_mouse = st.toggle("Mouse")
            with species3:
                all_rat = st.toggle("Rat")
            with species4:
                all_droso = st.toggle("Drosophila")
            with species5:
                all_zebra = st.toggle("Zebrafish")

            st.markdown('**🔹 :blue[Step 1.2]** Select regions for all genes:',
                        help='Checking a box allows you to check all the corresponding boxes for each gene. Warning: if you have manually checked boxes in the table, they will be reset.')

            region1, region2, region3, region4 = st.columns(4)

            all_prom = region1.toggle("Promoter")
            all_term = region2.toggle("Terminator")
            all_rna = region3.toggle("RNA")
            all_mrna = region4.toggle("mRNA")

            st.markdown('**🔹 :blue[Step 1.2]** Select genome version for all genes:',
                        help="Human: Current  ->  GrCh38, Previous  ->  GrCh37\n\n"
                             "Mouse: Current  ->  GRCm39, Previous  ->  GRCm38\n\n"
                             "Rat: Current  ->  GRCr8, Previous  ->  mRatBN7.2\n\n"
                             "Drosophila: Current  ->  BDGP6, Previous  ->  BDGP5\n\n"
                             "Zebrafish: Current  ->  GRCz11, Previous  ->  GRCz10\n\n")

            gr1, gr2 = st.columns(2)

            with gr1:
                gr_current = st.toggle("Current")
            with gr2:
                gr_previous = st.toggle("Previous")

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
            if all_rna:
                data_df["rna"] = True
            if all_mrna:
                data_df["mrna"] = True
            if gr_current:
                data_df['current'] = True
            if gr_previous:
                data_df['previous'] = True

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
                    ),
                    "rna": st.column_config.CheckboxColumn(
                        "RNA",
                        default=False,
                    ),
                    "mrna": st.column_config.CheckboxColumn(
                        "mRNA",
                        default=False,
                    ),
                    "current": st.column_config.CheckboxColumn(
                        "Current genome",
                        default=False,
                    ),
                    "previous": st.column_config.CheckboxColumn(
                        "Previous genome",
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

            if st.button("🧬 :blue[**Step 1.4**] Extract sequences", help="(~5sec/seq)", key='Advance'):
                response = requests.get(
                    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=nos2[Gene%20Name]+AND+human[Organism]&retmode=json&rettype=xml')

                ncbi_status = True if response.status_code == 200 else False

                if ncbi_status is True:
                    st.session_state['result_promoter_text'] = ""
                    with colprom1:
                        st.session_state['upstream'] = upstream_entry
                        upstream = int(upstream_entry)
                        downstream = int(downstream_entry)
                        pbar = st.progress(0,
                                           text='**:blue[Extract sequence...] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                        for i, gene_info in enumerate(data_dff.itertuples(index=False)):
                            gene_id = gene_info.Gene
                            if gene_id.isdigit() or gene_id.startswith('XM_') or gene_id.startswith(
                                    'NM_') or gene_id.startswith('XR_') or gene_id.startswith('NR_'):
                                for genome in genome_type:
                                    for search_type in search_types:
                                        if getattr(gene_info, f'{search_type}') and getattr(gene_info, f'{genome}'):
                                            seq_type = search_type.capitalize()

                                            pbar.progress((i + 1) / len(data_dff),
                                                          text=f'**:blue[Extract sequence... {seq_type} **{gene_id}** from **{species}**] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')

                                            all_variants, message = NCBIdna(gene_id, seq_type=seq_type.lower(),
                                                                            upstream=upstream, downstream=downstream, genome_version=gr.lower(),
                                                                            all_slice_forms=True if all_slice_forms else False).find_sequences()
                                            if message == "OK":
                                                for nm_id, data in all_variants.items():
                                                    exon_coords = data.get('exon_coords', [])
                                                    upstream_seq = data.get('upstream')

                                                    if not exon_coords:
                                                        print(
                                                            f"Erreur: 'exon_coords' est manquant ou vide pour {nm_id}.")
                                                        continue

                                                    dna_sequence = (
                                                        f">{nm_id} {data.get('gene_name', 'Unknown')} | "
                                                        f"{data.get('genomic_info', 'No Title')} {data.get('chraccver', 'No chraccver')} | "
                                                        f"Strand: {data.get('strand', 'Unknown')} | "
                                                        f"Type: {data.get('seq_type', 'Unknown')} | "
                                                    )

                                                    if data.get('seq_type') == 'promoter':
                                                        coord = (
                                                            f"TSS (on chromosome): {exon_coords[0][0] + 1} "
                                                            f"{f'| TSS (on sequence): {upstream_seq}' if upstream_seq is not None else ''}\n"
                                                        )
                                                    elif data.get('seq_type') == 'terminator':
                                                        coord = (
                                                            f"Gene end (on chromosome): {exon_coords[-1][1]} "
                                                            f"{f'| Gene end (on sequence): {upstream_seq}' if upstream_seq is not None else ''}\n"
                                                        )
                                                    else:
                                                        coord = (
                                                            f"TSS (on chromosome): {exon_coords[0][0] + 1} "
                                                            f"| Gene end (on chromosome): {exon_coords[-1][1]}\n"
                                                        )

                                                    sequence = f"{data.get('sequence', 'No sequence available')}\n"

                                                    result_promoter.append(dna_sequence + coord + sequence)

                                                    st.toast(
                                                        f"{data.get('seq_type', 'Unknown')} **{gene_id}** from **{species}** extracted",
                                                        icon='🧬'
                                                    )
                                            else:
                                                st.error(message)
                                                continue

                            else:
                                for genome in genome_type:
                                    for species in species_list:
                                        for search_type in search_types:
                                            if getattr(gene_info, f'{species}') and getattr(gene_info,
                                                                                            f'{search_type}') and getattr(
                                                gene_info, f'{genome}'):
                                                seq_type = search_type.capitalize()

                                                pbar.progress((i + 1) / len(data_dff),
                                                              text=f'**:blue[Extract sequence... {seq_type} **{gene_id}** from **{species.capitalize()}**] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')

                                                all_variants, message = NCBIdna(gene_id, species, seq_type=seq_type.lower(),
                                                                                upstream=upstream,
                                                                                downstream=downstream,
                                                                                genome_version=gr.lower(),
                                                                                all_slice_forms=True if all_slice_forms else False).find_sequences()
                                                if message == "OK":
                                                    for nm_id, data in all_variants.items():
                                                        exon_coords = data.get('exon_coords', [])
                                                        upstream_seq = data.get('upstream')

                                                        if not exon_coords:
                                                            print(
                                                                f"Erreur: 'exon_coords' est manquant ou vide pour {nm_id}.")
                                                            continue

                                                        dna_sequence = (
                                                            f">{nm_id} {data.get('gene_name', 'Unknown')} | "
                                                            f"{data.get('genomic_info', 'No Title')} {data.get('chraccver', 'No chraccver')} | "
                                                            f"Strand: {data.get('strand', 'Unknown')} | "
                                                            f"Type: {data.get('seq_type', 'Unknown')} | "
                                                        )

                                                        if data.get('seq_type') == 'promoter':
                                                            coord = (
                                                                f"TSS (on chromosome): {exon_coords[0][0] + 1} "
                                                                f"{f'| TSS (on sequence): {upstream_seq}' if upstream_seq is not None else ''}\n"
                                                            )
                                                        elif data.get('seq_type') == 'terminator':
                                                            coord = (
                                                                f"Gene end (on chromosome): {exon_coords[-1][1]} "
                                                                f"{f'| Gene end (on sequence): {upstream_seq}' if upstream_seq is not None else ''}\n"
                                                            )
                                                        else:
                                                            coord = (
                                                                f"TSS (on chromosome): {exon_coords[0][0] + 1} "
                                                                f"| Gene end (on chromosome): {exon_coords[-1][1]}\n"
                                                            )

                                                        sequence = f"{data.get('sequence', 'No sequence available')}\n"

                                                        result_promoter.append(dna_sequence + coord + sequence)

                                                        st.toast(
                                                            f"{data.get('seq_type', 'Unknown')} **{gene_id}** from **{species}** extracted",
                                                            icon='🧬'
                                                        )
                                                else:
                                                    st.error(message)
                                                    continue

                        result_promoter_text = "\n".join(result_promoter)
                        st.session_state['result_promoter_text'] = result_promoter_text
                        st.success(f"{seq_type} extraction complete !")
                        st.toast(f"{seq_type} extraction complete !", icon='😊')

                elif ncbi_status is False:
                    st.warning("⚠ NCBI servers are under maintenance or have an error")


def fasta(result_promoter_text=None):
    st.divider()
    promcol1, promcol2 = st.columns([0.75, 0.25], gap='small')

    with promcol2:
        st.markdown('')
        uploaded_files = st.file_uploader(
            "Upload FASTA/TXT files:",
            type=["fasta", "txt"],
            accept_multiple_files=True,
            help="Upload multiple FASTA/TXT files"
        )
        if uploaded_files:
            concatenated_sequences = ""
            for uploaded_file in uploaded_files:
                uploaded_content = uploaded_file.getvalue().decode("utf-8")
                concatenated_sequences += uploaded_content + "\n"
            st.session_state['result_promoter_text'] = concatenated_sequences.strip()

    with promcol1:
        st.markdown("🔹 :blue[**Step 2.1**] Sequences:", help='Copy: Click in sequence, CTRL+A, CTRL+C')

        dna_sequence = st.text_area(
            "🔹 :blue[**Step 2.1**] Sequences:",
            value="" if 'result_promoter_text' not in st.session_state
            else st.session_state['result_promoter_text'],
            label_visibility='collapsed',
            height=125
        )
        st.session_state['result_promoter_text'] = dna_sequence

        if '>' in dna_sequence:
            num_sequences = dna_sequence.count('>')
            st.markdown(f"**Number of sequences**: {num_sequences}")
        elif dna_sequence.strip():
            num_sequences = 1
            st.markdown("**Single sequence detected**")
        else:
            num_sequences = 0
            st.markdown("**No sequence entered**")

    current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    txt_output = f"{dna_sequence}"
    promcol2.download_button(
        label="💾 Download sequence (.fasta)",
        data=txt_output,
        file_name=f"Sequences_{current_date_time}.fasta",
        mime="text/plain"
    )

    return dna_sequence
