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
import streamlit as st

from tfinder import NCBIdna


def prom_extractor_page():
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

        # Verify if gene is available for all species
        if st.button('🔎 Check genes avaibility',
                     help='Sometimes genes do not have the same name in all species or do not exist.'):
            with st.spinner('Please wait...'):
                species_list = ['ID', 'Human', 'Mouse', 'Rat', 'Drosophila', 'Zebrafish']
                gene_disponibility_output = []
                pbar = st.progress(0,
                                   text='**:blue[Analyse genes...] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                for i, gene_id in enumerate(gene_ids):
                    pbar.progress(i / len(gene_ids),
                                  text=f'**:blue[Analyse genes... {gene_id}] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                    gene_disponibility_output.append(NCBIdna.analyse_gene(gene_id))
                    pbar.progress((i + 1) / len(gene_ids),
                                  text=f'**:blue[Analyse genes... {gene_id}] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')

                species_columns = ['Gene'] + species_list
                gene_disponibility_output = pd.DataFrame(gene_disponibility_output, columns=species_columns)

                st.session_state['gene_disponibility_output'] = gene_disponibility_output

        if 'gene_disponibility_output' in st.session_state:
            st.dataframe(st.session_state['gene_disponibility_output'], hide_index=True)

    with colprom2:
        tab1, tab2 = st.tabs(['Default', 'Advance'])

        with tab1:
            # Species
            st.markdown("🔹 :blue[**Step 1.2**] Species of gene names and sliced variants:")
            col1, col2, = st.columns(2)
            with col1:
                species = st.selectbox("🔹 :blue[**Step 1.2**] Select species of gene names:",
                                       ["Human", "Mouse", "Rat", "Drosophila", "Zebrafish"], index=0,
                                       label_visibility='collapsed')

            with col2:
                all_variants = st.toggle('All variant')

            # Upstream/Downstream Promoter
            st.markdown("🔹 :blue[**Step 1.3**] Regulatory region:")
            prom_term = st.radio("🔹 :blue[**Step 1.3**] Regulatory region:", ('Promoter', 'Terminator'),
                                 horizontal=True,
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

            if 'button_clicked' not in st.session_state:
                st.session_state.button_clicked = False

            # Run Promoter Finder
            if st.button(f"🧬 :blue[**Step 1.5**] Extract {prom_term}", help='(~5sec/gene)'):
                with st.spinner('Please wait...'):
                    with colprom1:

                        st.session_state.button_clicked = True

                        pbar = st.progress(0,
                                           text='**:blue[Extract sequence...] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                        for i, gene_id in enumerate(gene_ids):
                            pbar.progress(i / len(gene_ids),
                                          text=f'**:blue[Extract sequence... {gene_id}] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                            result_promoter_output = NCBIdna(gene_id, prom_term, upstream, downstream,
                                                             species,
                                                             all_slice_forms=True if all_variants else False).find_sequences()
                            if not str(result_promoter_output).startswith('P'):
                                pbar.progress((i + 1) / len(gene_ids),
                                              text=f'**:blue[Extract sequence... {gene_id}] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                                st.toast(f"{prom_term} **{gene_id}** from **{species}** extracted", icon='🧬')

                                result_promoter.append(result_promoter_output)
                            else:
                                st.error(result_promoter_output)
                                continue

                        result_promoter_text = "\n".join(result_promoter)

                        st.session_state['result_promoter_text'] = result_promoter_text

                        st.success(f"{prom_term} extraction complete !")
                        st.toast(f"{prom_term} extraction complete !", icon='😊')

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
                }
            )

            species_list = ['human', 'mouse', 'rat', 'drosophila', 'zebrafish']
            search_types = ['promoter', 'terminator']

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

            region1, region2 = st.columns(2)

            with region1:
                all_prom = st.toggle("Promoter")
            with region2:
                all_term = st.toggle("Terminator")

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

            if st.button("🧬 :blue[**Step 1.4**] Extract sequences", help="(~5sec/seq)", key='Advance'):
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
                            for search_type in search_types:
                                if getattr(gene_info, f'{search_type}'):
                                    prom_term = search_type.capitalize()

                                    pbar.progress((i + 1) / len(data_dff),
                                                  text=f'**:blue[Extract sequence... {prom_term} **{gene_id}** from **{species}**] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')

                                    result_promoter_output = NCBIdna(gene_id, prom_term, upstream,
                                                                     downstream).find_sequences()

                                    if not result_promoter_output.startswith('P'):
                                        st.toast(f'{prom_term} **{gene_id}** from **{species}** extracted',
                                                 icon='🧬')
                                        result_promoter.append(result_promoter_output)
                                        pass

                                    else:
                                        st.error(result_promoter_output)
                                        continue

                        else:
                            for species in species_list:
                                for search_type in search_types:
                                    if getattr(gene_info, f'{species}') and getattr(gene_info,
                                                                                    f'{search_type}'):
                                        prom_term = search_type.capitalize()

                                        pbar.progress((i + 1) / len(data_dff),
                                                      text=f'**:blue[Extract sequence... {prom_term} **{gene_id}** from **{species.capitalize()}**] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')

                                        result_promoter_output = NCBIdna(gene_id, prom_term, upstream,
                                                                         downstream,
                                                                         species).find_sequences()

                                        if not result_promoter_output.startswith('P'):
                                            st.toast(
                                                f'{prom_term} **{gene_id}** from **{species.capitalize()}** extracted',
                                                icon='🧬')
                                            result_promoter.append(result_promoter_output)
                                            pass

                                        else:
                                            st.error(result_promoter_output)
                                            continue

                    result_promoter_text = "\n".join(result_promoter)
                    st.session_state['result_promoter_text'] = result_promoter_text
                    st.success(f"{prom_term} extraction complete !")
                    st.toast(f"{prom_term} extraction complete !", icon='😊')

    # Promoter output state
    st.divider()
    promcol1, promcol2 = st.columns([0.9, 0.1], gap='small')
    with promcol1:
        st.markdown("🔹 :blue[**Step 2.1**] Sequences:", help='Copy: Click in sequence, CTRL+A, CTRL+C')
        if 'result_promoter_text' not in st.session_state:
            result_promoter_text = ''
            st.session_state['result_promoter_text'] = result_promoter_text
        dna_sequence = st.text_area("🔹 :blue[**Step 2.1**] Sequences:",
                                    value=st.session_state['result_promoter_text'],
                                    placeholder='If Step 1 not used, paste sequences here (FASTA required for multiple sequences).',
                                    label_visibility='collapsed', height=125)

    with promcol2:
        st.markdown('')
        st.markdown('')
        st.markdown('')
        current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        txt_output = f"{dna_sequence}"
        st.download_button(label="💾 Download (.fasta)", data=txt_output,
                           file_name=f"Sequences_{current_date_time}.fasta", mime="text/plain")
