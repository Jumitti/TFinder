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
import re
import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

import altair as alt
import pandas as pd
import requests
import streamlit as st
from stqdm import stqdm

from tfinder import IMO
from tfinder import NCBIdna


def email(excel_file, csv_file, txt_output, email_receiver, body, jaspar):
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

        attachment_csv = MIMEBase('application', 'octet-stream')
        attachment_csv.set_payload(csv_file)
        encoders.encode_base64(attachment_csv)
        attachment_csv.add_header('Content-Disposition', 'attachment',
                                  filename=f'Results_TFinder_{current_date_time}.csv')
        msg.attach(attachment_csv)

        if jaspar == 'PWM' or jaspar == 'Individual Motif':
            image = MIMEImage(st.session_state['weblogo'].read(), name=f'LOGOMAKER_{current_date_time}.jpg')
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


def result_table_output(source):
    score_range = source['Rel Score'].astype(float)
    ystart = score_range.min() - 0.02
    ystop = score_range.max() + 0.02

    source['Gene_Region'] = source['Gene'] + " " + source['Species'] + " " + source['Region']
    source['Beginning of sequences'] = source['Position']

    scale = alt.Scale(scheme='category10')
    color_scale = alt.Color("Gene_Region:N", scale=scale)

    gene_region_selection = alt.selection_point(fields=['Gene_Region'], on='click', bind='legend')

    y_dropdown = alt.binding_select(
        options=['Rel Score', 'Rel Score Adj'],
        name='(Y-axis) Relative Score: '
    )
    ycol_param = alt.param(value='Rel Score', bind=y_dropdown, name="y_axis")

    base_chart = alt.Chart(source).mark_circle().encode(
        x=alt.X('Beginning of sequences:Q', title='Position (bp)', axis=alt.Axis(labelAngle=0)),
        y=alt.Y('y:Q', axis=alt.Axis(title='Relative Score'), scale=alt.Scale(domain=[ystart, ystop])),
        color=alt.condition(gene_region_selection, color_scale, alt.value('lightgray')),
        tooltip=['Sequence', 'Position'] + (['Rel Position'] if "Rel Position" in source else []) + (
            ['Ch Position'] if "Ch Position" in source else []) + ['Rel Score'] + ['Score'] +
                    ['Rel Score Adj'] + ['Score Adj'] + (
                    ['p-value'] if 'p-value' in source else []) + ['Gene', 'Species', 'Region'],
        opacity=alt.condition(gene_region_selection, alt.value(0.8), alt.value(0.2))
    ).transform_calculate(
        y=f'datum[{ycol_param.name}]'
    ).properties(
        width=600, height=400
    )

    if "Rel Position" in source:
        secondary_axis = alt.Chart(source).mark_rule(opacity=0).encode(
            x=alt.X('Rel Position:Q', title='Position from TSS/Gene end (bp)',
                    axis=alt.Axis(orient='top')),
        ).interactive()

        chart = alt.layer(base_chart, secondary_axis).resolve_scale(x='independent', y='shared')
    else:
        chart = base_chart

    chart = chart.add_params(
        gene_region_selection, ycol_param
    ).interactive()

    st.altair_chart(chart, theme=None, use_container_width=True)







def aio_page():
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
            all_variants = col3.toggle(label='All variants')

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

            # Run Promoter Finder
            if st.button(f"🧬 :blue[**Step 1.5**] Extract {prom_term}", help='(~5sec/gene)'):
                response = requests.get(
                    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=nos2[Gene%20Name]+AND+human[Organism]&retmode=json&rettype=xml')

                ncbi_status = True if response.status_code == 200 else False

                if ncbi_status is True:
                    with st.spinner('Please wait...'):
                        with colprom1:

                            pbar = st.progress(0,
                                               text='**:blue[Extract sequence...] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                            for i, gene_id in enumerate(gene_ids):
                                pbar.progress(i / len(gene_ids),
                                              text=f'**:blue[Extract sequence... {gene_id}] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                                result_promoter_output, message = NCBIdna(gene_id, prom_term, upstream, downstream,
                                                                          species, gr.lower(),
                                                                          all_slice_forms=True if all_variants else False).find_sequences()
                                if message == "OK":
                                    pbar.progress((i + 1) / len(gene_ids),

                                                  text=f'**:blue[Extract sequence... {gene_id}] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')
                                    st.toast(f"{prom_term} **{gene_id}** from **{species}** extracted", icon='🧬')

                                    result_promoter.append(result_promoter_output)
                                else:
                                    st.error(message)
                                    continue

                            result_promoter_text = "\n".join(result_promoter)

                            st.session_state['result_promoter_text'] = result_promoter_text

                            st.success(f"{prom_term} extraction complete !")
                            st.toast(f"{prom_term} extraction complete !", icon='😊')

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
                    "current": [False] * len(gene_ids),
                    "previous": [False] * len(gene_ids),
                }
            )

            species_list = ['human', 'mouse', 'rat', 'drosophila', 'zebrafish']
            search_types = ['promoter', 'terminator']
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

            region1, region2 = st.columns(2)

            with region1:
                all_prom = st.toggle("Promoter")
            with region2:
                all_term = st.toggle("Terminator")

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
                                            prom_term = search_type.capitalize()

                                            pbar.progress((i + 1) / len(data_dff),
                                                          text=f'**:blue[Extract sequence... {prom_term} **{gene_id}** from **{species}**] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')

                                            result_promoter_output, message = NCBIdna(gene_id, prom_term, upstream,
                                                                                      downstream,
                                                                                      genome_version=genome).find_sequences()

                                            if message == "OK":
                                                st.toast(f"{prom_term} **{gene_id}** from **{species}** extracted",
                                                         icon='🧬')

                                                result_promoter.append(result_promoter_output)
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
                                                prom_term = search_type.capitalize()

                                                pbar.progress((i + 1) / len(data_dff),
                                                              text=f'**:blue[Extract sequence... {prom_term} **{gene_id}** from **{species.capitalize()}**] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**')

                                                result_promoter_output, message = NCBIdna(gene_id, prom_term, upstream,
                                                                                          downstream, species,
                                                                                          genome).find_sequences()

                                                if message == "OK":
                                                    st.toast(f"{prom_term} **{gene_id}** from **{species}** extracted",
                                                             icon='🧬')

                                                    result_promoter.append(result_promoter_output)
                                                else:
                                                    st.error(message)
                                                    continue

                        result_promoter_text = "\n".join(result_promoter)
                        st.session_state['result_promoter_text'] = result_promoter_text
                        st.success(f"{prom_term} extraction complete !")
                        st.toast(f"{prom_term} extraction complete !", icon='😊')

                elif ncbi_status is False:
                    st.warning("⚠ NCBI servers are under maintenance or have an error")

    # Promoter output state
    st.divider()
    st.subheader(':blue[Step 2] Binding Sites Finder')
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

    # Promoter detection information
    lines = dna_sequence
    dna_sequences = []
    if lines.startswith(("A", "T", "C", "G", "N", "a", "t", "c", "g", "n")):
        dna_sequence = lines.upper()
        isfasta = IMO.is_dna(dna_sequence)
        name = "n.d."
        species = "n.d"
        region = "n.d"
        strand = "n.d"
        tss_ch = 0
        dna_sequences.append((name, dna_sequence, species, region, strand, tss_ch))
    elif lines.startswith(">"):
        lines = dna_sequence.split("\n")
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.startswith(">"):
                species_prom = ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'Drosophila melanogaster',
                                'Danio rerio']
                promoter_name = line[1:]
                words = promoter_name.lstrip('>').split()
                pattern = r">(\w+)\s+(\w+)\s+\|"
                match = re.search(pattern, line)
                if match:
                    name = words[0] + ' ' + words[1]
                else:
                    name = words[0]
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
                dna_sequence = lines[i + 1].upper()
                isfasta = IMO.is_dna(dna_sequence)

                match = re.search(r"Strand:\s*(\w+)", line)
                if match:
                    strand = match.group(1).lower()
                    if strand not in ["plus", "minus"]:
                        strand = "n.d"
                else:
                    strand = "n.d"

                match = re.search(r"TSS \(on chromosome\):\s*(\d+)", line)
                if match:
                    tss_ch = match.group(1)
                else:
                    tss_ch = 0

                dna_sequences.append((name, dna_sequence, found_species, region, strand, tss_ch))
                i += 1
            else:
                i += 1
    elif not lines.startswith(("A", "T", "C", "G", "N", "a", "t", "c", "g", "n", "I", "i", "")):
        isfasta = True
    else:
        isfasta = False

    total_sequences_region_length = sum(len(dna_sequence) for _, dna_sequence, _, _, _, _ in dna_sequences)
    total_sequences = len(dna_sequences)

    # RE entry
    REcol1, REcol2 = st.columns([0.30, 0.70])
    with REcol1:
        st.markdown('🔹 :blue[**Step 2.2**] Motif type:')
        jaspar = st.radio('🔹 :blue[**Step 2.2**] Motif type:',
                          ('Individual Motif', 'JASPAR_ID', 'PWM'),
                          label_visibility='collapsed')
    if jaspar == 'JASPAR_ID':
        isUIPAC = True
        with REcol1:
            st.markdown("🔹 :blue[**Step 2.3**] JASPAR ID:")
            jaspar_id = st.text_input("🔹 :blue[**Step 2.3**] JASPAR ID:",
                                      value="MA0106.1" if 'JASPAR_ID_save' not in st.session_state
                                      else st.session_state['JASPAR_ID_save'],
                                      label_visibility='collapsed')
            st.session_state['JASPAR_ID_save'] = jaspar_id
            if jaspar_id:
                TF_name, TF_species, matrix, weblogo = IMO.matrix_extraction(jaspar_id)
                if TF_name != 'not found':
                    st.success(f"{TF_species} transcription factor {TF_name}")
                    with REcol2:
                        st.image(weblogo)
                    button = False
                    error_input_im = True
                else:
                    button = True
                    error_input_im = False
                    st.error('Wrong JASPAR_ID')
            else:
                button = True
                error_input_im = False
                st.warning('Please enter a JASPAR_ID')

    elif jaspar == 'PWM':
        with REcol1:
            st.markdown('🔹 :blue[**Step 2.2bis**] Matrix:')
            matrix_type = st.radio('🔹 :blue[**Step 2.2bis**] Matrix:', ('With FASTA sequences', 'With PWM'),
                                   label_visibility='collapsed')
        if matrix_type == 'With PWM':
            isUIPAC = True
            with REcol2:
                st.markdown("🔹 :blue[**Step 2.3**] Matrix:")
                matrix_str = st.text_area("🔹 :blue[**Step 2.3**] Matrix:",
                                          value="A [ 20.0 0.0 0.0 0.0 0.0 0.0 0.0 100.0 0.0 60.0 20.0 ]\nT [ 60.0 20.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]\nG [ 0.0 20.0 100.0 0.0 0.0 100.0 100.0 0.0 100.0 40.0 0.0 ]\nC [ 20.0 60.0 0.0 100.0 100.0 0.0 0.0 0.0 0.0 0.0 80.0 ]"
                                          if 'MATRIX_STR_save' not in st.session_state else st.session_state[
                                              'MATRIX_STR_save'],
                                          label_visibility='collapsed', height=125)
                with st.expander("How to build your PWM"):
                    st.write("**Row Structure:**\n\n"
                             "The PWM should contain four rows, each representing one of the four nucleotides: A (Adenine), T (Thymine), G (Guanine), and C (Cytosine).\n\n"
                             "**Columns Represent Motif Positions:**\n\n"
                             "Each column corresponds to a position in the DNA sequence or motif. The values in each column indicate the weight or score for each nucleotide at that position.\n\n"
                             "**Decimal Numbers:**\n\n"
                             "Use decimal numbers with a period (.) as the decimal separator. Ensure consistent formatting with at least one decimal place, even if the number is a whole (e.g., 20.0 rather than 20 or 20,0).\n\n"
                             "**Symmetry and Completeness:**\n\n"
                             "Ensure each row has the same number of columns (positions) so that the matrix is complete. In the example provided, each nucleotide row has 11 values, one for each position in the motif.")

                st.session_state['MATRIX_STR_save'] = matrix_str

                lines = matrix_str.split("\n")
                matrix = {}
                if len(lines) > 1:
                    for line in lines:
                        parts = line.split("[")
                        base = parts[0].strip()
                        values = [float(val.strip()) for val in parts[1][:-1].split()]
                        matrix[base] = values

                    try:
                        IMO.has_uniform_column_length(matrix_str)

                        weblogo = IMO.PWM_to_weblogo(matrix_str)
                        st.pyplot(weblogo.fig)
                        logo = io.BytesIO()
                        weblogo.fig.savefig(logo, format='png')
                        logo.seek(0)
                        st.session_state['weblogo'] = logo

                        error_input_im = True
                    except Exception as e:
                        error_input_im = False
                        REcol2.error(e)
                else:
                    error_input_im = False
                    REcol2.warning("Please input your PWM :)")
        else:
            with REcol1:
                st.markdown("🔹 :blue[**Step 2.3**] Sequences:",
                            help='Put FASTA sequences. Same sequence length required ⚠')
                individual_motif = st.text_area("🔹 :blue[**Step 2.3**] Sequences:",
                                                value=">seq1\nCTGCCGGAGGA\n>seq2\nAGGCCGGAGGC\n>seq3\nTCGCCGGAGAC\n>seq4\nCCGCCGGAGCG\n>seq5\nAGGCCGGATCG"
                                                if 'individual_motif_save' not in st.session_state else
                                                st.session_state['individual_motif_save'],
                                                label_visibility='collapsed')
                st.session_state['individual_motif_save'] = individual_motif
                individual_motif = individual_motif.upper()
            isUIPAC = True

            try:
                matrix, weblogo = IMO.individual_motif_pwm(individual_motif)
                matrix_str = ""
                for base, values in matrix.items():
                    values_str = " ".join([f"{val:.4f}" for val in values])
                    matrix_str += f"{base} [ {values_str} ]\n"
                with REcol2:
                    matrix_text = st.text_area('PWM', value=matrix_str, height=125,
                                               help='Copy to use later. Not editable.',
                                               disabled=True)
                    st.pyplot(weblogo.fig)
                    logo = io.BytesIO()
                    weblogo.fig.savefig(logo, format='png')
                    logo.seek(0)
                    st.session_state['weblogo'] = logo
                error_input_im = True
            except Exception as e:
                error_input_im = False
                REcol1.error(e)

    else:
        with REcol1:
            st.markdown("🔹 :blue[**Step 2.3**] Individual motif:", help="IUPAC authorized")
            IUPAC = st.text_input("🔹 :blue[**Step 2.3**] Individual motif (IUPAC authorized):",
                                  value="GGGRNYYYCC" if 'IUPAC_seq' not in st.session_state else
                                  st.session_state[
                                      'IUPAC_seq'],
                                  label_visibility='collapsed')
            st.session_state['IUPAC_seq'] = IUPAC
            IUPAC = IUPAC.upper()

        IUPAC_code = ['A', 'T', 'G', 'C', 'R', 'Y', 'M', 'K', 'W', 'S', 'B', 'D', 'H', 'V', 'N', '-', '.']

        if all(char in IUPAC_code for char in IUPAC):
            isUIPAC = True

            sequences = IMO.generate_iupac_variants(IUPAC, max_variant_allowed=1048576)

            if 'Too many' not in sequences:
                individual_motif = ""
                for i, seq in enumerate(sequences):
                    individual_motif += f">seq{i + 1}\n{seq}\n"

                try:
                    matrix, weblogo = IMO.individual_motif_pwm(individual_motif)
                    matrix_str = ""
                    for base, values in matrix.items():
                        values_str = " ".join([f"{val:.4f}" for val in values])
                        matrix_str += f"{base} [ {values_str} ]\n"
                    with REcol2:
                        matrix_text = st.text_area('PWM', value=matrix_str, height=125,
                                                   help='Copy to use later. Not editable.',
                                                   disabled=True)
                        st.pyplot(weblogo.fig)
                        logo = io.BytesIO()
                        weblogo.fig.savefig(logo, format='png')
                        logo.seek(0)
                        st.session_state['weblogo'] = logo
                    error_input_im = True
                except Exception as e:
                    error_input_im = False
                    REcol1.error(e)
            else:
                st.error(sequences)
                isUIPAC = False

        else:
            isUIPAC = False

    # TSS entry
    BSFcol1, BSFcol2, BSFcol3, BSFcol4 = st.columns([1, 1, 1, 2], gap="small")
    with BSFcol1:
        st.markdown("🔹 :blue[**Step 2.4**] Transcription Start Site (TSS)/gene end at (in bp):",
                    help="Distance of TSS and gene end from begin of sequences. If you use Step 1, it is positive value of upstream")
        if 'upstream' not in st.session_state:
            tss_ge_input = st.number_input(
                "🔹 :blue[**Step 2.4**] Transcription Start Site (TSS)/gene end at (in bp):",
                -10000, 10000, 0, label_visibility="collapsed")
        else:
            tss_ge_input = st.number_input(
                "🔹 :blue[**Step 2.4**] Transcription Start Site (TSS)/gene end at (in bp):",
                -10000, 10000, st.session_state['upstream'], label_visibility="collapsed")

    # Threshold pvalue

    with BSFcol2:
        st.markdown("🔹 :blue[**Step 2.5**] Relative Score threshold")
        auto_thre = st.toggle("Automatic threshold", value=True)
        threshold_entry = st.slider("🔹 :blue[**Step 2.5**] Relative Score threshold", 0.5, 1.0, 0.85,
                                    step=0.05,
                                    label_visibility="collapsed", disabled=auto_thre)
        if auto_thre is True:
            threshold_entry = 0

    with BSFcol3:
        st.markdown("🔹 :blue[**_Experimental_**] Calcul _p-value_", help='Experimental, take more times.')
        pvalue = st.toggle('_p-value_')
        if pvalue:
            if total_sequences > 10:
                st.markdown(
                    '⚠️Proportion of A, T, G, C imposed for the calculation of the p-value for more than 10 sequences. See "Resources" for more information')
                st.markdown('A 0.275 | C 0.225 | G 0.225 | T 0.275')
                calc_pvalue = 'ATGCPreset'
            else:
                pvalue_type = st.radio('Nucleotides proportion:', ['Sequence dependent', 'Imposed'],
                                       horizontal=True)
                if pvalue_type == 'Sequence dependent':
                    st.markdown(
                        '⚠️Proportion of A, T, G, C depending on the proportions in the sequence. See "Resources" for more information')
                    calc_pvalue = 'ATGCProportion'
                else:
                    st.markdown(
                        '⚠️Proportion of A, T, G, C imposed for the calculation of the p-value. See "Resources" for more information')
                    st.markdown('A 0.275 | C 0.225 | G 0.225 | T 0.275')
                    calc_pvalue = 'ATGCPreset'
        else:
            calc_pvalue = None

    with BSFcol4:
        with st.expander("Advanced settings", icon="🔹"):
            pseudocount_type = st.checkbox("Automatic pseudocount (to normalize PWM)", True,
                                           help="TFinder transforms PWM into PSSM log-odds. If a value is strictly equal "
                                                "to 0 in the PWM we cannot perform the transformation into log-odds "
                                                "(log(0) = -inf). To overcome this problem, we introduce a pseudocount "
                                                "calculated automatically as follows: √N ∗ bg[nucleotide] where N "
                                                "represents the total number of sequences used to construct the matrix "
                                                "and bg[nucleotide] represents the background of nucleotide frequencies.\n\n"
                                                "If you want to define it yourself, generally, the pseudocount is placed "
                                                "around 1. In our hands, a pseudocount of 0.2 for each nucleotide is "
                                                "sufficient (and gives the same results as "
                                                "[TFBSTools](https://bioconductor.org/packages/release/bioc/html/TFBSTools.html))."
                                                f"")
            if isUIPAC is True:
                pseudocount_auto, _, _ = IMO.transform_PWM(matrix)
                pc_col1, pc_col2, pc_col3, pc_col4 = st.columns(4, gap="small")
                pseudocount_A = pc_col1.number_input('A',
                                                     value=0.20 if pseudocount_type is False else pseudocount_auto['A'],
                                                     step=0.01, disabled=pseudocount_type)
                pseudocount_T = pc_col2.number_input('T',
                                                     value=0.20 if pseudocount_type is False else pseudocount_auto['T'],
                                                     step=0.01, disabled=pseudocount_type)
                pseudocount_G = pc_col3.number_input('G',
                                                     value=0.20 if pseudocount_type is False else pseudocount_auto['G'],
                                                     step=0.01, disabled=pseudocount_type)
                pseudocount_C = pc_col4.number_input('C',
                                                     value=0.20 if pseudocount_type is False else pseudocount_auto['C'],
                                                     step=0.01, disabled=pseudocount_type)
                pseudocount = {'A': pseudocount_A, 'T': pseudocount_T, 'G': pseudocount_G, 'C': pseudocount_C}

            st.divider()

            bgnf_type = st.checkbox("Automatic background nucleotide frequencies (for Score)", True,
                                    help="TFinder has 2 types of Score (and their respective normalization). For more "
                                         'information on the formula, refer to the article or “Resource”.\n\n'
                                         "The **first type (Score and Rel. Score)**: In order to compare the TFBS found and"
                                         " the TFBS requested independently of the heterogeneity of the nucleotide "
                                         "frequencies, we set for this score only the background frequencies at 0.25 per "
                                         "nucleotide.\n\n"
                                         "The 2nd type (Adj. Score and Rel. Adj. Score): This calculation is based on "
                                         "the same mathematical formula but the nucleotide frequencies of the background "
                                         "are calculated directly for each sequence analyzed. This allows you to obtain "
                                         "an Adj Score. (and Rel. Score Adj.) compensating for biases linked to the "
                                         "heterogeneity of nucleotide frequencies for a given sequence. For the "
                                         "calculation of the Adj Score. (and Rel. Score Adj.), you can let the software "
                                         "configure the background itself. However, it is important to note that "
                                         "(as for the p-value) small sequences (100-500bp) have higher heterogeneity "
                                         "(see Resources). In this case, you can configure the nucleotide frequencies "
                                         "of the background yourself.")
            bg_col1, bg_col2, bg_col3, bg_col4 = st.columns(4, gap="small")
            bgnf_A = bg_col1.number_input('A', value=0.25, min_value=0.00, max_value=1.00, step=0.01,
                                          disabled=bgnf_type)
            bgnf_T = bg_col2.number_input('T', value=0.25, min_value=0.00, max_value=1.00, step=0.01,
                                          disabled=bgnf_type)
            bgnf_G = bg_col3.number_input('G', value=0.25, min_value=0.00, max_value=1.00, step=0.01,
                                          disabled=bgnf_type)
            bgnf_C = bg_col4.number_input('C', value=0.25, min_value=0.00, max_value=1.00, step=0.01,
                                          disabled=bgnf_type)
            if bgnf_A + bgnf_G + bgnf_C + bgnf_T > 1 and bgnf_type is False:
                st.error(f"The sum of frequencies must equal 1.0. Current: {(bgnf_A + bgnf_G + bgnf_C + bgnf_T):.2f}")
            else:
                bgnf = {'A': bgnf_A, 'T': bgnf_T, 'G': bgnf_G, 'C': bgnf_C}

            st.divider()

            st.markdown("🔹 :blue[**_Experimental_**] Analyse all directions",
                        help='Directions: **original (+ →)**, **reverse-complement (- ←)**, reverse (+ ←), complement (- →)\n\n'
                             'Directions in bold are the default directions.')
            alldirection = st.toggle('All directions')
            if alldirection:
                st.markdown(
                    '⚠️Analyzes in the reverse (+ ←) and complement (- →) directions are generally not suitable for studying TFBS.')
                analyse = 4
            else:
                analyse = 2

    if tss_ge_input != 0:
        tss_ge_distance = int(tss_ge_input)
    else:
        tss_ge_distance = None

    threshold = float(threshold_entry)
    if jaspar == 'JASPAR_ID':
        pass
    else:
        if not isUIPAC or not error_input_im or isfasta or bgnf_A + bgnf_G + bgnf_C + bgnf_T > 1 and bgnf_type is False:
            button = True
            if not isUIPAC:
                st.error("Please use IUPAC code for Responsive Elements")
            elif isfasta:
                st.error("Please use only A, T, G, C, N in your sequence")
        else:
            button = False

    sequence_iteration = analyse * total_sequences_region_length
    num_random_seqs = 1000000
    if total_sequences <= 10:
        random_gen = total_sequences * num_random_seqs
    else:
        random_gen = num_random_seqs
    random_score = random_gen * analyse

    if pvalue:
        iteration = sequence_iteration + random_gen + random_score
    else:
        iteration = sequence_iteration

    st.markdown("")
    if st.button("🔹 :blue[**Step 2.6**] Click here to find motif in your sequences 🔎 🧬",
                 use_container_width=True,
                 disabled=button):

        with stqdm(total=iteration,
                   desc='**:blue[Analyse sequence...] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**',
                   mininterval=0.1) as progress_bar:
            individual_motif_occurrences, message = IMO.individual_motif_finder(
                dna_sequences, threshold, matrix, progress_bar, calc_pvalue, tss_ge_distance, alldirection,
                pseudocount if pseudocount_type is False else None,
                bgnf if bgnf_type is False else None)
        if message is True:
            st.session_state['individual_motif_occurrences'] = individual_motif_occurrences
            st.session_state['message'] = message
        elif message is False:
            st.error(individual_motif_occurrences)

    st.divider()
    try:
        if 'individual_motif_occurrences' in st.session_state:
            if len(st.session_state['individual_motif_occurrences']) > 0:
                current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                st.subheader(':blue[Results]')

                st.markdown('**Table**')
                tablecol1, tablecol2 = st.columns([0.75, 0.25])
                with tablecol1:
                    st.dataframe(st.session_state['individual_motif_occurrences'], hide_index=True)
                    csv_file = st.session_state['individual_motif_occurrences'].to_csv(index=False)
                    excel_file = io.BytesIO()
                    st.session_state['individual_motif_occurrences'].to_excel(excel_file, index=False,
                                                                              sheet_name='Sheet1')
                    excel_file.seek(0)

                with tablecol2:
                    st.success(f"Finding responsive elements done !")

                st.markdown("")
                st.markdown('**Graph**',
                            help='Zoom +/- with the mouse wheel. Drag while pressing the mouse to move the graph. Selection of a group by clicking on a point of the graph (double click de-selection). Double-click on a point to reset the zoom and the moving of graph.')

                result_table_output(st.session_state['individual_motif_occurrences'])

                with tablecol2:
                    st.download_button("💾 Download table (.xlsx)", excel_file,
                                       file_name=f'Results_TFinder_{current_date_time}.xlsx',
                                       mime="application/vnd.ms-excel", key='download-excel')
                    st.download_button(label="💾 Download table (.csv)", data=csv_file,
                                       file_name=f"Results_TFinder_{current_date_time}.csv", mime="text/csv")

                    try:
                        if st.session_state["LOCAL"] == "False":
                            email_receiver = st.text_input('Send results by email ✉',
                                                           value='', placeholder='Send results by email ✉',
                                                           label_visibility="collapsed")
                            if st.button("Send ✉"):
                                if jaspar == 'PWM':
                                    if matrix_type == 'With PWM':
                                        body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nPosition Weight Matrix:\n{matrix_str}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                                    if matrix_type == 'With FASTA sequences':
                                        body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nResponsive Elements:\n{individual_motif}\n\nPosition Weight Matrix:\n{matrix_text}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                                elif jaspar == 'JASPAR_ID':
                                    body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nJASPAR_ID: {jaspar_id} | Transcription Factor name: {TF_name}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                                else:
                                    body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nResponsive Elements:\n{IUPAC}\n\nPosition Weight Matrix:\n{matrix_text}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                                email(excel_file, csv_file, txt_output, email_receiver, body, jaspar)
                    except Exception:
                        print("You are in LOCAL")
            else:
                st.error(f"No consensus sequence found with the specified threshold")
    except Exception as e:
        print(e)
        st.write(e)
