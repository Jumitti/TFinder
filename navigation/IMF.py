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
import streamlit as st
from stqdm import stqdm

from tfinder import IMO


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

        if jaspar == 'PWM':
            if matrix_type == 'With FASTA sequences':
                image = MIMEImage(st.session_state['weblogo'].read(), name=f'LOGOMAKER_{current_date_time}.jpg')
                msg.attach(image)
        elif jaspar == 'Individual Motif':
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


def result_table_output(df):
    source = df
    score_range = source['Rel Score'].astype(float)
    ystart = score_range.min() - 0.02
    ystop = score_range.max() + 0.02
    source['Gene_Region'] = source['Gene'] + " " + source['Species'] + " " + source['Region']
    source['Beginning of sequences'] = source['Position']
    if 'Rel Position' in source:
        source['From TSS/gene end'] = source['Rel Position']
    scale = alt.Scale(scheme='category10')
    color_scale = alt.Color("Gene_Region:N", scale=scale)

    gene_region_selection = alt.selection_point(fields=['Gene_Region'], on='click', bind='legend')

    dropdown = alt.binding_select(
        options=['Beginning of sequences', 'From TSS/gene end' if "Rel Position" in source else []],
        name='(X-axis) Position from: ')

    xcol_param = alt.param(value='Beginning of sequences', bind=dropdown)

    chart = alt.Chart(source).mark_circle().encode(
        x=alt.X('x:Q').title('Position (bp)'),
        y=alt.Y('Rel Score:Q', axis=alt.Axis(title='Relative Score'),
                scale=alt.Scale(domain=[ystart, ystop])),
        color=alt.condition(gene_region_selection, color_scale, alt.value('lightgray')),
        tooltip=['Sequence', 'Position'] + (['Rel Position'] if "Rel Position" in source else []) + ['Rel Score'] + (
            ['p-value'] if 'p-value' in source else []) + ['Gene', 'Species', 'Region'],
        opacity=alt.condition(gene_region_selection, alt.value(0.8), alt.value(0.2))
    ).transform_calculate(x=f'datum[{xcol_param.name}]').properties(width=600,
                                                                    height=400).interactive().add_params(
        gene_region_selection, xcol_param)
    st.altair_chart(chart, theme=None, use_container_width=True)


def BSF_page():
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
        dna_sequences.append((name, dna_sequence, species, region))
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
                dna_sequences.append((name, dna_sequence, found_species, region))
                i += 1
            else:
                i += 1
    elif not lines.startswith(("A", "T", "C", "G", "N", "a", "t", "c", "g", "n", "I", "i", "")):
        isfasta = True
    else:
        isfasta = False

    total_sequences_region_length = sum(len(dna_sequence) for _, dna_sequence, _, _ in dna_sequences)
    total_sequences = len(dna_sequences)

    # RE entry
    REcol1, REcol2 = st.columns([0.30, 0.70])
    with REcol1:
        st.markdown('🔹 :blue[**Step 2.2**] Responsive elements type:')
        jaspar = st.radio('🔹 :blue[**Step 2.2**] Responsive elements type:',
                          ('Individual Motif', 'JASPAR_ID', 'PWM'),
                          label_visibility='collapsed')
    if jaspar == 'JASPAR_ID':
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
                st.markdown("🔹 :blue[**Step 2.3**] Matrix:",
                            help="Only PWM generated with our tools are allowed")
                matrix_str = st.text_area("🔹 :blue[**Step 2.3**] Matrix:",
                                          value="A [ 20.0 0.0 0.0 0.0 0.0 0.0 0.0 100.0 0.0 60.0 20.0 ]\nT [ 60.0 20.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]\nG [ 0.0 20.0 100.0 0.0 0.0 100.0 100.0 0.0 100.0 40.0 0.0 ]\nC [ 20.0 60.0 0.0 100.0 100.0 0.0 0.0 0.0 0.0 0.0 80.0 ]"
                                          if 'MATRIX_STR_save' not in st.session_state else st.session_state['MATRIX_STR_save'],
                                          label_visibility='collapsed', height=125)
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
                                                if 'individual_motif_save' not in st.session_state else st.session_state['individual_motif_save'],
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
            st.markdown("🔹 :blue[**Step 2.3**] Responsive element:", help="IUPAC authorized")
            IUPAC = st.text_input("🔹 :blue[**Step 2.3**] Responsive element (IUPAC authorized):",
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
    BSFcol1, BSFcol2, BSFcol3 = st.columns([2, 2, 2], gap="medium")
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
        if auto_thre:
            threshold_entry = 0
        else:
            threshold_entry = st.slider("🔹 :blue[**Step 2.5**] Relative Score threshold", 0.5, 1.0, 0.85,
                                        step=0.05,
                                        label_visibility="collapsed")
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

    if tss_ge_input != 0:
        tss_ge_distance = int(tss_ge_input)
    else:
        tss_ge_distance = None

    threshold = float(threshold_entry)
    if jaspar == 'JASPAR_ID':
        pass
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
            button = False

    sequence_iteration = 4 * total_sequences_region_length
    num_random_seqs = 1000000
    if total_sequences <= 10:
        random_gen = total_sequences * num_random_seqs
    else:
        random_gen = num_random_seqs
    random_score = random_gen * 4

    if pvalue:
        iteration = sequence_iteration + random_gen + random_score
    else:
        iteration = sequence_iteration

    st.markdown("")
    if st.button("🔹 :blue[**Step 2.6**] Click here to find motif in your sequences 🔎 🧬",
                 use_container_width=True,
                 disabled=button):
        with stqdm(total=iteration,
                   desc='**:blue[Extract sequence...] ⚠️:red[PLEASE WAIT UNTIL END WITHOUT CHANGING ANYTHING]**',
                   mininterval=0.1) as progress_bar:
            individual_motif_occurrences = IMO.individual_motif_finder(dna_sequences, threshold, matrix,
                                                                       progress_bar,
                                                                       calc_pvalue,
                                                                       tss_ge_distance)
        st.session_state['individual_motif_occurrences'] = individual_motif_occurrences

    st.divider()
    if 'individual_motif_occurrences' in st.session_state:
        if len(st.session_state['individual_motif_occurrences']) > 1:
            current_date_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            st.subheader(':blue[Results]')

            df = pd.DataFrame(st.session_state['individual_motif_occurrences'][1:],
                              columns=st.session_state['individual_motif_occurrences'][0])
            st.session_state['df'] = df

            st.markdown('**Table**')
            tablecol1, tablecol2 = st.columns([0.75, 0.25])
            with tablecol1:
                st.dataframe(df, hide_index=True)
                csv_file = df.to_csv(index=False)
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
                st.download_button(label="💾 Download table (.csv)", data=csv_file,
                                   file_name=f"Results_TFinder_{current_date_time}.csv", mime="text/csv")

                if st.session_state["LOCAL"] == "False":
                    email_receiver = st.text_input('Send results by email ✉',
                                                   value='', placeholder='Send results by email ✉',
                                                   label_visibility="collapsed")
                    if st.button("Send ✉"):
                        if jaspar == 'PWM':
                            if matrix_type == 'With PWM':
                                body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nPosition Weight Matrix:\n{matrix_text}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                            if matrix_type == 'With FASTA sequences':
                                body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nResponsive Elements:\n{individual_motif}\n\nPosition Weight Matrix:\n{matrix_text}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                        elif jaspar == 'JASPAR_ID':
                            body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nJASPAR_ID: {jaspar_id} | Transcription Factor name: {TF_name}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                        else:
                            body = f"Hello 🧬\n\nResults obtained with TFinder.\n\nResponsive Elements:\n{IUPAC}\n\nPosition Weight Matrix:\n{matrix_text}\n\nThis email also includes the sequences used in FASTA format and an Excel table of results.\n\nFor all requests/information, please refer to the 'Contact' tab on the TFinder website. We would be happy to answer all your questions.\n\nBest regards\nTFinder Team 🔎🧬"
                        email(excel_file, csv_file, txt_output, email_receiver, body, jaspar)
        else:
            st.error(f"No consensus sequence found with the specified threshold")
