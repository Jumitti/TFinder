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
import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

import altair as alt
import streamlit as st


class ResultDisplayExport:
    def __init__(self):
        ok = ok
    @staticmethod
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
                    image = MIMEImage(st.session_state['weblogo'].read(), name=f'LOGOMAKER_{current_date_time}.jpg')
                    msg.attach(image)
            elif jaspar == 'Manual sequence':
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

    @staticmethod
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
            tooltip=['Position'] + (['Rel Position'] if "Rel Position" in source else []) + ['Rel Score'] + (
                ['p-value'] if 'p-value' in source else []) + ['Sequence', 'Gene', 'Species', 'Region'],
            opacity=alt.condition(gene_region_selection, alt.value(0.8), alt.value(0.2))
        ).transform_calculate(x=f'datum[{xcol_param.name}]').properties(width=600,
                                                                        height=400).interactive().add_params(
            gene_region_selection, xcol_param)
        st.altair_chart(chart, theme=None, use_container_width=True)
