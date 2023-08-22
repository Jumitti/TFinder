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


def result_table_output(df, position_type):
    source = df
    score_range = source['Rel Score'].astype(float)
    ystart = score_range.min() - 0.02
    ystop = score_range.max() + 0.02
    source['Gene_Region'] = source['Gene'] + " " + source['Species'] + " " + source['Region']
    scale = alt.Scale(scheme='category10')
    color_scale = alt.Color("Gene_Region:N", scale=scale)
    gene_region_selection = alt.selection_point(fields=['Gene_Region'], on='click', bind='legend')

    if 'p-value' in source:
        ispvalue = True

    chart = alt.Chart(source).mark_circle().encode(
        x=alt.X('Rel Position:Q' if position_type == 'From TSS/gene end' else 'Position:Q',
                axis=alt.Axis(title='Relative position (bp)'), sort='ascending'),
        y=alt.Y('Rel Score:Q', axis=alt.Axis(title='Relative Score'),
                scale=alt.Scale(domain=[ystart, ystop])),
        color=alt.condition(gene_region_selection, color_scale, alt.value('lightgray')),
        tooltip=['Rel Position' if position_type == 'From TSS/gene end' else 'Position', 'Rel Score'] + (
            ['p-value'] if 'p-value' in source else []) + ['Sequence', 'Gene', 'Species', 'Region'],
        opacity=alt.condition(gene_region_selection, alt.value(0.8), alt.value(0.2))
    ).properties(width=600, height=400).interactive().add_params(gene_region_selection)
    st.altair_chart(chart, theme=None, use_container_width=True)