# Copyright (c) 2023 Minniti Julien

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import streamlit as st
import hydralit_components as hc
import requests
import pandas as pd
import altair as alt
import math
import pickle

from utils.components import footer_style, footer
from navigation.home import home_page
from navigation.resource import resource_page
from navigation.contact import contact_page
from navigation.allapp import allapp_page


import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

st.set_page_config(
        page_title='TFinder by Minniti Julien',
        initial_sidebar_state="expanded"
)

max_width_str = f"max-width: {75}%;"

st.markdown(f"""
        <style>
        .appview-container .main .block-container{{{max_width_str}}}
        </style>
        """,
        unsafe_allow_html=True,
    )

st.markdown("""
        <style>
               .block-container {
                    padding-top: 0rem;
                    padding-bottom: 0rem;
                    
                }
        </style>
        """, unsafe_allow_html=True)       

#Footer

st.markdown(footer_style, unsafe_allow_html=True) ## Footer

#NavBar

HOME = 'Home'
APPLICATION = 'Tools/Software'
RESOURCE = 'Resources'
CONTACT = 'Contact'

tabs = [
    HOME,
    APPLICATION,
    RESOURCE,
    CONTACT
]

option_data = [
   {'icon': "🏠", 'label':HOME},
   {'icon': "🖥️", 'label':APPLICATION},
   {'icon':"📑",'label':RESOURCE},
   {'icon':"✉️",'label':CONTACT}
   
]

over_theme = {'txc_inactive': 'black','menu_background':'#ECECEC','txc_active':'white','option_active':'#fab947'}
font_fmt = {'font-class':'h3','font-size':'50%'}

chosen_tab = hc.option_bar(
    option_definition=option_data,
    title='',
    key='PrimaryOptionx',
    override_theme=over_theme,
    horizontal_orientation=True)

if chosen_tab == HOME:
    home_page()

elif chosen_tab == APPLICATION: 
    allapp_page() 
    
elif chosen_tab == RESOURCE: 
    resource_page()
    
elif chosen_tab == CONTACT: 
    contact_page()
    
for i in range(4):
    st.markdown('#')
st.markdown(footer,unsafe_allow_html=True)

# Credit Eastereggs
try:
    with open("ratings.pkl", "rb") as file:
        ratings = pickle.load(file)
except FileNotFoundError:
    ratings = []
rating = st.sidebar.slider("Rate it 😊 (1-5 ⭐)", 1, 5, 5)
submit_button = st.sidebar.button("Submit Rating")
if submit_button:
    ratings.append(rating)
    with open("ratings.pkl", "wb") as file:
        pickle.dump(ratings, file)
    st.sidebar.success("Thank you for rating the application!")
average_rating = sum(ratings) / len(ratings) if ratings else 0
num_ratings = len(ratings)
st.sidebar.write(f"Average rating: {average_rating:.2f} ⭐ ({num_ratings} votes)")

#Help

st.sidebar.title("Help")
with st.sidebar.expander("Video tutorials"):
    st.write('coming soon')
    
with st.sidebar.expander("Promoter & Terminator Extractor"):
    st.subheader("Gene ID:")
    st.write("ENTREZ_GENE_ID of NCBI and gene names are allowed.")
    st.write("There is no limit to the number of gene names/ENTREZ_GENE_ID. Add them with a line break (like those displayed by default). You can mix ENTREZ_GENE_ID and gene names as long as they are of the same species.")
    st.subheader("Species:")
    st.write("Human, mouse, rat, drosophila and zebrafish are allowed.")
    st.write("If you use several ENTREZ_GENE_ID/gene names, make sure you select the correct species.")
    st.subheader("Promoter/Terminator and Upstream/Downstream")
    st.write("Distance to Transcription Start Site (TSS) or gene end in bp.")
    st.image("https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/whatisagene.png")
    st.subheader("Promoter & Terminator:")
    st.write('Use "Find promoter/extractor" button or paste your sequences. FASTA format allowed and required for multiple sequences.')
    st.write('FASTA format: All sequences must have the TSS at the same distance, otherwise you assume the inconsistency of the positions of found sequences')
    
with st.sidebar.expander("Binding Sites Finder"):
    st.subheader("Responsive element:")
    st.write('IUPAC code is authorized for manual sequences')
    st.write('You can generate a PWM with several sequences in FASTA format (same lenght required) or use a PWM already generated with our tools')
    st.write("For JASPAR option, use the JASPAR_ID of your transcription factor.")
    st.image("https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/IUPAC.png")
    st.subheader("Transcription Start Site (TSS) or gene end:")
    st.write('Distance to Transcription Start Site (TSS) or gene end in bp')
    st.write('Note: If you use Step 1 , it will be defined automatically.')
    st.subheader("Relative Score Threshold:")
    st.write('Eliminates responsive element with Relative Score < threshold')
    st.write('Relative Score = (TFBS-min)/(max-min) where TFBS is the score of the element found, min is the minimum score obtainable with the PWM and max is the maximum score obtainable with the PWM.')
    st.subheader('_p-value_')
    st.write('The p-value calculation takes time so it is optional. it represents the probability that a random generated sequence of the lenght of the PWM with the nucleotide proportions of the sequence has a score greater than or equal to the element found.')

url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=nos2[Gene%20Name]+AND+human[Organism]&retmode=json&rettype=xml'
url1 = 'https://jaspar.genereg.net/api/v1/matrix/'
response = requests.get(url)
response1 = requests.get(url1)

ncbi_status = "✅" if response.status_code == 200 else "❌"
jaspar_status = "✅" if response1.status_code == 200 else "❌"

data = {
    "NCBI": [ncbi_status],
    "JASPAR": [jaspar_status]
}

st.sidebar.title("Servers status")
df = pd.DataFrame(data, index=["Servers status"])

st.sidebar.table(df)
st.sidebar.markdown('✅: servers are reachable. ',help='You can use extract regions via NCBI/use the JASPAR_IDs')
st.sidebar.markdown('❌: servers are unreachable. ',help='You can still use TFinder if you have a sequence in FASTA format and a pattern to search in the sequence')


# Fonction pour envoyer l'e-mail
def send_email(sender_email, sender_password, receiver_email, subject, body):
    try:
        # Configuration du serveur SMTP de Gmail
        server = smtplib.SMTP_SSL("smtp.gmail.com", 465)
        server.ehlo()
        server.login(sender_email, sender_password)
        st.swrite(login)

        # Création du message
        message = MIMEMultipart()
        message["From"] = sender_email
        message["To"] = receiver_email
        message["Subject"] = subject
        message.attach(MIMEText(body, "plain"))

        # Envoi de l'e-mail
        server.sendmail(sender_email, receiver_email, message.as_string())

        # Fermeture du serveur SMTP
        server.quit()

        st.success("E-mail envoyé avec succès !")
    except Exception as e:
        st.error(f"Erreur lors de l'envoi de l'e-mail : {e}")

# Interface utilisateur Streamlit
st.title("Envoyer un e-mail avec Streamlit")

sender_email = "noreply.results.tfinder@gmail.com"
sender_password = "PaulineJulien201097@" 
receiver_email = st.text_input("Adresse e-mail du destinataire")
subject = st.text_input("Sujet de l'e-mail")
body = st.text_area("Corps de l'e-mail", height=200)

if st.button("Envoyer l'e-mail"):
    send_email(sender_email, sender_password, receiver_email, subject, body)

