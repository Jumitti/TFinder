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

import streamlit as st
import hydralit_components as hc
import requests
import pandas as pd
import pickle
from utils.components import footer_style, footer
from navigation.home import home_page
from navigation.resource import resource_page
from navigation.contact import contact_page
from navigation.allapp import allapp_page
import streamlit_analytics
from streamlit_modal import st_modal

st.set_page_config(
    page_title='TFinder by Minniti Julien',
    page_icon="./.streamlit/REF.png",
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

# Footer

st.markdown(footer_style, unsafe_allow_html=True)

# NavBar

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
    {'icon': "🏠", 'label': HOME},
    {'icon': "🖥️", 'label': APPLICATION},
    {'icon': "📑", 'label': RESOURCE},
    {'icon': "✉️", 'label': CONTACT}

]

over_theme = {'txc_inactive': 'black', 'menu_background': '#D6E5FA', 'txc_active': 'white', 'option_active': '#749BC2'}
font_fmt = {'font-class': 'h3', 'font-size': '50%'}

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
st.markdown(footer, unsafe_allow_html=True)

streamlit_analytics.start_tracking()

# Credit rating
st.sidebar.image("https://raw.githubusercontent.com/Jumitti/TFinder/main/img/REF.png")
try:
    with open("ratings.pkl", "rb") as file:
        ratings = pickle.load(file)
except FileNotFoundError:
    ratings = []
rating = st.sidebar.slider("Rate it 😊 (1-5 ⭐)", 1, 5, 5)
colrate1, colrate2 = st.sidebar.columns(2)
with colrate1:
    submit_button = st.button("Submit Rating")
if submit_button:
    st.balloons()
    ratings.append(rating)
    with open("ratings.pkl", "wb") as file:
        pickle.dump(ratings, file)
    st.toast("Thank you for rating the application!", icon='😍')
    st.balloons()
average_rating = sum(ratings) / len(ratings) if ratings else 0
num_ratings = len(ratings)
with colrate2:
    st.write(f"{average_rating:.2f} ⭐ ({num_ratings} votes)")

# Help
st.sidebar.title("Help")
with st.sidebar.expander("Video tutorials"):
    st.write('coming soon')

with st.sidebar.expander("Promoter & Terminator Extractor"):
    st.subheader("Gene ID:")
    st.write("ENTREZ_GENE_ID of NCBI and gene names are allowed.")
    st.write(
        "There is no limit to the number of gene names/ENTREZ_GENE_ID. Add them with a line break (like those displayed by default). You can mix ENTREZ_GENE_ID and gene names as long as they are of the same species.")
    st.subheader("Species:")
    st.write("Human, mouse, rat, drosophila and zebrafish are allowed.")
    st.write("If you use several ENTREZ_GENE_ID/gene names, make sure you select the correct species.")
    st.subheader("Promoter/Terminator and Upstream/Downstream")
    st.write("Distance to Transcription Start Site (TSS) or gene end in bp.")
    st.image("https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/whatisagene.png")
    st.subheader("Promoter & Terminator:")
    st.write(
        'Use "Find promoter/extractor" button or paste your sequences. FASTA format allowed and required for multiple sequences.')
    st.write(
        'FASTA format: All sequences must have the TSS at the same distance, otherwise you assume the inconsistency of the positions of found sequences')

with st.sidebar.expander("Binding Sites Finder"):
    st.subheader("Responsive element:")
    st.write('IUPAC code is authorized for manual sequences')
    st.write(
        'You can generate a PWM with several sequences in FASTA format (same lenght required) or use a PWM already generated with our tools')
    st.write("For JASPAR option, use the JASPAR_ID of your transcription factor.")
    st.image("https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/IUPAC.png")
    st.subheader("Transcription Start Site (TSS) or gene end:")
    st.write('Distance to Transcription Start Site (TSS) or gene end in bp')
    st.write('Note: If you use Step 1 , it will be defined automatically.')
    st.subheader("Relative Score Threshold:")
    st.write('Eliminates responsive element with Relative Score < threshold')
    st.write(
        'The Relative Score represents the Score calculated for each k-mer of the length of the PWM in the given sequence where each corresponding probability is added according to each nucleotide. This Score is then normalized to the maximum and minimum PWM Score.')
    st.subheader('_p-value_')
    st.write(
        'The p-value calculation takes time so it is optional. it represents the probability that a random generated sequence of the lenght of the PWM with the nucleotide proportions of the sequence has a score greater than or equal to the element found.')

st.sidebar.title("Servers status",
                 help='✅: servers are reachable. You can use extract regions via NCBI/use the JASPAR_IDs\n\n❌: servers are unreachable. You can still use TFinder if you have a sequence in FASTA format and a pattern to search in the sequence')

if st.sidebar.button("Check"):
    response = requests.get(
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=nos2[Gene%20Name]+AND+human[Organism]&retmode=json&rettype=xml')
    response1 = requests.get('https://jaspar.genereg.net/api/v1/matrix/MA0106.1')

    ncbi_status = "✅" if response.status_code == 200 else "❌"
    jaspar_status = "✅" if response1.status_code == 200 else "❌"

    st.session_state['ncbi_status'] = ncbi_status
    st.session_state['jaspar_status'] = jaspar_status

    data = {
        "NCBI": [ncbi_status],
        "JASPAR": [jaspar_status]
    }

    df = pd.DataFrame(data, index=["Servers status"])

    st.sidebar.table(df)

st.sidebar.title("More")
st.sidebar.markdown("[Report a bug 🐞](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=bug&projects=&template=bug_report.md&title=%5BBUG%5D)")
st.sidebar.markdown("[Need HELP 🆘](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=help+wanted&projects=&template=help.md&title=%5BHELP%5D)")
st.sidebar.markdown("[Have a question 🤔](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=question&projects=&template=question_report.md&title=%5BQUESTION%5D)")
st.sidebar.markdown("[Features request 💡](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=enhancement&projects=&template=feature_request.md&title=%5BFEATURE%5D)")
st.sidebar.markdown("[Want to talk ? 🙋🏼‍♂](https://github.com/Jumitti/TFinder/discussions)")

streamlit_analytics.stop_tracking()

views = streamlit_analytics.main.counts["total_pageviews"]
st.sidebar.markdown(f"Total users 👥: {int(views)}")

st.sidebar.markdown('TFinder use NCBI API and JASPAR ID. See Resources for more information')

with st_modal():
    st.write("Contenu de la pop-up")
    st.button("Fermer")
