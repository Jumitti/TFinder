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

import hydralit_components as hc
import platform
import pandas as pd
import requests
import streamlit as st
# import streamlit_analytics
from streamlit_modal import Modal
import streamlit_lottie
import time
import json

from navigation.allapp import allapp_page
from navigation.contact import contact_page
from navigation.home import home_page
from navigation.resource import resource_page
from utils.components import footer_style, footer

try:
    from streamlit import rerun as rerun
except ImportError:
    # conditional import for streamlit version <1.27
    from streamlit import experimental_rerun as rerun

import os


def load_lottiefile(filepath: str):
    with open(filepath, "r") as f:
        return json.load(f)


st.set_page_config(
    page_title='TFinder by Minniti Julien',
    page_icon="img/TFinder_logo_page.png",
    initial_sidebar_state="expanded",
    layout="wide"
)

if 'lottie' not in st.session_state:
    st.session_state.lottie = False

if not st.session_state.lottie:
    lottfinder = load_lottiefile(".streamlit/TFinder_logo_animated.json")
    st.lottie(lottfinder, speed=1.3, loop=False, quality="high", height=1520, width=1520)
    time.sleep(2)
    st.session_state.lottie = True
    rerun()

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
# RESOURCE = 'Resources'
CONTACT = 'Contact'

tabs = [
    HOME,
    APPLICATION,
    # RESOURCE,
    CONTACT,
]

option_data = [
    {'icon': "🏠", 'label': HOME},
    {'icon': "🖥️", 'label': APPLICATION},
    # {'icon': "📑", 'label': RESOURCE},
    {'icon': "✉️", 'label': CONTACT},
]

over_theme = {'txc_inactive': 'black', 'menu_background': '#D6E5FA', 'txc_active': 'white', 'option_active': '#749BC2'}
font_fmt = {'font-class': 'h3', 'font-size': '50%'}

chosen_tab = hc.option_bar(
    option_definition=option_data,
    title='',
    key='PrimaryOptionx',
    override_theme=over_theme,
    horizontal_orientation=True)

st.success("Hello everyone, TFinder is growing every day and we would like to know you a little more. "
           "We will not collect any data through Streamlit and it is difficult for us to know your uses and your feedback. "
           f"[HERE](https://airtable.com/appRn3TQqhuSFS8KO/pagm4Vau8lEFdRX3q/form) you will find a form to answer some of our questions if you wish. See you soon 😊")

st.success("If the application does not work, here are other deployments:\n"
           f"   - TFinder on [Streamlit](https://streamlit.io/): [https://tfinder-ipmc.streamlit.app/](https://tfinder-ipmc.streamlit.app/)\n"
           f"   - TFinder on [Health Universe](https://www.healthuniverse.com/): [https://apps.healthuniverse.com/nhu-dxv-ktj](https://apps.healthuniverse.com/nhu-dxv-ktj)\n"
           f"   - (BETA) TFinder: [https://tfinder-beta.streamlit.app/](https://tfinder-beta.streamlit.app/)\n")

st.warning(
    "The Resources tab no longer exists, as the documentation has been (partially) updated with the new updates. I invite you to discover it in the panel on the left or [HERE](https://jumitti.notion.site/tfinder?pvs=4)")

if chosen_tab == HOME:
    home_page()

elif chosen_tab == APPLICATION:
    allapp_page()

# elif chosen_tab == RESOURCE:
#     resource_page()

elif chosen_tab == CONTACT:
    contact_page()

for i in range(4):
    st.markdown('#')
st.markdown(footer, unsafe_allow_html=True)

# streamlit_analytics.start_tracking()

# Credit
st.logo("img/TFinder_logo_site.png")
st.sidebar.image("img/TFinder_logo_site.png")

# Help
st.sidebar.title("Help")
# with st.sidebar.expander("Video tutorials"):
#     st.write('coming soon')
st.sidebar.markdown("FULL DOCUMENTATION [HERE](https://jumitti.notion.site/tfinder?pvs=4)")

with st.sidebar.expander("Regulatory regions extractor"):
    st.subheader("Gene ID:")
    st.write("ENTREZ_GENE_ID of NCBI and gene names are allowed.")
    st.write(
        "There is no limit to the number of gene names/ENTREZ_GENE_ID. Add them with a line break "
        "(like those displayed by default). You can mix ENTREZ_GENE_ID and gene names as long as they "
        "are of the same species.")
    st.write("**Advance mode** allows you to select multiple species for genes")
    st.write("⚠️A **Check genes avaibility** button allows you to analyse if your gene is accessible for species"
             "and if ID is correct. Please use it. ")

    st.subheader("Species:")
    st.write("Human, mouse, rat, drosophila and zebrafish are allowed.")
    st.write("If you use several ENTREZ_GENE_ID/gene names, make sure you select the correct species.")
    st.write("⚠️Use **Check genes avaibility** button for checking species")
    st.write("**Advance mode** allows you to select multiple species for genes")

    st.subheader("Regulatory regions and Upstream/Downstream")
    st.write("Distance to Transcription Start Site (TSS) or gene end in bp.")
    st.image("img/whatisagene.png")

    st.subheader("Sequences to analyse:")
    st.write(
        'Use "Find promoter/extractor" button or paste your sequences. FASTA format allowed and required for multiple sequences.')
    st.write(
        'FASTA format: All sequences must have the TSS at the same distance, otherwise you assume the inconsistency of the positions of found sequences')

with st.sidebar.expander("Individual Motif Finder"):
    st.subheader("Responsive element:")
    st.write('For **Individual Motif**: IUPAC code is authorized')
    st.write(
        'For **PWM**: You can generate a PWM with several sequences in FASTA format or use a PWM already generated with our tools  (same length required)')
    st.write("For **JASPAR_ID** option, use the JASPAR_ID of your transcription factor.")
    st.image("img/IUPAC.png")
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
    with st.sidebar:
        with st.spinner('Please wait...'):
            response = requests.get(
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=nos2[Gene%20Name]+AND+human[Organism]&retmode=json&rettype=xml')
            response1 = requests.get('https://jaspar.elixir.no/api/v1/matrix/MA0106.1')

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
st.sidebar.markdown(
    "[Report a bug 🐞](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=bug&projects=&template=bug_report.md&title=%5BBUG%5D)")
st.sidebar.markdown(
    "[Need HELP 🆘](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=help+wanted&projects=&template=help.md&title=%5BHELP%5D)")
st.sidebar.markdown(
    "[Have a question 🤔](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=question&projects=&template=question_report.md&title=%5BQUESTION%5D)")
st.sidebar.markdown(
    "[Features request 💡](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=enhancement&projects=&template=feature_request.md&title=%5BFEATURE%5D)")

# streamlit_analytics.stop_tracking()
# views = streamlit_analytics.main.counts["total_pageviews"]


modal = Modal(key="TFinder Key", title="Disclaimers - Welcome to TFinder", padding=50, max_width=900)

if 'popup_closed' not in st.session_state:
    st.session_state.popup_closed = False

if not st.session_state.popup_closed:
    with modal.container():
        st.markdown('')
        st.markdown(
            'TFinder use [NCBI API](https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen)'
            ': More information [NCBI Website and Data Usage Policies and Disclaimers](https://www.ncbi.nlm.nih.gov/home/about/policies/)')
        st.markdown("TFinder use [JASPAR API](https://doi.org/10.1093/bioinformatics/btx804)")
        st.markdown('')
        st.markdown(
            'If you encounter a problem, please send an email to minniti@ipmc.cnrs.fr or minnitijulien06@gmail.com or use the [Issues](https://github.com/Jumitti/TFinder/issues) tab on GitHub')
        st.markdown(
            'Links are also available at the bottom of the left sidebar. You can contact us using the “Contact” tab too.')
        value = st.checkbox("By checking this box, you agree with data usage polices of NCBI and JASPAR")
        if value:
            st.button('Close')
            st.session_state.popup_closed = True
