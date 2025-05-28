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
from streamlit_modal import Modal
import streamlit_lottie
import time
import json

from navigation.allapp import allapp_page
from navigation.contact import contact_page
from navigation.home import home_page
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
CONTACT = 'Contact'

tabs = [
    HOME,
    APPLICATION,
    CONTACT,
]

option_data = [
    {'icon': "🏠", 'label': HOME},
    {'icon': "🖥️", 'label': APPLICATION},
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

if 'LOCAL' not in st.session_state:
    local_test = platform.processor()

    if local_test == "":
        st.session_state["LOCAL"] = 'False'
    else:
        print("Platform:", local_test)
        st.session_state["LOCAL"] = 'True'

if st.session_state["LOCAL"] == 'False':
    if st.secrets["message_from_god"] != "":
        st.warning(st.secrets["message_from_god"])
    if st.secrets['ncbi_error'] == "True":
        st.error("⚠ NCBI server maintenance, problems and slowdowns may be observed")

if chosen_tab == HOME:
    home_page()

elif chosen_tab == APPLICATION:
    allapp_page()

elif chosen_tab == CONTACT:
    contact_page()

for i in range(4):
    st.markdown('#')
st.markdown(footer, unsafe_allow_html=True)

# Credit
st.logo("img/TFinder_logo_site.png")
st.sidebar.image("img/TFinder_logo_site.png")

st.sidebar.title("Other links")
st.sidebar.markdown("If the application does not work, here are other deployments:\n"
                    f"   - TFinder on [Streamlit](https://streamlit.io/): [https://tfinder-ipmc.streamlit.app/](https://tfinder-ipmc.streamlit.app/)\n"
                    f"   - TFinder on [Health Universe](https://www.healthuniverse.com/): [https://apps.healthuniverse.com/pmb-xci-tsb](https://apps.healthuniverse.com/pmb-xci-tsb)\n"
                    f"   - (BETA) TFinder: [https://tfinder-beta.streamlit.app/](https://tfinder-beta.streamlit.app/)\n")

# Help
st.sidebar.title("Help")
st.sidebar.markdown("FULL DOCUMENTATION [HERE](https://jumitti.notion.site/tfinder?pvs=4)")

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

st.sidebar.title('Disclaimer')
st.sidebar.markdown(
    'TFinder use [NCBI API](https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen)'
    ': More information [NCBI Website and Data Usage Policies and Disclaimers](https://www.ncbi.nlm.nih.gov/home/about/policies/)')
if st.session_state['LOCAL'] == 'False':
    if st.secrets['ncbi_error'] == "True":
        st.sidebar.error("⚠ NCBI server maintenance, problems and slowdowns may be observed")
st.sidebar.markdown("TFinder use [JASPAR API](https://doi.org/10.1093/bioinformatics/btx804)")
st.sidebar.markdown('')
st.sidebar.markdown(
    'If you encounter a problem, please send an email to minniti@ipmc.cnrs.fr or minnitijulien06@gmail.com or use the [Issues](https://github.com/Jumitti/TFinder/issues) tab on GitHub')
st.sidebar.markdown(
    'Links are also available at the bottom of the left sidebar. You can contact us using the “Contact” tab too.')
st.sidebar.markdown("By using TFinder, you agree with data usage polices of NCBI and JASPAR")

st.sidebar.title("More")
st.sidebar.markdown(
    "[Report a bug 🐞](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=bug&projects=&template=bug_report.md&title=%5BBUG%5D)")
st.sidebar.markdown(
    "[Need HELP 🆘](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=help+wanted&projects=&template=help.md&title=%5BHELP%5D)")
st.sidebar.markdown(
    "[Have a question 🤔](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=question&projects=&template=question_report.md&title=%5BQUESTION%5D)")
st.sidebar.markdown(
    "[Features request 💡](https://github.com/Jumitti/TFinder/issues/new?assignees=&labels=enhancement&projects=&template=feature_request.md&title=%5BFEATURE%5D)")

if st.session_state["LOCAL"] == 'True':
    st.sidebar.markdown(f"TFinder Local Version")
else:
    unique_users = st.secrets['unique_users']
    st.sidebar.markdown(f"Unique users 👥: {unique_users}")
