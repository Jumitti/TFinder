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

st.set_page_config(
        page_title='REF by Minniti',
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
RESOURCE = 'Resource'
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
    st.write("How to extract promoter and find responsive elements")
    st.video('https://www.youtube.com/watch?v=lknbKbZCXuo')
    st.write("How to use FASTA sequences and find responsive elements")
    st.video('https://www.youtube.com/watch?v=QelVLLuNJqs')
    st.write("How to use JASPAR option")
    st.video('https://www.youtube.com/watch?v=DH8PBVqa860')
    
with st.sidebar.expander("Promoter & Terminator Extractor"):
    st.subheader("Gene ID:")
    st.write("ENTREZ_GENE_ID of NCBI and gene names are allowed.")
    st.write("There is no limit to the number of gene names/ENTREZ_GENE_ID. Add them with a line break (like those displayed by default). You can mix ENTREZ_GENE_ID and gene names as long as they are of the same species.")
    st.subheader("Species:")
    st.write("Human, mouse, rat, drosophila and zebrafish are allowed.")
    st.write("If you use several ENTREZ_GENE_ID/gene names, make sure you select the correct species.")
    st.subheader("Upstream/Downstream:")
    st.write("Distance to Transcription Start Site (TSS) in bp.")
    st.image("https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/whatisagene.png")
    st.subheader("Promoter & Terminator:")
    st.write('Use "Find promoter/extractor" button or paste your sequences. FASTA format allowed and required for multiple sequences.')
    st.write('FASTA format: All sequences must have the TSS at the same distance, otherwise you assume the inconsistency of the positions of found sequences')
    
with st.sidebar.expander("Binding Sites Finder"):
    st.subheader("Responsive element:")
    st.write("To use the JASPAR option, check the box and use the JASPAR_ID of your transcription factor.")
    st.write('If you want to use your responsive element, do not check the JASPAR option.')
    st.write('IUPAC allowed')
    st.image("https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/IUPAC.png")
    st.subheader("Transcription Start Site (TSS):")
    st.write('Distance to Transcription Start Site (TSS) in bp')
    st.write('Note: If you use Step 1 , it will be defined automatically.')
    st.subheader("Threshold:")
    st.write('Eliminates responsive element with homology < threshold or score < threshold')
    st.write('Note for JASPAR option: Score is normalized to the maximum PWM score of the requested transcription factor. The result is displayed as a percentage')
    st.write('Note without JASPAR option: Homology is calculated between the responsive element in the promoter and the responsive element requested. The calculation uses the Hamming distance, counts the number of differences and gives a percentage score homology.')
    
