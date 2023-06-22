import streamlit as st
from PIL import Image
import datetime

TITLE = 'COVID-19db linkage maps of cell surface proteins and transcription factors in immune cells'

def home_page():
    
    _, middle, _ = st.columns([0.1, 1, 0.1])
    
    with middle:
        # st.markdown('### ' + TITLE) 
        st.markdown(f"<h3 style='text-align: center; color: black;'>{TITLE}</h1>", unsafe_allow_html=True)  
        
        st.markdown('👉 Read our paper here: https://www.biorxiv.org/content/10.1101/2022.12.14.520411v1')
    