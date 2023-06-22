import streamlit as st
from PIL import Image
import datetime

def home_page():
    
    _, middle, _ = st.columns([1, 1, 0.1])
    
    with middle:
        st.header('Responsive Elements Finder 🧬🔎')
    