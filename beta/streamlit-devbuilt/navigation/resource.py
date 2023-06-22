import streamlit as st
from streamlit_pdf_viewer import pdf_viewer

def resource_page():
    st.success('resource')
    pdf_viewer("Promoter_finder_HELP.pdf")