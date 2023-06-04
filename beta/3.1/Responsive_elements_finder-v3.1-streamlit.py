import streamlit as st

#StreamLit

st.title('Responsive Elements Finder')

# Promoter Finder STREAMLIT

st.header('Promoter Finder')

# Gene ID
gene_id_entry = st.text_area("Gene ID:")

# Species
st.text("Species :")
species_combobox = st.selectbox("", ["Human", "Mouse", "Rat"], index=0)

# Responsive ELements Finder
st.header('Responsive Elements Finder')