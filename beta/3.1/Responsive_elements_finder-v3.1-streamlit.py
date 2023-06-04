import streamlit as st

#StreamLit

st.title('Responsive Elements Finder')

# Promoter Finder STREAMLIT

st.header('Promoter Finder')

# Gene ID
gene_id_entry = st.text_area("Gene ID:")

# Species
species_combobox = st.selectbox("Species:", ["Human", "Mouse", "Rat"], index=0)

#Upstream

upstream_entry = st.text_input("Upstream:", value="2000")

#Downstream

# Responsive ELements Finder
st.header('Responsive Elements Finder')