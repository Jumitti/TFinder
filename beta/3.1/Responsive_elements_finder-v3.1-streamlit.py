import streamlit as st

#StreamLit

st.title('Responsive Elements Finder')

# Promoter Finder STREAMLIT

st.header('Promoter Finder')

# Gene ID
gene_id_entry = st.text_area("Gene ID:", value="PRKN \n 5071")

# Species
species_combobox = st.selectbox("Species:", ["Human", "Mouse", "Rat"], index=0)

#Upstream
upstream_entry = st.text_input("Upstream:", value="2000")

#Downstream
downstream_entry = st.text_input("Downstream:", value="500")

#Run Promoter Finder
if st.button("Find promoter (~30sec/gene)"):
    get_sequence()
    
#Promoter
result_promoter = st.text_area("Promoter:")

# Responsive ELements Finder
st.header('Responsive Elements Finder')