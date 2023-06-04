import streamlit as st
import clipboard

#StreamLit

st.title('Responsive Elements Finder')

# Promoter Finder STREAMLIT

st.header('Promoter Finder')

# Gene ID
gene_id_entry = st.text_area("Gene ID:")

# Species
st.text("Species:")

# Sélection de l'espèce dans le menu déroulant
selected_species = st.selectbox("", ["Human", "Mouse", "Rat"], index=0)

# Affichage de l'espèce sélectionnée
st.write("Species :", selected_species)



# Responsive ELements Finder
st.header('Responsive Elements Finder')