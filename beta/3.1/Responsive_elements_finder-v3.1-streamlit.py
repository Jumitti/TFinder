import streamlit as st
import clipboard

#StreamLit

st.title('Responsive Elements Finder')

# Promoter Finder STREAMLIT

st.header('Promoter Finder')

# Gene ID
gene_id_entry = st.text_area("Gene ID:")

if st.button("Paste"):
    clipboard_content = clipboard.paste()

    gene_id_entry = st.text_area("Gene ID:", value=clipboard_content)

# Species



# Responsive ELements Finder
st.header('Responsive Elements Finder')