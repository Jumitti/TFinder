import streamlit as st

#StreamLit

st.title('Responsive Elements Finder')

# Promoter Finder STREAMLIT
st.header('Promoter Finder')

# Gene ID
if st.button("Paste"):
    clipboard_content = st.paste()

    gene_id_entry = st.text_area("Gene ID:", value=clipboard_content)

# Responsive ELements Finder
st.header('Responsive Elements Finder')