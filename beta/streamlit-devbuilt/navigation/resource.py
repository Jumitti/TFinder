import streamlit as st
import requests

def resource_page():
    url = "https://github.com/Jumitti/Responsive-Elements-Finder/blob/main/Promoter_finder_HELP.pdf"
    response = requests.get(url)

    # Vérifier que la requête a réussi
    if response.status_code == 200:
        # Enregistrer le fichier PDF localement
        with open("fichier.pdf", "wb") as f:
            f.write(response.content)
    st.success('resource')