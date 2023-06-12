import streamlit as st
import requests
import pandas as pd

def find_transcription_factors(sequence, species):
    url = "https://jaspar.genereg.net/api/v1/matrix/search/"
    data = {
        "format": "json",
        "seq": sequence,
        "tax_id": species
    }
    response = requests.post(url, json=data)
    if response.status_code == 200:
        response_data = response.json()
        transcription_factors = [entry['collection'] for entry in response_data['results']]
        return transcription_factors
    else:
        return None

# Interface utilisateur avec Streamlit
st.title("Recherche de facteurs de transcription")

sequence = st.text_input("Saisissez la séquence d'ADN", "ATCGATCGATCG")
species = st.selectbox("Sélectionnez l'espèce", ["Homo sapiens", "Mus musculus", "Drosophila melanogaster"])

if st.button("Rechercher"):
    transcription_factors = find_transcription_factors(sequence, species)
    if transcription_factors:
        df = pd.DataFrame(transcription_factors, columns=["Facteur de transcription"])
        st.table(df)
    else:
        st.warning("Aucun facteur de transcription trouvé pour la séquence donnée et l'espèce spécifiée.")
