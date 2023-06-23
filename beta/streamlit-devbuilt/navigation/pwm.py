import streamlit as st
import numpy as np
import json

def pwm_page():
    # Obtenir la matrice depuis la zone de texte
    matrix_text = st.text_area("Saisir la matrice")

    # Convertir la chaîne de caractères en une structure de données Python
    matrix = eval(matrix_text)

    # Conversion au format JSON
    json_matrix = json.dumps(matrix, indent=4)

    # Affichage du résultat
    st.code(json_matrix, language='json')

