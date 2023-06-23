import json
import streamlit as st

def pwm_page():
    # Obtenir la matrice depuis la zone de texte
    matrix_text = st.text_area("Saisir la matrice")

    # Analyser la chaîne de caractères pour extraire les valeurs de la matrice
    matrix_lines = matrix_text.split('\n')
    matrix_data = {}
    for line in matrix_lines:
        line = line.strip()
        if line:
            key, values = line.split('[', 1)
            values = values.replace(']', '').split()
            values = [float(value) for value in values]
            matrix_data[key.strip()] = values

    # Conversion au format JSON
    json_matrix = json.dumps(matrix_data, indent=4)

    # Affichage du résultat
    st.code(json_matrix, language='json')