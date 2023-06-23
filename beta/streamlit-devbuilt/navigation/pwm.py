import json
import streamlit as st
import ast

def pwm_page():
    # Obtenir la matrice depuis la zone de texte
    matrix_text = st.text_area("Saisir la matrice")

    # Convertir la chaîne de caractères en une structure de données Python
    try:
        matrix = ast.literal_eval(matrix_text)
    except (SyntaxError, ValueError):
        st.error("Erreur lors de l'analyse de la matrice. Assurez-vous qu'elle est correctement formatée.")
        matrix = None

    if matrix:
        # Conversion au format JSON
        json_matrix = json.dumps(matrix, indent=4)

        # Affichage du résultat
        st.code(json_matrix, language='json')
