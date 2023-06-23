import json
import streamlit as st

# Obtenir la matrice depuis la zone de texte
matrix_text = st.text_area("Saisir la matrice")

# Convertir la chaîne de caractères en une structure de données Python
matrix = eval(matrix_text)

# Conversion au format JSON
json_matrix = json.dumps(matrix, indent=4)

# Affichage du résultat
st.code(json_matrix, language='json')
