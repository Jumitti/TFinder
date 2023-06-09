import requests
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tabulate import tabulate

def search_sequence(jaspar_id, promoter_sequence, matrices):
    results = []
    max_scores = []

    for matrix_name, matrix in matrices.items():
        seq_length = len(matrix['A'])

        # Calcul du score maximum pour chaque position
        max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
        max_scores.append(max_score)

        for i in range(len(promoter_sequence) - seq_length + 1):
            seq = promoter_sequence[i:i + seq_length]
            score = calculate_score(seq, matrix)
            normalized_score = score / max_score  # Score normalisé
            results.append({'Matrix': matrix_name, 'Sequence': seq, 'Score': normalized_score})

    # Affichage des résultats dans un tableau avec 3 chiffres significatifs
    table_data = []
    for result in results:
        table_data.append([result['Matrix'], result['Sequence'], f"{result['Score']:.3g}"])

    table = tabulate(table_data, headers=['Matrix', 'Sequence', 'Score'], tablefmt='orgtbl')
    output_text.delete(1.0, tk.END)
    output_text.insert(tk.END, table)

    max_score_label.config(text=f"Score maximums: {', '.join([f'{score:.3g}' for score in max_scores])}")

def calculate_score(sequence, matrix):
    score = 0
    for i, base in enumerate(sequence):
        if base in {'A', 'C', 'G', 'T'}:
            base_score = matrix[base]
            score += base_score[i]
    return score

def transform_matrix(matrix):
    reversed_matrix = {base: list(reversed(scores)) for base, scores in matrix.items()}
    complement_matrix = {
        'A': matrix['T'],
        'C': matrix['G'],
        'G': matrix['C'],
        'T': matrix['A']
    }
    reversed_complement_matrix = {base: list(reversed(scores)) for base, scores in complement_matrix.items()}

    return {
        'Original': matrix,
        'Reversed': reversed_matrix,
        'Complement': complement_matrix,
        'Reversed Complement': reversed_complement_matrix
    }

def matrix_extraction():
    # Récupération des valeurs des champs d'entrée
    jaspar_id = entry_jaspar.get()
    promoter_sequence = entry_promoter.get()

    # Appel à l'API de Jaspar pour récupérer la matrice de fréquence
    url = f"https://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
    response = requests.get(url)
    if response.status_code == 200:
        response_data = response.json()
        matrix = response_data['pfm']
    else:
        messagebox.showerror("Erreur", f"Erreur lors de la récupération de la matrice de fréquence : {response.status_code}")
        return

    # Transformation de la matrice de fréquence
    matrices = transform_matrix(matrix)

    # Recherche de la séquence référence dans chaque matrice
    search_sequence(jaspar_id, promoter_sequence, matrices)

# Création de la fenêtre principale
window = tk.Tk()
window.title("Recherche de séquence Jaspar")
window.geometry("400x350")

# Frame pour les champs d'entrée et le bouton
input_frame = ttk.Frame(window)
input_frame.pack(pady=20)

# Champ d'entrée pour l'ID Jaspar
label_jaspar = ttk.Label(input_frame, text="ID Jaspar:")
label_jaspar.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
entry_jaspar = ttk.Entry(input_frame, width=30)
entry_jaspar.grid(row=0, column=1, padx=5, pady=5)

# Champ d'entrée pour le promoteur
label_promoter = ttk.Label(input_frame, text="Promoteur:")
label_promoter.grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
entry_promoter = ttk.Entry(input_frame, width=30)
entry_promoter.grid(row=1, column=1, padx=5, pady=5)

# Bouton de recherche
search_button = ttk.Button(window, text="Rechercher", command=matrix_extraction)
search_button.pack()

# Étiquette pour afficher les scores maximums
max_score_label = ttk.Label(window, text="Score maximums: ")
max_score_label.pack()

# Zone de texte pour afficher les résultats
output_text = tk.Text(window, width=40, height=10)
output_text.pack(pady=20)

# Lancement de la boucle principale
window.mainloop()
