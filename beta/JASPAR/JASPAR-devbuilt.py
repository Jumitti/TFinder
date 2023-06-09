import requests
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tabulate import tabulate

def search_sequence():
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

    # Recherche de la séquence référence dans le promoteur
    results = []
    seq_length = len(matrix['A'])  # Utiliser la longueur de n'importe quelle clé, car elles doivent toutes avoir la même longueur
    for i in range(len(promoter_sequence) - seq_length + 1):
        seq = promoter_sequence[i:i+seq_length]
        score = calculate_score(seq, matrix)
        results.append({'Sequence': seq, 'Score': score})

    # Affichage des résultats dans un tableau
    table_data = []
    for result in results:
        table_data.append([result['Sequence'], result['Score']])

    table = tabulate(table_data, headers=['Sequence', 'Score'], tablefmt='orgtbl')
    output_text.delete(1.0, tk.END)
    output_text.insert(tk.END, table)

def calculate_score(sequence, matrix):
    score = 0
    for i, base in enumerate(sequence):
        if base in {'A', 'C', 'G', 'T'}:
            base_scores = matrix[base]
            score += base_scores[i]
    return score


# Création de la fenêtre principale
window = tk.Tk()
window.title("Recherche de séquence Jaspar")
window.geometry("400x300")

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
search_button = ttk.Button(window, text="Rechercher", command=search_sequence)
search_button.pack()

# Zone de texte pour afficher les résultats
output_text = tk.Text(window, width=40, height=10)
output_text.pack(pady=20)

# Lancement de la boucle principale
window.mainloop()
