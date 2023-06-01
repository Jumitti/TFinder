import os
import pandas as pd
import pyperclip
import requests
import tkinter as tk
import tkinter.messagebox as messagebox
import webbrowser

from tkinter import ttk
from PIL import Image, ImageTk
from tabulate import tabulate
from tkinter import filedialog

#Copy/Paste Button
def copy_sequence():
    sequence = sequence_text.get('1.0', tk.END).strip()
    pyperclip.copy(sequence)
    messagebox.showinfo("Copy", "The sequence has been copied to the clipboard.")

def paste_sequence():
    sequence = window.clipboard_get()
    text_promoter.delete("1.0", "end")
    text_promoter.insert("1.0", sequence)
    messagebox.showinfo("Paste", "The sequence has been pasted.")

def load_fasta():
    file_path = filedialog.askopenfilename(filetypes=[("FASTA Files", "*.fasta")])
    if file_path:
        with open(file_path, 'r') as file:
            fasta_content = file.read()
            sequence = parse_fasta(fasta_content)
            text_promoter.delete("1.0", "end")
            text_promoter.insert("1.0", sequence)
            messagebox.showinfo("Load FASTA", "FASTA sequence loaded successfully.")

# Reverse complement
def reverse_complement(sequence):
    text_statut.delete("1.0", "end")
    statut = "Reverse complement..."
    text_statut.insert("1.0", statut)
    window.update_idletasks()
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
    text_statut.delete("1.0", "end")
    statut = "Reverse complement -> Done"
    text_statut.insert("1.0", statut)
    window.update_idletasks()
    return complement_sequence

# Generation of all responsive elements
def generate_variants(sequence):
    text_statut.delete("1.0", "end")
    statut = "Generate responsive elements variants..."
    text_statut.insert("1.0", statut)
    window.update_idletasks()    
    variants = []

    # Original sequence
    variants.append(sequence)

    # Reverse sequence
    variants.append(sequence[::-1])

    # Complementary sequence
    complement_sequence = "".join(reverse_complement(base) for base in sequence)
    variants.append(complement_sequence)
    complement_mirror_sequence = complement_sequence[::-1]
    variants.append(complement_mirror_sequence)
    
    text_statut.delete("1.0", "end")
    statut = "Generate responsive elemnts variants -> Done"
    text_statut.insert("1.0", statut)
    window.update_idletasks()   
    return variants

# IUPAC code
def generate_iupac_variants(sequence):
    text_statut.delete("1.0", "end")
    statut = "Generate IUPAC variants..."
    text_statut.insert("1.0", statut)
    window.update_idletasks()   
    iupac_codes = {
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "M": ["A", "C"],
        "K": ["G", "T"],
        "W": ["A", "T"],
        "S": ["C", "G"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"]
    }

    sequences = [sequence]
    for i, base in enumerate(sequence):
        if base.upper() in iupac_codes:
            new_sequences = []
            for seq in sequences:
                for alternative in iupac_codes[base.upper()]:
                    new_sequence = seq[:i] + alternative + seq[i + 1:]
                    new_sequences.append(new_sequence)
            sequences = new_sequences
    
    text_statut.delete("1.0", "end")
    statut = "Generate IUPAC variants -> Done"
    text_statut.insert("1.0", statut)
    window.update_idletasks()
    return sequences

# Responsive Elements Finder (consensus sequence)
def find_sequence_consensus():
    global table
    table = []
    text_result.delete("1.0", "end")
    sequence_consensus_input = entry_sequence.get()
    text_statut.delete("1.0", "end")
    statut = "Find responsive elements..."
    text_statut.insert("1.0", statut)
    window.update_idletasks()   
    tis_value = int(entry_tis.get())

    # Transform with IUPAC code
    sequence_consensus = generate_iupac_variants(sequence_consensus_input)

    threshold = float(threshold_entry.get())

    found_positions = []

    # Responsive elements finder
    lines = text_promoter.get("1.0", "end-1c").split("\n")
    promoters = []

    first_line = lines[0]

    if first_line.startswith(("A", "T", "C", "G")):
        shortened_promoter_name = "n.d."
        promoter_region = first_line
        promoters.append((shortened_promoter_name, promoter_region))
    else :
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.startswith(">"):
                promoter_name = line[1:]
                shortened_promoter_name = promoter_name[:10] if len(promoter_name) > 10 else promoter_name
                promoter_region = lines[i+1]
                promoters.append((shortened_promoter_name, promoter_region))
                i += 2  # Passer à la ligne suivante
            else:
                i += 1

    # Affichage des noms et des séquences correspondantes
    for shortened_promoter_name, promoter_region in promoters:
        print("Nom:", shortened_promoter_name)
        print("Séquence:", promoter_region)
        for consensus in sequence_consensus:
            variants = generate_variants(consensus)
            max_mismatches = len(variants[0]) // 4  # Mismatches authorized
            for variant in variants:
                variant_length = len(variant)

                for i in range(len(promoter_region) - variant_length + 1):
                    sequence = promoter_region[i:i + variant_length]

                    mismatches = sum(a != b for a, b in zip(sequence, variant))  # Mismatches with Hamming distance

                    if mismatches <= max_mismatches:
                        # Eliminates short responsive elements that merge with long ones
                        similar_position = False
                        for position, _, _, _, _ in found_positions:
                            if abs(i - position) <= 1:
                                similar_position = True
                                break

                        if not similar_position:
                            homology_percentage = (variant_length - mismatches) / variant_length * 100  # % Homology

                            found_positions.append((i, sequence, variant, mismatches, homology_percentage))

        # Sort positions in ascending order
        found_positions.sort(key=lambda x: x[0])

        # Creating a results table
        if len(found_positions) > 0:
            for position, sequence, variant, mismatches, homology_percentage in found_positions:
                start_position = max(0, position - 3)
                end_position = min(len(promoter_region), position + len(sequence) + 3)
                sequence_with_context = promoter_region[start_position:end_position]

                sequence_parts = []
                for j in range(start_position, end_position):
                    if j < position or j >= position + len(sequence):
                        sequence_parts.append(sequence_with_context[j - start_position].lower())
                    else:
                        sequence_parts.append(sequence_with_context[j - start_position].upper())

                sequence_with_context = ''.join(sequence_parts)
                tis_position = position - tis_value

                row = [position, tis_position, sequence_with_context, homology_percentage, variant, shortened_promoter_name]
                table.append(row)

            table.sort(key=lambda x: (x[5], float(x[3])), reverse=False)

            # Filter results based on threshold
            filtered_table = [row for row in table if float(row[3]) >= threshold]
            
            filtered_table = sorted(filtered_table, key=lambda x: (x[5], -float(x[3])))

            if len(filtered_table) > 0:
                result = tabulate(filtered_table, headers=header, tablefmt="pipe")
            else:
                result = "No consensus sequence found with the specified threshold."
        else:
            result = "No consensus sequence found in the promoter region."
            
        text_statut.delete("1.0", "end")
        statut = "Find sequence -> Done"
        text_statut.insert("1.0", statut)
        window.update_idletasks()
        
        text_result.delete("1.0", "end")
        text_result.insert("1.0", result)

#Def table
table = []
header = ["Position", "Position (TIS)", "Sequence", "% Homology", "Ref seq", "Prom."]

#Export to excel
def export_to_excel():
    file_path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel Files", "*.xlsx")])
    
    if file_path:
        try:
            df = pd.DataFrame(table, columns=header)
            df.to_excel(file_path, index=False)
            messagebox.showinfo("Export Successful", f"Table exported to {file_path}")
        except Exception as e:
            messagebox.showerror("Export Failed", f"An error occurred while exporting the table:\n{str(e)}")
    else:
        messagebox.showwarning("Export Cancelled", "Export operation was cancelled by the user.")

#HELP
def show_help_PDF():
    script_directory = os.path.dirname(os.path.abspath(__file__))
    pdf_path = os.path.join(script_directory, "Promoter_finder_HELP.pdf")
    webbrowser.open(pdf_path)

#Logo
script_dir = os.path.dirname(os.path.abspath(__file__))
image_path = os.path.join(script_dir, "REF.png")

# Github
def open_site():
    url = "https://github.com/Jumitti/Responsive-Elements-Finder"
    webbrowser.open(url)
    
# Create TK windows
window = tk.Tk()
window.title("Responsive Elements Finder")
logo_image = tk.PhotoImage(file=image_path)
window.iconphoto(True, logo_image)

#How to use
help_button = tk.Button(window, text="How to use", command=show_help_PDF)
help_button.place(x=10, y=10)

# Github
button = tk.Button(window, text="Github", command=open_site)
button.place(x=10, y=40)

# Credit
credit_label = tk.Label(window, text="By MINNITI Julien")
credit_label.place(x=10, y=70)

# Section "Responsive Elements finder"
section_responsive_finder = tk.LabelFrame(window, text="Responsive Elements Finder")
section_responsive_finder.grid(row=0, column=2, padx=10, pady=10)

# Promoter entry
label_promoter = tk.Label(section_responsive_finder, text="Promoter")
label_promoter.grid(row=0, column=0)
text_promoter = tk.Text(section_responsive_finder, height=2, width=50)
text_promoter.grid(row=1, column=0)
paste_button = tk.Button(section_responsive_finder, text="Past", command=paste_sequence)
paste_button.grid(row=2, column=0)

# RE entry
label_sequence = tk.Label(section_responsive_finder, text="Responsive element (IUPAC authorize) :")
label_sequence.grid(row=4, column=0)
entry_sequence = tk.Entry(section_responsive_finder)
entry_sequence.grid(row=5, column=0)

# TIS entry
label_tis = tk.Label(section_responsive_finder, text="Transcription Initiation Site (bp)")
label_tis.grid(row=6, column=0)
label_tis_info = tk.Label(section_responsive_finder, text="(distance from start of promoter or 'Upstream' if you use the Promoter Finder)")
label_tis_info.grid(row=7, column=0)
entry_tis = tk.Entry(section_responsive_finder)
entry_tis.insert(0, "0")  # $"0" default
entry_tis.grid(row=8, column=0)

# Threshold
threshold_label = tk.Label(section_responsive_finder, text="Threshold (%)")
threshold_label.grid(row=9, column=0)
threshold_entry = tk.Entry(section_responsive_finder)
threshold_entry.insert(80, "80")  # $"80" default
threshold_entry.grid(row=10, column=0)

# Find RE
button_search = tk.Button(section_responsive_finder, text="Find responsive elements", command=find_sequence_consensus)
button_search.grid(row=11, column=0)

# RE output
label_result = tk.Label(section_responsive_finder, text="Responsive elements")
label_result.grid(row=12, column=0)
text_result = tk.Text(section_responsive_finder, height=16, width=100)
text_result.grid(row=13, column=0)

# Création du bouton Export to Excel
export_button = tk.Button(section_responsive_finder, text="Export to Excel", command=export_to_excel)
export_button.grid(row=14, column=0)

# Section "Statut"
section_statut = tk.LabelFrame(window, text="Statut")
section_statut.grid(row=15, column=2, padx=10, pady=10)

# Statut output
text_statut = tk.Text(section_statut, height=1, width=100)
text_statut.grid(row=16, column=0)

# Configure grid weights
window.grid_rowconfigure(0, weight=1)
window.grid_columnconfigure(1, weight=1)
window.grid_columnconfigure(2, weight=1)
section_responsive_finder.grid_rowconfigure(12, weight=1)

# Lancement de la boucle principale
window.mainloop()