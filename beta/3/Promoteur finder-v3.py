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

#Gene informations
def get_gene_info(gene_id, species):
    try:
        # Request for gene informations
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json&rettype=xml&species={species}"
        response = requests.get(url)

        if response.status_code == 200:
            response_data = response.json()

            # Extraction of gene informations
            gene_info = response_data['result'][str(gene_id)]

            return gene_info

        else:
            raise Exception(f"Error during extraction of gene information : {response.status_code}")

    except Exception as e:
        raise Exception(f"Error : {str(e)}")

#Promoter finder
def get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream):
    try:
        # Request for DNA sequence
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&rettype=fasta&retmode=text"
        response = requests.get(url)

        if response.status_code == 200:
            # Extraction of DNA sequence
            dna_sequence = response.text.split('\n', 1)[1].replace('\n', '')

            # Calculating extraction coordinates on the chromosome
            if chrstop > chrstart:
                start = chrstart - upstream
                end = chrstart + downstream
                sequence = dna_sequence[start:end]
            else:
                start = chrstart + upstream
                end = chrstart - downstream
                sequence = dna_sequence[end:start]
                sequence = reverse_complement(sequence)
            return sequence

        else:
            raise Exception(f"An error occurred while retrieving the DNA sequence : {response.status_code}")

    except Exception as e:
        raise Exception(f"Error : {str(e)}")

# Copy/Paste Button
def copy_sequence():
    sequence = result_text.get('1.0', tk.END).strip()
    pyperclip.copy(sequence)
    messagebox.showinfo("Copy", "The sequence has been copied to the clipboard.")

def paste_sequence():
    sequence = window.clipboard_get()
    text_promoter.delete("1.0", "end")
    text_promoter.insert("1.0", sequence)
    messagebox.showinfo("Paste", "The sequence has been pasted.")

# Display gene and promoter
def get_sequence():
    gene_ids = gene_id_entry.get("1.0", tk.END).strip().split("\n")
    total_gene_ids = len(gene_ids)
    species = species_combobox.get()
    upstream = int(upstream_entry.get())
    downstream = int(downstream_entry.get())
    result_text.delete("1.0", tk.END)
    for i, gene_id in enumerate(gene_ids, start=1):
        try:
            number_gene_id = i
            
            # Gene information retrieval
            text_statut.delete("1.0", "end")
            text_statut.insert("1.0", f"Find gene information... ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()
            gene_info = get_gene_info(gene_id, species)
            gene_name = gene_info['name']
            text_statut.delete("1.0", "end")
            text_statut.insert("1.0", f"Find {gene_name} information -> Done ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()
            
            text_statut.delete("1.0", "end")
            text_statut.insert("1.0", f"Extract {gene_name} promoter... ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()
            chraccver = gene_info['genomicinfo'][0]['chraccver']
            chrstart = gene_info['genomicinfo'][0]['chrstart']
            chrstop = gene_info['genomicinfo'][0]['chrstop']

            # Promoter retrieval
            dna_sequence = get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream)
            text_statut.delete("1.0", "end")
            text_statut.insert("1.0", f"Extract {gene_name} promoter -> Done ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()
            
            # Append the result to the result_text
            result_text.insert(tk.END, f">{gene_name} | {chraccver} | TIS: {chrstart}\n{dna_sequence}\n\n")
            text_statut.delete("1.0", "end")
            text_statut.insert("1.0", f"Extract promoter -> Done ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()
        except Exception as e:
            result_text.insert(tk.END, f"Error retrieving gene information for ID: {gene_id}\nError: {str(e)}\n")

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
window.title("Promoter finder")
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

# Section "Promoter finder"
section_promoter_finder = tk.LabelFrame(window, text="Promoter Finder")
section_promoter_finder.grid(row=0, column=1, padx=10, pady=10)

# Gene ID entry
gene_id_label = tk.Label(section_promoter_finder, text="Gene ID:")
gene_id_label.grid(row=0, column=0)
gene_id_entry = tk.Text(section_promoter_finder, height=5, width=20)
gene_id_entry.grid(row=1, column=0)


# Species selection
species_label = tk.Label(section_promoter_finder, text="Species:")
species_label.grid(row=2, column=0)
species_combobox = ttk.Combobox(section_promoter_finder, values=["Human", "Mouse", "Rat"])
species_combobox.grid(row=3, column=0)

# Upstream/downstream entry
upstream_label = tk.Label(section_promoter_finder, text="Upstream (bp):")
upstream_label.grid(row=4, column=0)
upstream_entry = tk.Entry(section_promoter_finder)
upstream_entry.insert(2000, "2000")  # $"2000" default
upstream_entry.grid(row=5, column=0)

downstream_label = tk.Label(section_promoter_finder, text="Downstream (bp):")
downstream_label.grid(row=6, column=0)
downstream_entry = tk.Entry(section_promoter_finder)
downstream_entry.insert(500, "500")  # $"500" default
downstream_entry.grid(row=7, column=0)

# Search
search_button = tk.Button(section_promoter_finder, text="Find promoter  (CAN BE STUCK ! Don't worry, just wait)", command=get_sequence)
search_button.grid(row=8, column=0)

# Section "Results"
section_results = tk.LabelFrame(window, text="Results")
section_results.grid(row=15, column=1, padx=10, pady=10)

# Results output
result_text = tk.Text(section_results, height=10, width=50)
result_text.grid(row=0, column=0)
copy_button = tk.Button(section_results, text="Copy", command=copy_sequence)
copy_button.grid(row=1, column=0)

# Section "Statut"
section_statut = tk.LabelFrame(window, text="Statut")
section_statut.grid(row=16, column=1, padx=10, pady=10)

# Statut output
text_statut = tk.Text(section_statut, height=1, width=50)
text_statut.grid(row=0, column=0)

# Configure grid weights
window.grid_rowconfigure(0, weight=1)
window.grid_columnconfigure(1, weight=1)
window.grid_columnconfigure(2, weight=1)
section_promoter_finder.grid_rowconfigure(13, weight=1)

# Lancement de la boucle principale
window.mainloop()