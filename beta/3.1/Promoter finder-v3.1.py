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

# Convert gene to ENTREZ_GENE_ID
def convert_gene_names_to_entrez_ids(gene_names):
    try:
        entrez_ids = []
        for gene_id in gene_names:
            species = species_combobox.get()
            # Requête de recherche de gène par nom et espèce
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_id}[Gene%20Name]+AND+{species}[Organism]&retmode=json&rettype=xml"
            response = requests.get(url)

            if response.status_code == 200:
                response_data = response.json()

                if response_data['esearchresult']['count'] == '0':
                    raise Exception(f"No gene found for name: {gene_name}")

                else:
                    gene_id = response_data['esearchresult']['idlist'][0]
                    entrez_ids.append(gene_id)

            else:
                raise Exception(f"Error during gene search: {response.status_code}")

        return entrez_ids

    except Exception as e:
        raise Exception(f"Error: {str(e)}")

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
        # Calculating extraction coordinates on the chromosome
        if chrstop > chrstart:
            start = chrstart - upstream
            end = chrstart + downstream
            
            # Request for DNA sequence
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&from={start}&to={end}&rettype=fasta&retmode=text"
            response = requests.get(url)
            
            if response.status_code == 200:
                # Extraction of DNA sequence
                dna_sequence = response.text.split('\n', 1)[1].replace('\n', '')
                
                sequence = dna_sequence
            
            else:
                raise Exception(f"An error occurred while retrieving the DNA sequence : {response.status_code}")
        
        else:
            start = chrstart + upstream
            end = chrstart - downstream
            
            # Request for DNA sequence
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&from={end}&to={start}&rettype=fasta&retmode=text"
            response = requests.get(url)
            
            if response.status_code == 200:
                # Extraction of DNA sequence
                dna_sequence = response.text.split('\n', 1)[1].replace('\n', '')
                
                
                sequence = dna_sequence
                sequence = reverse_complement(sequence)
                
            else:
                raise Exception(f"An error occurred while retrieving the DNA sequence : {response.status_code}")

        return sequence

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
    species = species_combobox.get()
    gene_ids = gene_id_entry.get("1.0", tk.END).strip().split("\n")
    total_gene_ids = len(gene_ids)
    upstream = int(upstream_entry.get())
    downstream = int(downstream_entry.get())
    result_text.delete("1.0", tk.END)
    
    for i, gene_id in enumerate(gene_ids, start=1):
        try:
            number_gene_id = i
            
            # Convertir le nom du gène en ENTREZ_GENE_ID si nécessaire
            gene_names = []  # Déclaration de gene_names
            if not gene_id.isdigit():
                gene_id = gene_id.strip("'\'")
                gene_names.append(gene_id)
                gene_id = convert_gene_names_to_entrez_ids(gene_names)
                gene_entrez_id = gene_id.copy()
                
            else:
                gene_id = gene_id.strip("'\"")
                gene_entrez_id = [gene_id]
            
        except Exception as e:
            result_text.insert(tk.END, f"Error retrieving gene information for ID: {gene_id}\nError: {str(e)}\n")
                
        # Gene information retrieval
        for gene_id in gene_entrez_id:
            text_statut.delete("1.0", "end")
            text_statut.insert("1.0", f"Find gene information... ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()
            gene_info = get_gene_info(gene_id, species)
            gene_name = gene_info['name']
            text_statut.delete("1.0", "end")
            text_statut.insert("1.0", f"Find {gene_name} information -> Done ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()

            chraccver = gene_info['genomicinfo'][0]['chraccver']
            chrstart = gene_info['genomicinfo'][0]['chrstart']
            chrstop = gene_info['genomicinfo'][0]['chrstop']

            # Promoter retrieval
            text_statut.delete("1.0", "end")
            text_statut.insert("1.0", f"Extract {gene_name} promoter... ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()
            dna_sequence = get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream)
            text_statut.delete("1.0", "end")
            text_statut.insert("1.0", f"Extract {gene_name} promoter -> Done ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()

            # Append the result to the result_text
            result_text.insert(tk.END, f">{gene_name} | {species} | {chraccver} | TSS: {chrstart}\n{dna_sequence}\n\n")
            text_statut.delete("1.0", "end")
            text_statut.insert("1.0", f"Extract promoter -> Done ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()

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
species_combobox.current(0)
species_combobox.grid(row=3, column=0)

# Upstream distance entry
upstream_label = tk.Label(section_promoter_finder, text="Upstream Distance:")
upstream_label.grid(row=0, column=1)
upstream_entry = tk.Entry(section_promoter_finder, width=10)
upstream_entry.grid(row=1, column=1)
upstream_entry.insert(tk.END, "2000")

# Downstream distance entry
downstream_label = tk.Label(section_promoter_finder, text="Downstream Distance:")
downstream_label.grid(row=2, column=1)
downstream_entry = tk.Entry(section_promoter_finder, width=10)
downstream_entry.grid(row=3, column=1)
downstream_entry.insert(tk.END, "200")

# Promoter Finder button
find_button = tk.Button(section_promoter_finder, text="Find Promoters", command=get_sequence)
find_button.grid(row=4, column=0, columnspan=2, pady=10)

# Copy/Paste section
copy_paste_frame = tk.Frame(window)
copy_paste_frame.grid(row=1, column=1, padx=10)

# Copy Button
copy_button = tk.Button(copy_paste_frame, text="Copy", command=copy_sequence)
copy_button.grid(row=0, column=0, padx=5)

# Paste Button
paste_button = tk.Button(copy_paste_frame, text="Paste", command=paste_sequence)
paste_button.grid(row=0, column=1, padx=5)

# Result section
result_frame = tk.LabelFrame(window, text="Promoter Sequences")
result_frame.grid(row=0, column=2, rowspan=2, padx=10, pady=10)

# Result text
result_text = tk.Text(result_frame, height=20, width=50)
result_text.pack(side=tk.LEFT, fill=tk.BOTH)
result_scroll = tk.Scrollbar(result_frame, orient=tk.VERTICAL)
result_scroll.config(command=result_text.yview)
result_scroll.pack(side=tk.RIGHT, fill=tk.Y)
result_text.config(yscrollcommand=result_scroll.set)

# Status section
status_frame = tk.LabelFrame(window, text="Status")
status_frame.grid(row=2, column=1, columnspan=2, padx=10, pady=10)

# Status text
text_statut = tk.Text(status_frame, height=2, width=50)
text_statut.pack(side=tk.LEFT, fill=tk.BOTH)
scrollbar_statut = tk.Scrollbar(status_frame, orient=tk.VERTICAL)
scrollbar_statut.config(command=text_statut.yview)
scrollbar_statut.pack(side=tk.RIGHT, fill=tk.Y)
text_statut.config(yscrollcommand=scrollbar_statut.set)

window.mainloop()
