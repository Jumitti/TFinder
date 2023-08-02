# Copyright (c) 2023 Minniti Julien

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of TFinder and associated documentation files, to deal
# in TFinder without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of TFinder, and to permit persons to whom TFinder is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of TFinder.

# TFINDER IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH TFINDER OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import requests
import tkinter as tk
from tkinter import ttk
import tkinter.messagebox as messagebox
import pyperclip
from PIL import Image, ImageTk
import webbrowser
import os
from tabulate import tabulate

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

#Reverse complement
def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
    return complement_sequence

#Copy/Paste Button
def copy_sequence():
    sequence = sequence_text.get('1.0', tk.END).strip()
    pyperclip.copy(sequence)
    messagebox.showinfo("Copy", "The sequence has been copied to the clipboard.")

def paste_sequence():
    sequence = window.clipboard_get()
    text_promoter.delete("1.0", "end")
    text_promoter.insert("1.0", sequence)

#Display gene and promoter
def get_sequence():
    gene_id = gene_id_entry.get()
    species = species_combobox.get()
    upstream = int(upstream_entry.get())
    downstream = int(downstream_entry.get())

    try:
        # Gene information retrieval
        gene_info = get_gene_info(gene_id, species)
        gene_name = gene_info['name']
        chraccver = gene_info['genomicinfo'][0]['chraccver']
        chrstart = gene_info['genomicinfo'][0]['chrstart']
        chrstop = gene_info['genomicinfo'][0]['chrstop']

        # Display informations
        gene_name_label.config(text=f"Gene ID : {gene_name}")
        gene_chr_label.config(text=f"Chromosome : {chraccver}")
        chrstart_label.config(text=f"Transcription Initiation Site : {chrstart}")

        # Promoter retrieval
        dna_sequence = get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream)

        # Display promoter
        sequence_text.delete('1.0', tk.END)
        sequence_text.insert(tk.END, dna_sequence)

    except Exception as e:
        gene_name_label.config(text="Gene ID :")
        gene_chr_label.config(text="Chromosome :")
        chrstart_label.config(text=f"Transcription Initiation Site : {chrstart}")
        sequence_text.delete('1.0', tk.END)
        messagebox.showerror("Error", f"Error : {str(e)}")

#Generation of all responsive elements
def generate_variants(sequence):
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

    # Sequence without first or last nucleotide (and reverse/complement)
    if len(sequence) > 1:
        variants.append(sequence[1:])
        variants.append(sequence[::-1][1:])
        variants.append(complement_sequence[1:])
        variants.append(complement_mirror_sequence[1:])

    if len(sequence) > 1:
        variants.append(sequence[:-1])
        variants.append(sequence[::-1][:-1])
        variants.append(complement_sequence[:-1])
        variants.append(complement_mirror_sequence[:-1])

    return variants

# IUPAC code
def generate_iupac_variants(sequence):
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

    return sequences

# Responsive Elements Finder (consensus sequence)
def find_sequence_consensus():
    promoter_region = text_promoter.get("1.0", "end-1c")
    sequence_consensus_input = entry_sequence.get()
    tis_value = int(entry_tis.get())

    # Transform with IUPAC code
    sequence_consensus = generate_iupac_variants(sequence_consensus_input)

    found_positions = []

    # Responsive elements finder
    for consensus in sequence_consensus:
        variants = generate_variants(consensus)

        for variant in variants:
            variant_length = len(variant)

            for i in range(len(promoter_region) - variant_length + 1):
                sequence = promoter_region[i:i + variant_length]

                if sequence.upper() == variant.upper():
                    # Eliminates short responsives elements that merge with long ones 
                    similar_position = False
                    for position, _ in found_positions:
                        if abs(i - position) <= 3:
                            similar_position = True
                            break

                    if not similar_position:
                        found_positions.append((i, sequence))

    # Sort positions in ascending order
    found_positions.sort(key=lambda x: x[0])

    # Creating a results table
    if len(found_positions) > 0:
        table = []
        header = ["Position", "Position (TIS)", "Sequence", "Lenght"]

        for position, sequence in found_positions:
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
            row = [position, tis_position, sequence_with_context, len(sequence)]
            table.append(row)

        result = tabulate(table, headers=header, tablefmt="pipe")

    else:
        result = "No consensus sequence found in the promoter region."

    text_result.delete("1.0", "end")
    text_result.insert("1.0", result)

#HELP
def show_help_image():
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
help_label = tk.Label(window, text="How to use", cursor="hand2")
help_label.place(x=10, y=10)
help_label.bind("<Button-1>", lambda event: show_help_image())

#Credit
credit_label = tk.Label(window, text="By MINNITI Julien")
credit_label.place(x=10, y=40)

#Github
button = tk.Button(window, text="Github", command=open_site)
button.place(x=10 , y=70)


# Section "Promoter finder"
section_promoter_finder = tk.LabelFrame(window, text="Promoter Finder")
section_promoter_finder.pack(padx=10, pady=10)

# Gene ID entry
gene_id_label = tk.Label(section_promoter_finder, text="Gene ID:")
gene_id_label.pack()
gene_id_entry = tk.Entry(section_promoter_finder)
gene_id_entry.pack()

# Species selection
species_label = tk.Label(section_promoter_finder, text="Species:")
species_label.pack()
species_combobox = ttk.Combobox(section_promoter_finder, values=["Human", "Mouse", "Rat"])
species_combobox.pack()

# Upstream/downstream entry
upstream_label = tk.Label(section_promoter_finder, text="Upstream (bp):")
upstream_label.pack()
upstream_entry = tk.Entry(section_promoter_finder)
upstream_entry.pack()

downstream_label = tk.Label(section_promoter_finder, text="Downstream (bp):")
downstream_label.pack()
downstream_entry = tk.Entry(section_promoter_finder)
downstream_entry.pack()

#Search
search_button = tk.Button(section_promoter_finder, text="Find promoter  (CAN BE STUCK ! Don't worry, just wait)", command=get_sequence)
search_button.pack()

# Gene informations
gene_name_label = tk.Label(section_promoter_finder, text="Gene name :")
gene_name_label.pack()

gene_chr_label = tk.Label(section_promoter_finder, text="Chromosome :")
gene_chr_label.pack()

chrstart_label = tk.Label(section_promoter_finder, text="Transcription Initiation Site (on the chromosome) :")
chrstart_label.pack()

# Promoter sequence
promoter_label = tk.Label(section_promoter_finder, text="Promoter")
promoter_label.pack()
sequence_text = tk.Text(section_promoter_finder, height=3, width=50)
sequence_text.pack()
copy_button = tk.Button(section_promoter_finder, text="Copy", command=copy_sequence)
copy_button.pack()

# Section "Responsive Elements finder"
section_responsive_finder = tk.LabelFrame(window, text="Responsive Elements Finder")
section_responsive_finder.pack(padx=10, pady=10)

# Promoter entry
label_promoter = tk.Label(section_responsive_finder, text="Promoter")
label_promoter.pack()
text_promoter = tk.Text(section_responsive_finder, height=3, width=50)
text_promoter.pack()
paste_button = tk.Button(section_responsive_finder, text="Past", command=paste_sequence)
paste_button.pack()

# RE entry
label_sequence = tk.Label(section_responsive_finder, text="Responsive element (IUPAC authorize) :")
label_sequence.pack()
entry_sequence = tk.Entry(section_responsive_finder)
entry_sequence.pack()

#TIS entry
label_tis = tk.Label(section_responsive_finder, text="Transcription Initiation Site (bp)")
label_tis.pack()
label_tis = tk.Label(section_responsive_finder, text="(distance from start of promoter or 'Upstream' if you use the Promoter Finder above)")
label_tis.pack()
entry_tis = tk.Entry(section_responsive_finder)
entry_tis.insert(0, "0")  # $"0" default
entry_tis.pack()

#Find RE
button_search = tk.Button(section_responsive_finder, text="Find responsive elements", command=find_sequence_consensus)
button_search.pack()

#RE output
label_result = tk.Label(section_responsive_finder, text="Responsive elements")
label_result.pack()
text_result = tk.Text(section_responsive_finder, height=10, width=100)
text_result.pack()

# Lancement de la boucle principale
window.mainloop()
