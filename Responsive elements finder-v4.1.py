import itertools
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

#Promoter Finder

# Convert gene to ENTREZ_GENE_ID
def convert_gene_names_to_entrez_ids(gene_names):
    try:
        entrez_ids = []
        for gene_id in gene_names:
            species = species_combobox.get()
            # Request for ENTREZ_GENE_ID
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

# Display gene and promoter
def get_sequence():
    species = species_combobox.get()
    gene_ids = gene_id_entry.get("1.0", tk.END).strip().split("\n")
    total_gene_ids = len(gene_ids)
    upstream = int(upstream_entry.get())
    downstream = int(downstream_entry.get())
    result_promoter.delete("1.0", tk.END)
    
    for i, gene_id in enumerate(gene_ids, start=1):
        try:
            number_gene_id = i
            
            # gene name to ENTREZ_GENE_ID
            gene_names = []
            if not gene_id.isdigit():
                gene_id = gene_id.strip("'\'")
                gene_names.append(gene_id)
                gene_id = convert_gene_names_to_entrez_ids(gene_names)
                gene_entrez_id = gene_id.copy()
                
            else:
                gene_id = gene_id.strip("'\"")
                gene_entrez_id = [gene_id]
            
        except Exception as e:
            result_promoter.insert(tk.END, f"Error retrieving gene information for ID: {gene_id}\nError: {str(e)}\n")
                
        # Gene information retrieval
        for gene_id in gene_entrez_id:
            text_status.delete("1.0", "end")
            text_status.insert("1.0", f"Find gene information... ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()
            gene_info = get_gene_info(gene_id, species)
            gene_name = gene_info['name']
            text_status.delete("1.0", "end")
            text_status.insert("1.0", f"Find {gene_name} information -> Done ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()

            chraccver = gene_info['genomicinfo'][0]['chraccver']
            chrstart = gene_info['genomicinfo'][0]['chrstart']
            chrstop = gene_info['genomicinfo'][0]['chrstop']

            # Promoter retrieval
            text_status.delete("1.0", "end")
            text_status.insert("1.0", f"Extract {gene_name} promoter... ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()
            dna_sequence = get_dna_sequence(chraccver, chrstart, chrstop, upstream, downstream)
            text_status.delete("1.0", "end")
            text_status.insert("1.0", f"Extract {gene_name} promoter -> Done ({number_gene_id}/{total_gene_ids})")
            window.update_idletasks()

            # Append the result to the result_promoter
            result_promoter.insert(tk.END, f">{gene_name} | {species} | {chraccver} | TSS: {chrstart}\n{dna_sequence}\n\n")
            window.update_idletasks()
    
    text_status.delete("1.0", "end")
    text_status.insert("1.0", f"Extract promoter -> Done ({number_gene_id}/{total_gene_ids})")
    global entry_tis_var
    entry_tis_var.set(upstream)
    messagebox.showinfo("Promoter", "Promoters region extracted.")

# Responsive-Elements-Finder

# Responsive Elements Finder (consensus sequence)
def find_sequence_consensus():
    if use_jaspar.get() == 1 :
        print('ok')
        matrix_extraction()
        
    else:
        global table
        table = []
        text_result.delete("1.0", "end")
        sequence_consensus_input = entry_sequence.get()
        text_status.delete("1.0", "end")
        status = "Find responsive elements..."
        text_status.insert("1.0", status)
        window.update_idletasks()   
        tis_value = int(entry_tis.get())

        # Transform with IUPAC code
        sequence_consensus = generate_iupac_variants(sequence_consensus_input)

        threshold = float(threshold_entry.get())

        # Promoter input type
        lines = result_promoter.get("1.0", "end-1c")
        promoters = []
        
        first_line = lines
        if first_line.startswith(("A", "T", "C", "G")):
            shortened_promoter_name = "n.d."
            promoter_region = lines
            promoters.append((shortened_promoter_name, promoter_region))
        else :
            lines = result_promoter.get("1.0", "end-1c").split("\n")
            i = 0
            while i < len(lines):
                line = lines[i]
                if line.startswith(">"):
                    promoter_name = line[1:]
                    shortened_promoter_name = promoter_name[:10] if len(promoter_name) > 10 else promoter_name
                    promoter_region = lines[i+1]
                    promoters.append((shortened_promoter_name, promoter_region))
                    i += 2
                else:
                    i += 1
                    
        # REF
        for j, (shortened_promoter_name, promoter_region) in enumerate(promoters, start=1):
            
            found_positions = []
            
            total_promoter = len(promoters)
            
            pattern = "\|/-\|/-"
            cycle = itertools.cycle(pattern)
            
            for consensus in sequence_consensus:
                variants = generate_variants(consensus)
                for variant in variants:
                    
                    variant_length = len(variant)
                    
                    status_char = next(cycle)                 
                    text_status.delete("1.0", "end")
                    text_status.insert("1.0", f"Find responsive element in {shortened_promoter_name}...({j}/{total_promoter}) {status_char}")
                    window.update_idletasks()        

                    for i in range(len(promoter_region) - variant_length + 1):
                        sequence = promoter_region[i:i + variant_length]

                        mismatches = sum(a != b for a, b in zip(sequence, variant))  # Mismatches
                        
                        homology_percentage = (variant_length - mismatches) / variant_length * 100  # % Homology
                        
                        # Find best homology sequence
                        better_homology = False
                        for position, _, _, _, best_homology_percentage in found_positions:
                            if abs(i - position) < 1 and homology_percentage <= best_homology_percentage:
                                better_homology = True
                                break

                        if not better_homology:
                           
                            best_homology_percentage = (variant_length - mismatches) / variant_length * 100  # % Homology
                            
                            found_positions.append((i, sequence, variant, mismatches, best_homology_percentage))                            

            # Sort positions in ascending order
            found_positions.sort(key=lambda x: x[0])

            # Creating a results table
            if len(found_positions) > 0:
                for position, sequence, variant, mismatches, best_homology_percentage in found_positions:
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

                    row = [position, tis_position, sequence_with_context, best_homology_percentage, variant, shortened_promoter_name]
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
                
            text_status.delete("1.0", "end")
            text_status.insert("1.0", f"Find sequence -> Done ({total_promoter}/{total_promoter})")
            window.update_idletasks()
            
            text_result.delete("1.0", "end")
            text_result.insert("1.0", result)

# IUPAC code
def generate_iupac_variants(sequence):
    text_status.delete("1.0", "end")
    status = "Generate IUPAC variants..."
    text_status.insert("1.0", status)
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
    
    text_status.delete("1.0", "end")
    status = "Generate IUPAC variants -> Done"
    text_status.insert("1.0", status)
    window.update_idletasks()
    return sequences

# Generation of all responsive elements
def generate_variants(sequence):
    text_status.delete("1.0", "end")
    status = "Generate responsive elements variants..."
    text_status.insert("1.0", status)
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
    
    text_status.delete("1.0", "end")
    status = "Generate responsive elements variants -> Done"
    text_status.insert("1.0", status)
    window.update_idletasks()   
    return variants

#Find with JASPAR
def search_sequence(jaspar_id, matrices):
    text_result.delete("1.0", "end")
    results = []
    max_scores = []
    threshold = float(threshold_entry.get())
    tis_value = int(entry_tis.get())

    for matrix_name, matrix in matrices.items():
        seq_length = len(matrix['A'])

        # Max score per matrix
        max_score = sum(max(matrix[base][i] for base in matrix.keys()) for i in range(seq_length))
        max_scores.append(max_score)
        
        # Promoter input type
        lines = result_promoter.get("1.0", "end-1c")
        promoters = []
        
        first_line = lines
        if first_line.startswith(("A", "T", "C", "G")):
            shortened_promoter_name = "n.d."
            promoter_region = lines
            promoters.append((shortened_promoter_name, promoter_region))
        else :
            lines = result_promoter.get("1.0", "end-1c").split("\n")
            i = 0
            while i < len(lines):
                line = lines[i]
                if line.startswith(">"):
                    promoter_name = line[1:]
                    shortened_promoter_name = promoter_name[:10] if len(promoter_name) > 10 else promoter_name
                    promoter_region = lines[i+1]
                    promoters.append((shortened_promoter_name, promoter_region))
                    i += 2
                else:
                    i += 1
                    
        # REF
        for j, (shortened_promoter_name, promoter_region) in enumerate(promoters, start=1):
            
            found_positions = []
            
            total_promoter = len(promoters)
            
            pattern = "\|/-\|/-"
            cycle = itertools.cycle(pattern)
                    
            status_char = next(cycle)                 
            text_status.delete("1.0", "end")
            text_status.insert("1.0", f"Find responsive element in {shortened_promoter_name}...({j}/{total_promoter}) {status_char}")
            window.update_idletasks()        

            for i in range(len(promoter_region) - seq_length + 1):
                seq = promoter_region[i:i + seq_length]
                score = calculate_score(seq, matrix)
                normalized_score = (score / max_score)*100
                
                if normalized_score >= threshold:  
                    tis_position = i - tis_value
                    
                    sequence_with_case = promoter_region[i-3:i].lower() + seq.upper() + promoter_region[i+seq_length:i+seq_length+3].lower()
                    
                    results.append({
                        'Matrix': matrix_name,
                        'Sequence': sequence_with_case,
                        'Score': normalized_score,
                        'Position': i,
                        'TIS Position': tis_position,
                        'Promoter Name': shortened_promoter_name
                    })

    
    sorted_results = sorted(results, key=lambda x: (x['Score'], x['Promoter Name']), reverse=True)

    
    table_data = []
    for result in sorted_results:
        table_data.append([
            result['Position'],
            result['TIS Position'],
            result['Sequence'],
            f"{result['Score']:.3g}",
            result['Promoter Name']
        ])
    
    text_status.delete("1.0", "end")
    text_status.insert("1.0", f"Find sequence -> Done ({total_promoter}/{total_promoter})")
    window.update_idletasks()
    
    table = tabulate(table_data, headers=['Position', 'TIS Position', 'Sequence', 'Score %', 'Promoter Name'], tablefmt='orgtbl')
    text_result.delete("1.0", "end")
    text_result.insert("1.0", table)

#Extract JASPAR matrix
def matrix_extraction():
    jaspar_id = entry_sequence.get()
    url = f"https://jaspar.genereg.net/api/v1/matrix/{jaspar_id}/"
    response = requests.get(url)
    if response.status_code == 200:
        response_data = response.json()
        matrix = response_data['pfm']
    else:
        messagebox.showerror("Erreur", f"Erreur lors de la récupération de la matrice de fréquence : {response.status_code}")
        return

    # Transform matrix in revers, complement, reverse-complement
    matrices = transform_matrix(matrix)

    # search sequence
    search_sequence(jaspar_id,matrices)

#Transform JASPAR matrix
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

#Calculate score with JASPAR
def calculate_score(sequence, matrix):
    score = 0
    for i, base in enumerate(sequence):
        if base in {'A', 'C', 'G', 'T'}:
            base_score = matrix[base]
            score += base_score[i]
    return score


#Def table
table = []
header = ["Position", "Position (TIS)", "Sequence", "% Homology", "Ref seq", "Prom."]

# Common functions

# Reverse complement
def reverse_complement(sequence):
    text_status.delete("1.0", "end")
    status = "Reverse complement..."
    text_status.insert("1.0", status)
    window.update_idletasks()
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
    text_status.delete("1.0", "end")
    status = "Reverse complement -> Done"
    text_status.insert("1.0", status)
    window.update_idletasks()
    return complement_sequence

# Copy/Paste Button
def copy_sequence():
    sequence = result_promoter.get('1.0', tk.END).strip()
    pyperclip.copy(sequence)
    messagebox.showinfo("Copy", "The sequence has been copied to the clipboard.")

def paste_sequence():
    sequence = window.clipboard_get()
    result_promoter.delete("1.0", "end")
    result_promoter.insert("1.0", sequence)
    
def paste_gene():
    gene = window.clipboard_get()
    gene_id_entry.delete("1.0", "end")
    gene_id_entry.insert("1.0", gene)

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

# Github
def open_site():
    url = "https://github.com/Jumitti/Responsive-Elements-Finder"
    webbrowser.open(url)

# Github
def open_site_webui():
    url = "https://responsive-elements-finder2.streamlit.app/"
    webbrowser.open(url)

# TK window

#Logo
script_dir = os.path.dirname(os.path.abspath(__file__))
image_path = os.path.join(script_dir, "REF.png")

# Create TK windows
window = tk.Tk()
window.title("Responsive Elements Finder")
logo_image = tk.PhotoImage(file=image_path)
window.iconphoto(True, logo_image)

# Section "Promoter finder"
section_promoter_finder = tk.LabelFrame(window, text="Promoter Finder")
section_promoter_finder.grid(row=0, column=0, padx=10, pady=10)

# Gene ID entry
gene_id_label = tk.Label(section_promoter_finder, text="Gene ID:")
gene_id_label.grid(row=0, column=0)
gene_id_entry = tk.Text(section_promoter_finder, height=5, width=20)
gene_id_entry.grid(row=1, column=0)
paste_button = tk.Button(section_promoter_finder, text="Paste", command=paste_gene)
paste_button.grid(row=2, column=0)


# Species selection
species_label = tk.Label(section_promoter_finder, text="Species:")
species_label.grid(row=3, column=0)
species_combobox = ttk.Combobox(section_promoter_finder, values=["Human", "Mouse", "Rat","Drosophila", "Zebrafish"])
species_combobox.current(0)
species_combobox.grid(row=4, column=0)

# Upstream/downstream entry
upstream_label = tk.Label(section_promoter_finder, text="Upstream (bp):")
upstream_label.grid(row=5, column=0)
upstream_entry = tk.Entry(section_promoter_finder)
upstream_entry.insert(2000, "2000")  # $"2000" default
upstream_entry.grid(row=6, column=0)

downstream_label = tk.Label(section_promoter_finder, text="Downstream (bp):")
downstream_label.grid(row=7, column=0)
downstream_entry = tk.Entry(section_promoter_finder)
downstream_entry.insert(500, "500")  # $"500" default
downstream_entry.grid(row=8, column=0)

# Search
search_button = tk.Button(section_promoter_finder, text="Find promoter  (CAN BE STUCK ! Don't worry, just wait ~30sec/gene)", command=get_sequence)
search_button.grid(row=9, column=0)

# Section "Promoter"
section_promoter = tk.LabelFrame(section_promoter_finder, text="Promoter Finder")
section_promoter.grid(row=10, column=0, padx=10, pady=10)

# Promoter output
result_promoter = tk.Text(section_promoter, height=10, width=50)
result_promoter.pack()
copy_button = tk.Button(section_promoter, text="Copy", command=copy_sequence)
copy_button.pack(side="left", fill="x", expand=True)
paste_button = tk.Button(section_promoter, text="Paste", command=paste_sequence)
paste_button.pack(side="left", fill="x", expand=True)

# Section help and GitHub

# Help & Credit
section_help = tk.LabelFrame(window, text="Help & Credit")
section_help.grid(row=1, column=0, padx=10, pady=10)

#How to use
help_button = tk.Button(section_help, text="How to use", command=show_help_PDF)
help_button.pack(side="left", fill="x", expand=True)

# WebUI
button = tk.Button(section_help, text="WebUI", command=open_site_webui)
button.pack(side="left", fill="x", expand=True)

# Github
button = tk.Button(section_help, text="Github by MINNITI Julien", command=open_site)
button.pack(side="left", fill="x", expand=True)

# Section "Responsive Elements finder"
section_responsive_finder = tk.LabelFrame(window, text="Responsive Elements Finder")
section_responsive_finder.grid(row=0, column=1, padx=10, pady=10)

# RE entry
label_sequence = tk.Label(section_responsive_finder, text="Responsive element (IUPAC authorize) :")
label_sequence.grid(row=1, column=0)
entry_sequence = tk.Entry(section_responsive_finder)
entry_sequence.insert(0, "ATGCCGTA")  # default
entry_sequence.grid(row=2, column=0)
use_jaspar = tk.IntVar()
checkbox_jaspar = tk.Checkbutton(section_responsive_finder, text="JASPAR", variable=use_jaspar)
checkbox_jaspar.grid(row=3, column=0)

# TIS entry
label_tis = tk.Label(section_responsive_finder, text="Transcription Initiation Site (bp)")
label_tis.grid(row=4, column=0)
label_tis_info = tk.Label(section_responsive_finder, text="(distance from start of promoter or 'Upstream' if you use the Promoter Finder)")
label_tis_info.grid(row=5, column=0)
entry_tis_var = tk.StringVar()
entry_tis_var.set("0")  # Valeur par défaut
entry_tis = tk.Entry(section_responsive_finder, textvariable=entry_tis_var)
entry_tis.grid(row=6, column=0)

# Threshold
threshold_label = tk.Label(section_responsive_finder, text="Threshold (%)")
threshold_label.grid(row=7, column=0)
threshold_entry = tk.Entry(section_responsive_finder)
threshold_entry.insert(80, "80")  # $"80" default
threshold_entry.grid(row=8, column=0)

# Find RE
button_search = tk.Button(section_responsive_finder, text="Find responsive elements", command=find_sequence_consensus)
button_search.grid(row=9, column=0)

# RE output
label_result = tk.Label(section_responsive_finder, text="Responsive elements")
label_result.grid(row=10, column=0)
text_result = tk.Text(section_responsive_finder, height=16, width=100)
text_result.grid(row=11, column=0)

# Export to Excel
export_button = tk.Button(section_responsive_finder, text="Export to Excel", command=export_to_excel)
export_button.grid(row=12, column=0)

# Section "status"
section_status = tk.LabelFrame(window, text="Status")
section_status.grid(row=1, column=1, padx=10, pady=10)

# status output
text_status = tk.Text(section_status, height=1, width=100)
text_status.grid(row=0, column=0)

# Configure grid weights
window.grid_rowconfigure(0, weight=1)
window.grid_columnconfigure(1, weight=1)
window.grid_columnconfigure(2, weight=1)
section_promoter_finder.grid_rowconfigure(13, weight=1)
section_responsive_finder.grid_rowconfigure(12, weight=1)

# Main loop
window.mainloop()