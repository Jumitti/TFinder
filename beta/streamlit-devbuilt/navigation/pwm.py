import streamlit as st
import numpy as np
from weblogo import *
from Bio import motifs
import matplotlib.pyplot as plt
from io import StringIO

def pwm_page():
    def calculate_pwm(sequences):
        num_sequences = len(sequences)
        sequence_length = len(sequences[0])
        pwm = np.zeros((4, sequence_length))
        for i in range(sequence_length):
            counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
            for sequence in sequences:
                nucleotide = sequence[i]
                if nucleotide in counts:
                    counts[nucleotide] += 1
            pwm[0, i] = counts['A'] / num_sequences *100
            pwm[1, i] = counts['T'] / num_sequences *100
            pwm[2, i] = counts['G'] / num_sequences *100
            pwm[3, i] = counts['C'] / num_sequences *100

        return pwm

    def parse_fasta(fasta_text):
        sequences = []
        current_sequence = ""

        for line in fasta_text.splitlines():
            if line.startswith(">"):
                if current_sequence:
                    sequences.append(current_sequence)
                current_sequence = ""
            else:
                current_sequence += line

        if current_sequence:
            sequences.append(current_sequence)

        return sequences

    st.subheader("🧮 PWM generator")

    fasta_text = st.text_area("Put FASTA sequences. Same sequence length required ⚠️", height=300)
    
    # Fonction pour convertir les séquences au format FASTA en objet Bio.motifs.Motif
    def parse_sequences(sequences):
        records = []

        current_sequence_name = None
        current_sequence_instances = []

        for line in sequences:
            if line.startswith(">"):
                if current_sequence_name is not None:
                    current_sequence = motifs.create(current_sequence_instances, alphabet=motifs.IUPAC.unambiguous_dna, name=current_sequence_name)
                    records.append(current_sequence)
                    current_sequence_instances = []
                current_sequence_name = line[1:]
            else:
                current_sequence_instances.append(line)

        if current_sequence_name is not None:
            current_sequence = motifs.create(current_sequence_instances, alphabet=motifs.IUPAC.unambiguous_dna, name=current_sequence_name)
            records.append(current_sequence)

        return records
        
    # Fonction pour générer le weblogo à partir des motifs
    def generate_weblogo(motifs):
        m = motifs.combine(motifs)
        handle = StringIO()
        motifs.write(m, handle, format='jaspar')
        handle.seek(0)

        logo = motifs.read(handle, format='jaspar')
        logo.format = 'png'

        return logo

    if st.button('Generate PWM'):
        if fasta_text:
            
            # Votre input de séquences
            sequences_input = [
                '>p53',
                'CTGCCGGAGGA',
                '>PS1',
                'AGGCCGGAGGC',
                '>PS2',
                'TCGCCGGAGAC',
                '>CCNA2',
                'CCGCCGGAGCG',
                '>CCNB1',
                'AGGCCGGATCG'
            ]
            
            sequences = parse_fasta(fasta_text)
            sequences = [seq.upper() for seq in sequences]
            
            # Conversion des séquences en objets Motif
            sequences_weblogo = parse_sequences(fasta_text)
            
            # Génération du weblogo
            weblogo = generate_weblogo(sequences_weblogo)
            
            # Affichage du weblogo
            plt.figure(figsize=(6, 3))
            weblogo.plot(ax=plt.gca(), ic_scale=False)
            plt.xticks([])
            plt.yticks([])
            plt.show()

            if len(sequences) > 0:
                pwm = calculate_pwm(sequences)

                st.subheader("PWM: ")
                st.info("⬇️ Select and copy")
                bases = ['A', 'T', 'G', 'C']
                pwm_text = ""
                for i in range(len(pwm)):
                    base_name = bases[i]
                    base_values = pwm[i]

                    base_str = base_name + " ["
                    for value in base_values:
                        base_str += "\t" + format(value) + "\t" if np.isfinite(value) else "\t" + "NA" + "\t"

                    base_str += "]\n"
                    pwm_text += base_str

                st.text_area("PWM résultante", value=pwm_text)
                
                

            else:
                st.warning("You forget FASTA sequences :)")



