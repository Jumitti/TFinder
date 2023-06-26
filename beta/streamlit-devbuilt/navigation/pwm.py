import streamlit as st
import numpy as np
from Bio import SeqIO
from weblogo import *

def pwm_page():
    def calculate_pwm(sequences):
        sequence_length = len(sequences[0])
        num_sequences = len(sequences)
        
        # Vérifier que toutes les séquences ont la même longueur
        for sequence in sequences[1:]:
            if len(sequence) != sequence_length:
                st.warning("Sequence lengths are not consistent.")
                return None
        
        pwm = np.zeros((4, sequence_length))
        for i in range(sequence_length):
            counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
            for sequence in sequences:
                nucleotide = sequence[i]
                if nucleotide in counts:
                    counts[nucleotide] += 1
            pwm[0, i] = counts['A'] / num_sequences * 100
            pwm[1, i] = counts['T'] / num_sequences * 100
            pwm[2, i] = counts['G'] / num_sequences * 100
            pwm[3, i] = counts['C'] / num_sequences * 100

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

    if st.button('Generate PWM'):
        if fasta_text:
            
            sequences = parse_fasta(fasta_text)
            sequences = [seq.upper() for seq in sequences]

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
            
            sequences = []
            for line in fasta_text.splitlines():
                if line.startswith(">"):
                    if sequences:
                        sequences.append(current_sequence)
                    current_sequence = ""
                else:
                    current_sequence += line
            if current_sequence:
                sequences.append(current_sequence)

            # Création de l'objet LogoData à partir des séquences
            data = LogoData.from_seqs(sequences)

            # Configuration des options du logo Web
            options = LogoOptions()
            options.title = "WebLogo"
            options.fineprint = "Logo generated using Biopython and weblogo"
            options.color_scheme = classic

            # Génération du logo Web
            format = LogoFormat(data, options)
            png = png_formatter(data, format)

            # Affichage du logo Web
            st.image(png, use_column_width=True)



