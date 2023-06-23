import streamlit as st
import numpy as np

def pwm_page():
    def calculate_pwm(sequences):
        num_sequences = len(sequences)
        sequence_length = len(sequences[0])

        # Initialiser la matrice PWM avec des zéros
        pwm = np.zeros((4, sequence_length))

        # Compter les occurrences de chaque nucléotide pour chaque position
        for i in range(sequence_length):
            counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
            for sequence in sequences:
                nucleotide = sequence[i]
                if nucleotide in counts:
                    counts[nucleotide] += 1

            # Calculer les fréquences relatives des nucléotides
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

    # Streamlit app
    st.title("Calculatrice de PWM")

    fasta_text = st.text_area("Saisir les séquences au format FASTA", height=300)

    if fasta_text:
        sequences = parse_fasta(fasta_text)
        sequences = [seq.upper() for seq in sequences]

        if len(sequences) > 0:
            pwm = calculate_pwm(sequences)

            st.header("PWM résultante")

            # Afficher la matrice au format souhaité
            matrix = {
                'A': list(pwm[0]),
                'T': list(pwm[1]),
                'G': list(pwm[2]),
                'C': list(pwm[3])
            }

            for base, values in matrix.items():
                st.write(base + ":")
                st.write(values)

        else:
            st.warning("Aucune séquence valide n'a été trouvée.")

