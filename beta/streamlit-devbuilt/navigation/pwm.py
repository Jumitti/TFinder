import streamlit as st
import numpy as np

def pwm_page():
    def calculate_pwm(sequences):
        sequence_length = len(sequences[0])
        num_sequences = len(sequences)
        base_counts = np.zeros((4, sequence_length))

        # Count occurrences of each base at each position
        for sequence in sequences:
            for i, base in enumerate(sequence):
                if base == 'A':
                    base_counts[0][i] += 1
                elif base == 'T':
                    base_counts[1][i] += 1
                elif base == 'G':
                    base_counts[2][i] += 1
                elif base == 'C':
                    base_counts[3][i] += 1

        # Normalize the counts to frequencies
        base_frequencies = base_counts / num_sequences

        # Calculate the position weight matrix
        pwm = np.log2(base_frequencies)

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
            bases = ['A', 'C', 'G', 'T']
            for i in range(len(pwm)):
                base_name = bases[i]
                base_values = pwm[i]

                base_str = base_name + " ["
                for value in base_values:
                    base_str += "\t" + str(int(value)) + "\t" if np.isfinite(value) else "\t" + "NA" + "\t"


                base_str += "]"
                st.write(base_str)
        else:
            st.warning("Aucune séquence valide n'a été trouvée.")
