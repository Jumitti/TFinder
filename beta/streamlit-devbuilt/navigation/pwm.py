import streamlit as st
import numpy as np

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

    st.title("Calculatrice de PWM")

    fasta_text = st.text_area("Put FASTA sequences. Same sequence length required ⚠️", height=300)
    
    if st.button('Generate PWM'):
        if fasta_text:
            sequences = parse_fasta(fasta_text)
            sequences = [seq.upper() for seq in sequences]

            if len(sequences) > 0:
                pwm = calculate_pwm(sequences)

                st.header("PWM: ")
                st.info("⬇️ Select and copy")
                bases = ['A', 'T', 'G', 'C']
                for i in range(len(pwm)):
                    base_name = bases[i]
                    base_values = pwm[i]

                    base_str = base_name + " ["
                    for value in base_values:
                        base_str += "\t" + format(value) + "\t" if np.isfinite(value) else "\t" + "NA" + "\t"

                    base_str += "]"
                    st.write(base_str)

            else:
                st.warning("Aucune séquence valide n'a été trouvée.")
