import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

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

                # Créer le logo de séquence en utilisant Matplotlib
                fig, ax = plt.subplots()
                logo_heights = np.sum(pwm, axis=0)
                for i, base in enumerate(bases):
                    heights = pwm[i]
                    y_start = np.sum(logo_heights[:i])
                    for j, height in enumerate(heights):
                        if np.isfinite(height):
                            rect = Polygon(
                                np.array(
                                    [
                                        [j, y_start],
                                        [j + 1, y_start],
                                        [j + 1, y_start + height],
                                        [j, y_start + height],
                                    ]
                                ),
                                facecolor=base,
                                edgecolor="black",
                            )
                            ax.add_patch(rect)

                ax.set_xlim(0, len(sequences[0]))
                ax.set_ylim(0, np.sum(logo_heights))
                ax.set_xticks(np.arange(len(sequences[0])) + 0.5)
                ax.set_xticklabels(range(1, len(sequences[0]) + 1))
                ax.set_ylabel("Bits")
                ax.set_title("Sequence Logo")

                st.pyplot(fig)

            else:
                st.warning("You forget FASTA sequences :)")

pwm_page()
