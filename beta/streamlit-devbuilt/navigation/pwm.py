import streamlit as st
from weblogo import LogoData, LogoFormat, LogoOptions
from io import StringIO

def pwm_page():
    def calculate_pwm(sequences):
        num_sequences = len(sequences)
        sequence_length = len(sequences[0])
        pwm = [[0, 0, 0, 0] for _ in range(sequence_length)]

        for i in range(sequence_length):
            counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
            for sequence in sequences:
                nucleotide = sequence[i]
                if nucleotide in counts:
                    counts[nucleotide] += 1
            pwm[i][0] = counts['A'] / num_sequences
            pwm[i][1] = counts['T'] / num_sequences
            pwm[i][2] = counts['G'] / num_sequences
            pwm[i][3] = counts['C'] / num_sequences

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
                        base_str += "\t" + format(value) + "\t" if value != 0 else "\t" + "NA" + "\t"

                    base_str += "]\n"
                    pwm_text += base_str

                st.text_area("PWM résultante", value=pwm_text)

                # Générer le logo de séquence avec weblogo
                pwm_data = LogoData.from_counts("ACGT", pwm)
                format_options = LogoFormat(pwm_data, LogoOptions())
                logo_format = StringIO()
                format_options.write(logo_format, "png")
                logo_format.seek(0)

                st.image(logo_format, use_column_width=True)

            else:
                st.warning("You forget FASTA sequences :)")