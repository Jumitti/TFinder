import streamlit as st
from weblogo import *
import json
from Bio import SeqIO
from io import StringIO

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
    
    json_text = st.text_area("Enter sequence matrix in JSON format", height=300)

    def convert_json_to_fasta(json_text):
        try:
            matrix = json.loads(json_text)
            nucleotides = list(matrix.keys())
            seq_len = len(matrix[nucleotides[0]])

            fasta_text = ""
            for i in range(seq_len):
                seq = ""
                for nuc in nucleotides:
                    seq += nuc * int(matrix[nuc][i])
                fasta_text += f">{i}\n{seq}\n"
            return fasta_text
        except Exception as e:
            st.warning("Invalid JSON format")
            return ""

    def generate_weblogo(fasta_text):
        seqs = SeqIO.parse(StringIO(fasta_text), 'fasta')
        seq_list = list(seqs)
        if len(seq_list) == 0 or len(seq_list[0]) == 0:
            st.warning("No sequences found in the input")
            return

        counts = [seq.counts() for seq in seq_list]
        data = LogoData.from_counts(seq_list[0].alphabet, counts)

        options = LogoOptions()
        options.title = "WebLogo"
        options.fineprint = "Logo generated using Biopython and weblogo"
        options.color_scheme = classic

        format = LogoFormat(data, options)
        png = png_formatter(data, format)

        return png

    if st.button('Generate PWM'):
        if json_text:
            fasta_text = convert_json_to_fasta(json_text)
            weblogo1 = generate_weblogo(fasta_text)
            st.image(weblogo1, use_column_width=True)

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



