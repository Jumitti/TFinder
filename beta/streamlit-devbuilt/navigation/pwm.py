import streamlit as st
from weblogo import *
import json
from Bio import SeqIO
from io import StringIO
import requests
from PIL import Image
from io import BytesIO
from weblogo import LogoData, LogoFormat, LogoOptions, pdf_formatter
import weblogolib
import corebio
import os

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

            def create_web_logo(output_file=None, out_format=None, title=None, units='probability',
                                alphabet=corebio.seq.unambiguous_protein_alphabet):
                fasta_input = st.text_input("Séquences (format FASTA)", "")

                if fasta_input == "":
                    st.warning("Veuillez entrer les séquences au format FASTA")
                    return

                sequences = []
                current_sequence = ""
                current_header = ""

                for line in fasta_input.splitlines():
                    if line.startswith(">"):
                        if current_sequence != "":
                            sequences.append((current_header, current_sequence))
                        current_header = line[1:]
                        current_sequence = ""
                    else:
                        current_sequence += line

                if current_sequence != "":
                    sequences.append((current_header, current_sequence))

                if out_format is None:
                    extension = os.path.splitext(output_file)[1] if output_file is not None else ''
                    out_format = extension[1:] if extension else 'png'

                seqs = corebio.seq.SeqList([corebio.seq.Seq(s, alphabet) for _, s in sequences], alphabet)
                seqs.alphabet = alphabet
                data = weblogolib.LogoData.from_seqs(seqs)
                options = weblogolib.LogoOptions()
                if title is not None:
                    options.logo_title = title
                options.unit_name = units
                options.show_fineprint = False
                
                if out_format == 'png':
                    options.resolution = 400.0
                
                format = weblogolib.LogoFormat(data, options)
                
                formatters = {
                    'png': weblogolib.png_formatter,
                    'svg': weblogolib.svg_formatter
                }
                
                image = formatters[out_format](data, format)
                if output_file is None:
                    return image
                
                with open(output_file, 'wb') as fout:
                    fout.write(image)

            # Utilisation de la fonction create_web_logo
            create_web_logo()



