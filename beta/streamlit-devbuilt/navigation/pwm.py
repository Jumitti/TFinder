import streamlit as st
from weblogo import *
import json
from Bio import SeqIO
from io import StringIO
import requests
from PIL import Image
from io import BytesIO
from weblogo import LogoData, LogoFormat, LogoOptions, pdf_formatter

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

            # Fonction pour convertir les séquences FASTA en une matrice
            def fasta_to_matrix(fasta_text):
                matrix = {}

                # Création d'un objet fichier virtuel à partir de la chaîne de texte FASTA
                fasta_file = StringIO(fasta_text)

                # Analyse des séquences FASTA
                records = SeqIO.parse(fasta_file, "fasta")
                for record in records:
                    sequence = str(record.seq)
                    matrix[record.id] = list(sequence)

                return matrix

            # Vérification de la présence des séquences FASTA
            if fasta_text != "":
                # Conversion des séquences FASTA en une matrice
                matrix = fasta_to_matrix(fasta_text)
            else:
                st.warning("Veuillez saisir les séquences FASTA.")

            # Génération du logo Web si la matrice de séquences est disponible
            if "matrix" in locals():
                # Conversion de la matrice en une liste de séquences
                sequences = ["".join(matrix[key]) for key in matrix]

                # Création des données du logo à partir des séquences
                data = LogoData.from_seqs(sequences)

                # Options de configuration du logo
                options = LogoOptions()
                options.title = "Logo Web"
                options.show_errorbars = False

                # Format du logo
                format_weblogo = LogoFormat(data, options)

                # Chemin de sortie du fichier PDF
                output_path = "logo.png"

                # Génération du logo au format PNG
                with open(output_path, "wb") as output_file:
                    output_format = format_weblogo.formatter(format_weblogo, format_weblogo)
                    output_format.write(output_file)

                # Affichage du logo dans Streamlit
                st.image(output_path)



