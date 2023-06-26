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

                # Analyse des séquences FASTA
                records = SeqIO.parse(fasta_text, "fasta")
                for record in records:
                    sequence = str(record.seq)
                    matrix[record.id] = list(sequence)

                return matrix

            # Lecture des séquences FASTA depuis un fichier ou une zone de texte
            fasta_file = st.file_uploader("Sélectionnez un fichier FASTA", type=["fasta"])
            fasta_text = st.text_area("Saisissez les séquences FASTA")

            # Vérification de la présence des séquences FASTA
            if fasta_file is not None:
                # Lecture des séquences FASTA depuis le fichier
                matrix = fasta_to_matrix(fasta_file)
            elif fasta_text != "":
                # Lecture des séquences FASTA depuis la zone de texte
                matrix = fasta_to_matrix(fasta_text)
            else:
                st.warning("Veuillez charger un fichier FASTA ou saisir les séquences FASTA.")

            # Génération du logo Web si la matrice de séquences est disponible
            if "matrix" in locals():
                # Création des données du logo à partir de la matrice de séquences
                data = LogoData.from_counts(matrix)

                # Options de configuration du logo
                options = LogoOptions()
                options.logo_title = "Logo Web"
                options.show_errorbars = False

                # Format du logo
                format = LogoFormat(data, options)

                # Chemin de sortie du fichier PDF
                output_path = "logo.pdf"

                # Génération du logo au format PDF
                with open(output_path, "wb") as output_file:
                    pdf_formatter(data, format, output_file)

                # Affichage du logo dans Streamlit
                st.image(output_path)



