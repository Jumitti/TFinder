# Copyright (c) 2023 Minniti Julien

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of TFinder and associated documentation files, to deal
# in TFinder without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of TFinder, and to permit persons to whom TFinder is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of TFinder.

# TFINDER IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH TFINDER OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import io

import streamlit as st

from tfinder import IMO


def pwm_page():
    col1, col2 = st.columns(2)
    individual_motif = col1.text_area("🔹 :blue[**Step 2.3**] Sequences:",
                                      value=">seq1\nCTGCCGGAGGA\n>seq2\nAGGCCGGAGGC\n>seq3\nTCGCCGGAGAC\n>seq4\nCCGCCGGAGCG\n>seq5\nAGGCCGGATCG"
                                      if 'individual_motif_save' not in st.session_state else st.session_state[
                                          'individual_motif_save'], height=125,
                                      help='Put FASTA sequences. Same sequence length required ⚠')
    st.session_state['individual_motif_save'] = individual_motif
    individual_motif = individual_motif.upper()

    try:
        matrix, weblogo = IMO.individual_motif_pwm(individual_motif)
        matrix_str = ""
        for base, values in matrix.items():
            values_str = " ".join([f"{val:.4f}" for val in values])
            matrix_str += f"{base} [ {values_str} ]\n"
        matrix_text = col2.text_area('PWM', value=matrix_str, height=125,
                                     help='Copy to use later. Not editable.',
                                     disabled=True)
        st.pyplot(weblogo.fig)
        logo = io.BytesIO()
        weblogo.fig.savefig(logo, format='png')
        logo.seek(0)
        st.session_state['weblogo'] = logo
    except Exception as e:
        st.error(e)
