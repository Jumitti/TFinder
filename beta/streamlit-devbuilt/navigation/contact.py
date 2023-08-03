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

import streamlit as st

def contact_page():
    a, b = st.columns([0.9, 1.4])
    a.image('https://www.ipmc.cnrs.fr/fichiers/images/photo_1.jpg')
    a.markdown("""<span style="font-size:20px;">IPMC - Institut de Pharmacologie Moléculaire et Cellulaire</span>""", unsafe_allow_html=True)

    b.markdown("""<div style="width: 100%"><iframe src="https://www.google.com/maps/embed?pb=!1m18!1m12!1m3!1d2888.412189945306!2d7.052291366411674!3d43.61877920790173!2m3!1f0!2f0!3f0!3m2!1i1024!2i768!4f13.1!3m3!1m2!1s0x12cc2b06e4c3ff13%3A0x5790dfd8a00e8aac!2sIPMC%20-%20Institut%20de%20Pharmacologie%20Mol%C3%A9culaire%20et%20Cellulaire!5e0!3m2!1sfr!2sfr!4v1687462568391!5m2!1sfr!2sfr" width="700" height="450" style="border:0;" allowfullscreen="" loading="lazy" referrerpolicy="no-referrer-when-downgrade"></iframe></div>""", unsafe_allow_html=True) 
    
    a.markdown('##### If you have any questions or feedback, please feel free to contact us:')
    a.markdown("""<span style="font-size:18px;">PhD. Minniti Julien<br>PhD and Software developer<br>✉️ minniti@ipmc.cnrs.fr</span>""", unsafe_allow_html=True)
    a.markdown("""<span style="font-size:18px;">Dr. Duplan Eric<br>Research Engineer<br>✉️ duplan@ipmc.cnrs.fr</span>""", unsafe_allow_html=True)
    a.markdown("""<span style="font-size:18px;">Dr. Alves da Costa Cristine<br>Research Director<br>✉️ acosta@ipmc.cnrs.fr</span>""", unsafe_allow_html=True) 