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
from PIL import Image
import datetime

def home_page():
    st.header('TFinder 🧬🔎')
    colhome1, colhome2 = st.columns ([1.7, 0.3])
    with colhome1:
        st.image('https://github.com/Jumitti/TFinder/blob/main/img/coordinates.png?raw=true', caption='Figure 8: Chromosic coordinates on a genetic map from NCBI')
    with colhome2:
        st.text('TFinder is a Python easy-to-use web tool for identification of putative Transcription Factor Binding Sites (TFBS) in a sequence. It allows extracting directly the promoter or terminal regions of a gene via the NCBI API for 5 different species, with no limit on the number of genes. The reference pattern (ex: a TFBS) accepts both IUPAC codes and JASPAR entries. It is also possible to generate and to use a Position Weight Matrix (PWM). Finally, the data are presented in tabular form, along with a graph showing the relevance of the TFBSs found as a function of their relative position on the sequence.')
    