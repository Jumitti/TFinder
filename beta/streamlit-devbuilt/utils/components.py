# Copyright (c) 2023 Minniti Julien

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from io import TextIOWrapper
import re
from typing import Tuple
import uuid
import numpy as np
import stringdb
from PIL import Image
import extra_streamlit_components as stx
import streamlit as st
from pyvis.network import Network
import streamlit.components.v1 as components
from matplotlib import colors
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


footer_style = f"""
    <style>
        MainMenu {{visibility: hidden;}}
        footer {{visibility: hidden;}}
        footer:after {{
            content:'Copyright @2023 Minniti Julien'; 
            visibility: visible;
            display: block;
            position: relative;
            # background-color: red;
            padding: 5px;
            top: 2px;
        }}
    </style>
"""


hide_menu_style = """
        <style>
        #MainMenu {visibility: hidden;}
        </style>
        """
# st.markdown(hide_menu_style, unsafe_allow_html=True)

icons = """<style>
        img {
            max-width: 100%;
        
        }
        
        .headerStyle {
            transition: transform .2s;
        }
        
        .headerStyle:hover {
            
             transform: scale(1.5);
            transition: 0.2s;
        }
        .image1 { 
            padding: 10px;
             transition: transform .2s;
        }
        .image2 
        { 
            padding: 10px;
             transition: transform .2s;
        }
        .image1:hover {  
            ##border: 4px solid green;
            ##border-radius: 15px;
             transform: scale(1.5);
            transition: 0.2s;
        }

        .image2:hover {  
            ##border: 4px solid green;
            ##border-radius: 15px;
             transform: scale(1.5);
            transition: 0.2s;
        }
        
        a:link,
        a:visited {
            color: blue;
            background-color: transparent;
            text-decoration: underline;
        }

        a2:hover {
            border-style: solid;
            },
        a:active {
            color: red;
            background-color: transparent;
            text-decoration: underline;
        }
    
        .footer {
            position: fixed;
            width: 100%;
            background-color: white;
            color: black;
            display: block;
            text-align: center;
            padding: 100px;
            top: 0px;
        }
</style
<div class="footer">
        </a>
        <a href="https://github.com/Jumitti"><img class="image2" src="https://cdn-icons-png.flaticon.com/512/25/25231.png"
                alt="github" width="60" height="60">
</div>
""" ## Footer
footer="""<style>
img {
            max-width: 100%;
        
        }
        
        .headerStyle {
            transition: transform .2s;
        }
        
        .headerStyle:hover {
            
             transform: scale(1.5);
            transition: 0.2s;
        }
        .image1 { 
            padding: 10px;
             transition: transform .2s;
        }
        .image2 
        { 
            padding: 10px;
             transition: transform .2s;
        }
        .image1:hover {  
            ##border: 4px solid green;
            ##border-radius: 15px;
             transform: scale(1.5);
            transition: 0.2s;
        }

        .image2:hover {  
            ##border: 4px solid green;
            ##border-radius: 15px;
             transform: scale(1.5);
            transition: 0.2s;
        }
        
        a:link,
        a:visited {
            color: blue;
            background-color: transparent;
            text-decoration: underline;
        }

        a2:hover {
            border-style: solid;
            },
        a:active {
            color: red;
            background-color: transparent;
            text-decoration: underline;
        }
.footer {
position: relative;
width: 100%;
left: 0;
bottom: 0;
background-color: transparent;
margin-top: auto;
color: #F66833;
padding: 24px;
text-align: center;
}
</style>
<div class="footer">
<p style="font-size:  13px">Copyright (c) 2023 Minniti Julien.</p>
<p style="font-size: 13px">This software is distributed under an MIT licence. Please consult the LICENSE file for more details.</p>
<a href="https://www.ipmc.cnrs.fr/cgi-bin/site.cgi"><img class="image2" src="https://www.ipmc.cnrs.fr/fichiers/images/logo_ipmc.jpg"alt="github" width="80" height=57"></a>
<a href="https://www.cnrs.fr/fr"><img class="image2" src="https://www.cnrs.fr/themes/custom/cnrs/logo.svg"alt="github" width="70" height="70"></a>
<a href="https://www.inserm.fr/"><img class="image2" src="https://www.sfr-biosciences.fr/medias/images/accueil/310842-logo-inserm.jpg"alt="github" width="70" height="70"></a>
<a href="https://univ-cotedazur.fr/"><img class="image2" src="https://www6.inrae.fr/var/internet6_national_consortium_biocontrole/storage/images/les-membres/annuaire-des-membres-du-consortium-biocontrole/universite-cote-d-azur/39475-1-fre-FR/Universite-Cote-d-Azur_inra_image.png"alt="github" width="70" height="70"></a>
</a><a href="https://github.com/Jumitti"><img class="image2" src="https://cdn-icons-png.flaticon.com/512/25/25231.png" alt="github" width="70" height="70"></a>
"""

# """        """

footer_style = """
    <style>
        MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
    </style>
"""

color_style = """
    <style>
    :root {
      --primary-color: blue;
    }
    </style>
"""