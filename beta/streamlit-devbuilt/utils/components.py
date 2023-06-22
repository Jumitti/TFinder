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
            content:'Copyright @2023 www.osmanbeyoglulab.com'; 
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
        <a href="https://scholar.google.com/citations?user=YzCsmdgAAAAJ&hl=en&inst=7044483460239389945"><img class="image2" src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/c7/Google_Scholar_logo.svg/512px-Google_Scholar_logo.svg.png"
                alt="github" width="60" height="60">
        </a>
        <a href="https://www.instagram.com/osmanbeyoglulab/"><img class="image2" src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/e7/Instagram_logo_2016.svg/132px-Instagram_logo_2016.svg.png"
                alt="github" width="60" height="60">
        </a>
        <a href="https://github.com/osmanbeyoglulab"><img class="image2" src="https://cdn-icons-png.flaticon.com/512/25/25231.png"
                alt="github" width="60" height="60">
        </a>
        <a href="https://twitter.com/hosmanbeyoglu?lang=en"><img class="image2" src="https://upload.wikimedia.org/wikipedia/commons/f/f2/Logo_Twitter.png"
            alt="twitter" width="65" height="60"></a>
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
background-color: white;
margin-top: auto;
color: black;
padding: 24px;
text-align: center;
}
</style>
<div class="footer">
<p style="font-size:  13px">© 2023 Osmanbeyoglulab.com. All rights reserved.</p>
<a href="https://hillman.upmc.com/"><img class="image2" src="https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcS7c7pXIkFgMPVM2csVE6MenUFLEsgF5beCeMJzogkPkXPC4xEo9OTHf6nVpqsb3PvisRk&usqp=CAU"alt="github" width="70" height=50"></a>
<a href="https://www.pitt.edu/"><img class="image2" src="https://upload.wikimedia.org/wikipedia/en/thumb/f/fb/University_of_Pittsburgh_seal.svg/300px-University_of_Pittsburgh_seal.svg.png"alt="github" width="45" height="45"></a>
</a><a href="https://github.com/osmanbeyoglulab"><img class="image2" src="https://cdn-icons-png.flaticon.com/512/25/25231.png" alt="github" width="45" height="45"></a>
<a href="https://twitter.com/hosmanbeyoglu?lang=en"><img class="image2" src="https://upload.wikimedia.org/wikipedia/commons/6/6f/Logo_of_Twitter.svg"alt="twitter" width="45" height="40"></a>
<a href="https://scholar.google.com/citations?user=YzCsmdgAAAAJ&hl=en&inst=7044483460239389945"><img class="image2" src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/c7/Google_Scholar_logo.svg/512px-Google_Scholar_logo.svg.png"alt="github" width="45" height="45"></a>
<a href="https://www.instagram.com/osmanbeyoglulab/"><img class="image2" src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/e7/Instagram_logo_2016.svg/132px-Instagram_logo_2016.svg.png"alt="github" width="45" height="45"></a>
</div>
"""

# """        """

footer_style = """
    <style>
        MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
    </style>
"""