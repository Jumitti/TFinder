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
import hydralit_components as hc
import requests
import pandas as pd
import altair as alt
import math
import pickle

from navigation.aio import aio_page
from navigation.pwm import pwm_page
from navigation.BSF import BSF_page
from navigation.prom_extractor import prom_extractor_page

#NavBar

AIO = 'All in One'
PROM = 'Region Extractor'
REF = 'Binding Sites Finder'
PWM = 'Position Weight Matrix'

tabs = [
    AIO,
    PROM,
    REF,
    PWM
]

option_data = [
   {'icon': "🧬🔎🧮", 'label':AIO},
   {'icon': "🧬", 'label':PROM},
   {'icon': "🔎", 'label':REF},
   {'icon':"🧮",'label':PWM}
   
]

over_theme = {'txc_inactive': 'black','menu_background':'#CEE6F3','txc_active':'white','option_active':'#91C8E4'}
font_fmt = {'font-class':'h3','font-size':'50%'}

def allapp_page():
    page = hc.option_bar(
        option_definition=option_data,
        title='',
        override_theme=over_theme,
        horizontal_orientation=True)

    if page == AIO:
        aio_page()

    elif page == PROM: 
        prom_extractor_page()

    elif page == REF: 
        BSF_page()
        
    elif page == PWM: 
        pwm_page()
        
    for i in range(4):
        st.markdown('#')
        

