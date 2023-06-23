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
PROM = 'Promoter Extractor'
REF = 'REF'
PWM = 'PWM'

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

over_theme = {'txc_inactive': 'black','menu_background':'#ECECEC','txc_active':'white','option_active':'#fab947'}
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
        
    for i in range(6):
        st.markdown('#')
        

