import streamlit as st
import hydralit_components as hc
import requests
import pandas as pd
import altair as alt
import math
import pickle

from utils.components import footer_style, footer
from navigation.home import home_page
from navigation.REF import REF_page
from navigation.resource import resource_page
from navigation.contact import contact_page
from navigation.pwm import pwm_page

def application_page():
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
       {'icon': "🏠", 'label':AIO},
       {'icon': "🖥️", 'label':PROM},
       {'icon': "🖥️", 'label':REF},
       {'icon':"📑",'label':PWM}
       
    ]

    over_theme = {'txc_inactive': 'black','menu_background':'#ECECEC','txc_active':'white','option_active':'#fab947'}
    font_fmt = {'font-class':'h3','font-size':'50%'}

    chosen_tab = hc.option_bar(
        option_definition=option_data,
        title='',
        key='PrimaryOptionx',
        override_theme=over_theme,
        horizontal_orientation=True)

    if chosen_tab == AIO:
        REF_page()

    elif chosen_tab == PROM: 
        st.succes('PROM')

    elif chosen_tab == REF: 
        st.succes('REF'
        
    elif chosen_tab == PWM: 
        pwm_page()

