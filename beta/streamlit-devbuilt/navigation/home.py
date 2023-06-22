import streamlit as st
from PIL import Image
import datetime

import streamlit_analytics


def load_home_image():
    return Image.open('./assets/images/home.png')

TITLE = 'COVID-19db linkage maps of cell surface proteins and transcription factors in immune cells'



def home_page():

    from texts.descriptions import Desc
    img = load_home_image()
    
    _, middle, _ = st.columns([0.1, 1, 0.1])
    
    with middle:
        
        streamlit_analytics.start_tracking()    
        # st.markdown('### ' + TITLE) 
        st.markdown(f"<h3 style='text-align: center; color: black;'>{TITLE}</h1>", unsafe_allow_html=True)  
        
        streamlit_analytics.stop_tracking()
        
        
        with open('./clock.time', 'r') as f:
            last_updated_on = f.readlines()[0]
        
        st.caption(last_updated_on)    
        st.image(img)  
        
        st.markdown('👉 Read our paper here: https://www.biorxiv.org/content/10.1101/2022.12.14.520411v1')    
        
        
        st.markdown(f"<p style='text-align: justify; color: black;'>{Desc.home}</h4>", unsafe_allow_html=True)  
        
        
        # st.write(Desc.home)
        
        # st.markdown('---')
        
        # st.caption(f'Example usage of the {TITLE}')
        # video_file = open('./assets/video/video.mp4', 'rb')
        # video_bytes = video_file.read()
        # st.video(video_bytes)
    
        
if __name__ == '__main__':
    t = datetime.datetime.now()
    with open('./clock.time', 'w') as f:
        f.write(f'Last updated at {t.time():%H:%M} on {t.date():%Y-%m-%d}\n')
        print(f'Last updated at {t.time():%H:%M} on {t.date():%Y-%m-%d}\n')
    
    
    