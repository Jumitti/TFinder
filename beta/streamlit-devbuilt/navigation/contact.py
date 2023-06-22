import streamlit as st

def contact_page():
    st.markdown("<h2 style='text-align: center; color: black;'>Contact Us</h1>", unsafe_allow_html=True)  
    
    a, b = st.columns([0.8, 1.5])
    a.image('https://www.ipmc.cnrs.fr/fichiers/images/photo_1.jpg')
    a.markdown("""<span style="font-size:20px;">University of Pittsburgh, UPMC Hillman Cancer Center, Assembly Building</span>""", unsafe_allow_html=True)

    b.markdown("""<div style="width: 100%"><iframe width="100%" height="700" frameborder="0" scrolling="no" marginheight="0" marginwidth="0" src="https://maps.google.com/maps?width=100%25&amp;height=600&amp;hl=en&amp;q=5051%20Centre%20Avenue%20Pittsburgh,%20PA%2015213+(My%20Business%20Name)&amp;t=p&amp;z=14&amp;ie=UTF8&amp;iwloc=B&amp;output=embed"><a href="https://www.maps.ie/distance-area-calculator.html">measure distance on map</a></iframe></div>""", unsafe_allow_html=True)
    # b.write('Hillman Cancer Center, 5051 Centre Avenue Pittsburgh, PA 15213') 
    
    a.markdown('##### If you have any questions or feedback, please feel free to contact us:')
    a.markdown("""<span style="font-size:18px;">Hatice Osmanbeyoglu<br>Principal Investigator<br>✉️ osmanbeyogluhu@pitt.edu</span>""", unsafe_allow_html=True)
    a.markdown("""<span style="font-size:18px;">Koushul Ramjattun<br>Graduate Researcher<br>✉️ kor11@pitt.edu</span>""", unsafe_allow_html=True)
    a.markdown("""<span style="font-size:18px;">Xiaojun Ma<br>Senior Software Engineer<br>✉️ xim33@pitt.edu</span>""", unsafe_allow_html=True)