import streamlit as st

def contact_page():
    st.markdown("<h2 style='text-align: center; color: black;'>Contact Us</h1>", unsafe_allow_html=True)  
    
    a, b = st.columns([0.8, 1.5])
    a.image('https://www.ipmc.cnrs.fr/fichiers/images/photo_1.jpg')
    a.markdown("""<span style="font-size:20px;">University of Pittsburgh, UPMC Hillman Cancer Center, Assembly Building</span>""", unsafe_allow_html=True)

    b.markdown("""<div style="width: 100%"><iframe width="100%" height="700" frameborder="0" scrolling="no" marginheight="0" marginwidth="0" src="https://www.google.com/maps/place/IPMC+-+Institut+de+Pharmacologie+Mol%C3%A9culaire+et+Cellulaire/@43.6187171,7.0529351,15z/data=!4m14!1m7!3m6!1s0x12cc2b06e4c3ff13:0x5790dfd8a00e8aac!2sIPMC+-+Institut+de+Pharmacologie+Mol%C3%A9culaire+et+Cellulaire!8m2!3d43.6187171!4d7.0529351!16s%2Fg%2F11ddys_h2r!3m5!1s0x12cc2b06e4c3ff13:0x5790dfd8a00e8aac!8m2!3d43.6187171!4d7.0529351!16s%2Fg%2F11ddys_h2r?entry=ttu"><a href="https://www.maps.ie/distance-area-calculator.html">measure distance on map</a></iframe></div>""", unsafe_allow_html=True)
    # b.write('Hillman Cancer Center, 5051 Centre Avenue Pittsburgh, PA 15213') 
    
    a.markdown('##### If you have any questions or feedback, please feel free to contact us:')
    a.markdown("""<span style="font-size:18px;">Hatice Osmanbeyoglu<br>Principal Investigator<br>✉️ osmanbeyogluhu@pitt.edu</span>""", unsafe_allow_html=True)
    a.markdown("""<span style="font-size:18px;">Koushul Ramjattun<br>Graduate Researcher<br>✉️ kor11@pitt.edu</span>""", unsafe_allow_html=True)
    a.markdown("""<span style="font-size:18px;">Xiaojun Ma<br>Senior Software Engineer<br>✉️ xim33@pitt.edu</span>""", unsafe_allow_html=True)