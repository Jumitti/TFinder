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

footer_style = f"""
    <style>
        MainMenu {{visibility: hidden;}}
        footer {{visibility: hidden;}}
        footer:after {{ 
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
</style>
"""

footer = """<style>
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
color: #163172;
padding: 24px;
text-align: center;
}
</style>
<div class="footer">
<p style="font-size:  13px">Copyright (c) 2023 Minniti Julien.</p>
<p style="font-size: 13px">This software is distributed under an MIT licence. Please consult the LICENSE file for more details.</p>
<p style="font-size: 13px">Please cite as: "Julien Minniti, Eric Duplan, Frédéric Checler et al. TFinder: a Python web tool for predicting Transcription Factor Binding Sites, PREPRINT (Version 1) available at Research Square <a href="https://doi.org/10.21203/rs.3.rs-3782387/v1">[10.21203/rs.3.rs-3782387/v1]</a>"</p>
<a href="https://www.ipmc.cnrs.fr/cgi-bin/site.cgi"><img class="image2" src="https://raw.githubusercontent.com/Jumitti/TFinder/main/img/logo%20ipmc.png"alt="github" width="80" height=57"></a>
<a href="https://www.cnrs.fr/fr"><img class="image2" src="https://www.cnrs.fr/themes/custom/cnrs/logo.svg"alt="github" width="70" height="70"></a>
<a href="https://www.inserm.fr/"><img class="image2" src="https://raw.githubusercontent.com/Jumitti/TFinder/main/img/logo%20inserm.png"alt="github" width="100" height="70"></a>
<a href="https://univ-cotedazur.fr/"><img class="image2" src="https://raw.githubusercontent.com/Jumitti/TFinder/main/img/logo%20univ%20cote%20dazur.png"alt="github" width="70" height="70"></a>
</a><a href="https://github.com/Jumitti"><img class="image2" src="https://cdn-icons-png.flaticon.com/512/25/25231.png" alt="github" width="70" height="70"></a>
"""

color_style = """
    <style>
    :root {
      --primary-color: blue;
    }
    </style>
"""