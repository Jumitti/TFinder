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


def get_ppi_edge_list(gene, neighbors, size):
    genes = [gene]
    string_ids = stringdb.get_string_ids(genes)
    enrichment_df = stringdb.get_interaction_partners(string_ids.queryItem, limit=size)
    extras = list(enrichment_df['preferredName_B'].values[:neighbors])
    string_ids = stringdb.get_string_ids(extras+genes)
    enrichment_df = stringdb.get_interaction_partners(string_ids.queryItem, limit=size)

    return string_ids, enrichment_df[['preferredName_A',	'preferredName_B', 'score']].values


def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, color='black', fontsize=35,
            ha="left", va="center", transform=ax.transAxes)
    ax.tick_params(axis='x', labelsize=30)
    
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})

## https://matplotlib.org/stable/gallery/color/named_colors.html
palette_mapper = {
    'healthy': colors.to_hex('seagreen'),
    'mild': colors.to_hex('salmon'),
    'moderate': colors.to_hex('lightskyblue'),
    'critical': colors.to_hex('red'),
    'severe': colors.to_hex('magenta'),
    'stable': colors.to_hex('orange'),
    'progressive': colors.to_hex('blueviolet'),
}

severity_order = [
    'healthy',
    'stable',
    'mild',
    'progressive',
    'moderate',
    'severe',
    'critical'
]

def st_url(text, link, tags='######'):
    st.write(tags+f" [{text}]({link})")

def totags(tags):
    t = '<div>'
    for tag, color in zip(tags, ['red', 'blue', 'green', 'orange']):
        t += f"<span class='highlight {color}'>{tag}</span>\t"
    t += '</div>'
    return t

def make_palette(column, groups):
    return [colors.to_rgb(column.color_picker(g, palette_mapper.get(g.lower(), colors.to_hex('black')))) for g in groups]    


def fig2img(fig):
    """Convert a Matplotlib figure to a PIL Image and return it"""
    import io
    buf = io.BytesIO()
    fig.savefig(buf)
    buf.seek(0)
    img = Image.open(buf)
    return img

def ridge_plot(df_plt: pd.DataFrame, xlabel: str, column: st.container, group: str='patient_group'):
    
    group_sizes = df_plt.groupby('patient_group').size() 
    df_plt = df_plt.loc[df_plt['patient_group'].isin(group_sizes[group_sizes > 1].index)]
    
    
    groups_ = np.sort(np.unique(df_plt[group].values))
    groups_ = np.array([i.lower() for i in groups_])
    
    groups = []
    for g in severity_order:
        if g in groups_:
            groups.append(g)
        
    df_plt['patient_group'] = df_plt['patient_group'].str.lower()
    df_plt['patient_group'] = pd.Categorical(df_plt['patient_group'], categories=groups, ordered=True) 
    df_plt = pd.concat([group.sort_values('Pair_correlations', ascending=False) for name, group in df_plt.groupby('patient_group')])
    
    df_plt['y'] = 0
                    
    if column is not None:
        palette = make_palette(column=column, groups=groups)
    else:
        palette = sns.color_palette('pastel')
        
        
    g = sns.FacetGrid(df_plt, 
        palette = palette,
        row=group, hue=group, 
        aspect=5)
        
    g.map_dataframe(sns.kdeplot, x="Pair_correlations", fill=True, alpha=1)
    g.map_dataframe(sns.kdeplot, x="Pair_correlations", color='black')
    
    # g.map_dataframe(sns.scatterplot, y='y', x="Pair_correlations", color='black', s=150, marker='o')
    
    g.map(label, group)
    g.fig.subplots_adjust(hspace=-.25)
    
    g.set_titles("")
    g.set(yticks=[], xlabel=xlabel)
    g.despine( left=True)
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    plt.xlabel(xlabel, fontsize=36, labelpad=25)
    plt.gcf().subplots_adjust(bottom=0.15)
    
    plt.tight_layout()
    img = fig2img(g.fig)
    
    return img

def create_st_button(link_text, link_url, hover_color="#4EC5F1", st_col=None):    
    button_uuid = uuid.uuid1()
    button_id = re.sub("\d+", "", str(button_uuid))
    
    button_css = f"""
        <style>
            #{button_id} {{
                background-color: none;
                color: rgb(38, 39, 48);
                padding: 0.25em 0.38em;
                position: relative;
                text-decoration: none;
                border-radius: 4px;
                border-width: 2px;
                border-style: solid;
                border-color: rgb(13, 242, 201);
                border-image: initial;
            }}
            #{button_id}:hover {{
                border-color: {hover_color};
                color: {hover_color};
            }}
            #{button_id}:active {{
                box-shadow: none;
                background-color: {hover_color};
                color: white;
                }}
        </style> """

    html_str = f'<a href="{link_url}" target="_blank" id="{button_id}";>{link_text}</a><br></br>'

    if st_col is None:
        st.markdown(button_css + html_str, unsafe_allow_html=True)
    else:
        st_col.markdown(button_css + html_str, unsafe_allow_html=True)


def make_tab(name, desc=''):
    return stx.TabBarItemData(id=name, title=name, description=desc)


def display_ppi_data(option: str, disabled: bool = False) -> Tuple[TextIOWrapper, pd.DataFrame]:
    if not disabled:
        st.markdown(f'#### Protein-Protein Interactions for {option} from StringDB')
        st.caption('Threshold of significance to include an interaction: 400')
        st.caption('Learn more: https://string-db.org/help/api/')
        left_col, right_col = st.columns(2)
        size = left_col.slider('Interactors', 2, 100, 10)
        neighbors = right_col.slider('Neighbors', 1, 10, 3)
    else:
        size = 10
        neighbors = 3
        

    net = Network()
        
    try:
        caption, pairs = get_ppi_edge_list(option, neighbors, size)
        preferred_name = caption[caption['queryItem']==option].preferredName.values[0]

        for i, j, v in pairs:
            net.add_node(j, label=j)
            
        for n in net.nodes:
            if n['label'] == preferred_name:
                n['color'] = 'red'
                 
        for pair in pairs:
            net.add_edge(pair[0], pair[1], value=pair[2]/10)
    except:
        if not disabled: st.warning('No PPI Information available from stringdb')
        return None, None
    # net.show_buttons(filter_=['physics'])
    
    net.show('./utils/ppi.html')
        
    HtmlFile = open("./utils/ppi.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    if not disabled: components.html(source_code, height = 650)
    
    
    if not disabled: 
        st.markdown(f'#### Pathways & drugs associated with {option} derived from Wiki-CORONA')
    
    try:
        df = pd.read_html(f'http://severus.dbmi.pitt.edu/corona/index.php/search?q={option}')[0]
        df.columns = ['Interactant Symbol', 'Name', 'Associated Pathways', 'Binding Drugs', 'Associated Diseases']
        if not disabled: st.dataframe(df)
    except:
        if not disabled: st.warning('No interactors found from Wiki-CORONA')
        df = None
        
        
    return HtmlFile, df
    
def download_button(desc, data, left_col, right_col, label='Download as CSV'):
    right_col.download_button(
        label=label,
        data=data,
        file_name='projD.csv',
        mime='text/csv',
        key=uuid.uuid1()
    )
        
    left_col.write(desc)
    
    
def local_css(file_name='./assets/style.css'):
    with open(file_name) as f:
        st.markdown('<style>{}</style>'.format(f.read()), unsafe_allow_html=True)
        
        
def resize(img: Image, factor: float = 2) -> Image:
    x, y = img.size
    x = int(x/factor)
    y = int(y/factor)
    img = img.resize((x, y))
    return img