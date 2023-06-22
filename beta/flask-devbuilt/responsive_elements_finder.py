from flask import Flask
from flask import render_template
import requests
import pandas as pd
import altair as alt
import math
import pickle

app = Flask(__name__, template_folder='templates', static_folder='static')

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/process', methods=['POST'])
def process():
    gene_id_entry = request.form['gene_id_entry']
    species_combobox = request.form['species_combobox']
    prom_term = request.form['prom_term']
    upstream_entry = request.form['upstream_entry']
    downstream_entry = request.form['downstream_entry']

    # Run Promoter Finder
    if prom_term == 'Promoter':
        try:
            result_promoter = find_promoters(gene_id_entry, species_combobox, int(upstream_entry), int(downstream_entry))
            message = "Promoters extraction complete!"
            return render_template('result.html', result_promoter=result_promoter, message=message)
        except Exception as e:
            return render_template('error.html', error=str(e))
    else:
        try:
            result_promoter = find_promoters(gene_id_entry, species_combobox, int(upstream_entry), int(downstream_entry))
            message = "Terminators extraction complete!"
            return render_template('result.html', result_promoter=result_promoter, message=message)
        except Exception as e:
            return render_template('error.html', error=str(e))

if __name__ == '__main__':
    app.run()