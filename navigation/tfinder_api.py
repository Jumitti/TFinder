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


def tfinder_api():
    st.divider()
    st.markdown("<h3 style='text-align: center; color: black;'>Documentation for TFinder python package</h1>",
                unsafe_allow_html=True)
    st.markdown('How to use:')

    st.divider()
    st.markdown("<h3 style='text-align: center; color: black;'>DNA Region Extractor</h1>",
                unsafe_allow_html=True)

    with st.expander('Analyse gene availability'):
        st.markdown(
            'Some genes do not have the same name in different species. It can also happen that the gene ID is incorrect.')
        analyse_gene = 'NCBIdna.analyse_gene(gene_id)'
        st.code(analyse_gene)

        st.divider()
        st.markdown('**Parameter**')
        st.markdown('**gene_id** (list):')
        st.markdown(
            'Support only one by one gene_id. See example. Analyse if ID GENE is valid or if NAME GENE exist for Human, Mouse, Rat, Drosophila, Zebrafish')

        st.divider()
        st.markdown('**Return**')
        st.markdown('(list)')
        st.markdown('Order of results:')
        st.markdown('gene_id | ID | Human | Mouse | Rat | Drosophila | Zebrafish')

        st.divider()
        st.markdown('**Example**')
        st.markdown('Code:')
        example = '''gene_ids = ['4843', 'PRKN']
results = []

for gene_id in gene_ids:
    results.append(NCBIdna.analyse_gene(gene_id))

print(results)'''
        st.code(example)
        st.markdown('Result:')
        st.code("['4843', '✅', 'n.d', 'n.d', 'n.d', 'n.d', 'n.d'], ['PRKN', 'n.d', '✅', '✅', '✅', '❌', '✅']")

    with st.expander('Extraction of regulatory gene sequences (All in one)'):
        st.markdown('Recovery of promoter or terminal sequences of a given gene for 5 different species.')
        find_sequences = 'NCBIdna(gene_id, prom_term, upstream, downstream, species=None).find_sequences()'
        st.code(find_sequences)

        st.divider()
        st.markdown('**Parameters**')
        st.markdown('**gene_id** (list):')
        st.markdown('It should be **gene name** or **gene ID**. Support only one by one gene_id.')
        st.markdown('**prom_term** ("promoter" or "terminator"):')
        st.markdown(
            'Choose between "promoter" or "terminator". It defines origin of extraction (TSS for promoter and gene end for terminator). Chromosomal coordinates are display in results')
        st.markdown('**upstream** (integer):')
        st.markdown('Upstream defines the number of base pairs (bp) before the TSS/gene end ')
        st.markdown('**downstream** (integer):')
        st.markdown('Downstream defines the number of base pairs (bp) after the TSS/gene end"')
        st.markdown('**species** (str):')
        st.markdown(
            'Any species of https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/. **Note:** if gene ID is use, species is useless')

        st.divider()
        st.markdown('**Return**')
        st.markdown('(str)')
        st.markdown('FASTA output with gene name, species, accession reference of chromosome used, '
                    'regulatory region (prom/term), chromosomal coordinate of TSS or gene end, distance of TSS or gene '
                    'end on the sequence and sequence extracted ')

        st.divider()
        st.markdown('**Example**')
        exemple = '''gene_ids = ['4843', 'PRKN']
upstream = 1000
downstream = 500

results = []
for gene_id in gene_ids:
    if gene_id.isdigit():
        prom_term = 'Promoter'
        results.append(NCBIdna(gene_id, prom_term, upstream, downstream).find_sequences())
    else:
        prom_term = 'Terminator'
        species = 'Macaca mulatta'
        results.append(NCBIdna(gene_id, prom_term, upstream, downstream,species).find_sequences())

print(results)'''
        st.code(exemple)
        st.markdown("**Result**")
        result = '''>NOS2 | Homo sapiens | NC_000017.11 | Promoter | Gene end (on chromosome): 27756765 | Gene end (on sequence): 1000
TTTGAGAGGCTGAAGTGGGCAGATCACTTGAGCTTCAGAGTTCGAGACCAGCATGGACAACATGGTGAAACCCAGTCTCTACCAAAAACACAAAAATATTAGCTGGGTGTGGTGGTGCATGCCTGTAGTCCCAGCTACTCAGGAGGCTGAGGTGGGAGGATCGCTTGAGCCTGGGAGGCAGAAGTTGCAATGAGCAGAGATCGTGCCACTCCGCTCCAGTCTTGGTGACAGAATGAGACTCCATCTCAAAAATAAATAAATAAATAAAATAAATGAAATGAAATTATAAGAAATTACCACTTTTTCATGTAAGAAGTGATCATTTCCATTATAAGGGAAGGAATTTAATCCTACCTGCCATTCCACCAAAGCTTACCTAGTGCTAAAGGATGAGGTGTTAGTAAGACCAACATCTCAGAGGCCTCTCTGTGCCAATAGCCTTCCTTCCTTTCCCTTCCAAAAACCTCAAGTGACTAGTTCAGAGGCCTGTCTGGAATAATGGCATCATCTAATATCACTGGCCTTCTGGAACCTGGGCATTTTCCAGTGTGTTCCATACTGTCAATATTCCCCCAGCTTCCTGGACTCCTGTCACAAGCTGGAAAAGTGAGAGGATGGACAGGGATTAACCAGAGAGCTCCCTGCTGAGGAAAAAATCTCCCAGATGCTGAAAGTGAGGCCATGTGGCTTGGCCAAATAAAACCTGGCTCCGTGGTGCCTCTATCTTAGCAGCCACCCTGCTGATGAACTGCCACCTTGGACTTGGGACCAGAAAGAGGTGGGTTGGGTGAAGAGGCACCACACAGAGTGATGTAACAGCAAGATCAGGTCACCCACAGGCCCTGGCAGTCACAGTCATAAATTAGCTAACTGTACACAAGCTGGGGACACTCCCTTTGGAAACCAAAAAAAAAAAAAAAAAAAAGAGACCTTTATGCAAAAACAACTCTCTGGATGGCATGGGGTGAGTATAAATACTTCTTGGCTGCCAGTGTGTTCATAACTTTGTAGCGAGTCGAAAACTGAGGCTCCGGCCGCAGAGAACTCAGCCTCATTCCTGCTTTAAAATCTCTCGGCCACCTTTGATGAGGGGACTGGGCAGTTCTAGACAGTCCCGAAGTTCTCAAGGCACAGGTCTCTTCCTGGTTTGACTGTCCTTACCCCGGGGAGGCAGTGCAGCCAGCTGCAAGGTGAGTTGCCTTCATTTCTGGGGAAGCGGCTGTTTTGAGAGGGTTTTGTTTCTTCCTCTTTGAGAAGGCTCAGAAAATTGTGGGAATTTTCTGCCTACAGAGAGAAGGTGTTGGAAAGTCTAGGTAAAAAAATGCCGACATGAGGTTTGGTCATTTGAACATGGCATCTTGGTCAGATTTCTTTTCTGCAAATATTAGCTGTGGTTGTTACATGACAAGGAAAAACTTTCTGAAAAGTTCAAAATGGGAGTTTGTATCTGCTAAGAATTTTTTTTAACAGAGATTATCTTATCAGTCCTTATAACATGTAA

>PRKN | Macaca mulatta | NC_041757.1 | Terminator | Gene end (on chromosome): 9397206 | Gene end (on sequence): 1000
TGTCTGAAGATCTGCCGACTCCAGCATGGCCCCATGGTGACAACAGACCTGCGACAGGAAGCCCAAAGCTCACATCGAAATGGTAGAGAGATCAAAGTCTCTATCGTAAGGGAAAAAAAGAGAGGTCACAGGCATGAGCCCCGGCAGCCAGTGGCTTGTGTCCACACGGAGTCCAGACCCTGATCAGGCCCTGACTTAGTATCGCTGGCAATCCCACTAAATTGCAGTTCCTTACACTGGCCCGATGCAACAAATCAGGTGGCTCCCTTCTGTCACAGGGACCACACAGTGTTTCCCATCATCCATAGCTTTCTTCCTGATGATGTTTGCATGATTGCGCCTTCCCAGCCTGCATGCTGCATTGGGCTTGCAGTGCCTGAACAAGGTTTGCTCCCAAGAGCTCAGGCACCCTAGGAACCCCTGTTAGACTATTAGACTGTCCAGCTTGGTCTCCTTCCCTTTCTTGGTGGTGGTCTTTTCCCTTTCCAGAATAGAACAGTGATTCTTAAAATAAGTTCGAGCAGGCTGGGCACGGTGGCTCATGCCTGTAATCCCAGCACTTTGGGAGGCTGAGGTGGGTGGATCACGAGGTCAGGAGTTCAAGACCAGCCTGGCCAAGACAGCGAAACCCCATCTCTACTAAAAATACAAAGTTAGCTGGGCGCGGTGGCAGGCGCCTGTAATCCCAGCTTCTCGGGAGGCTGAGGCAGGAGAATCACTTGAACCCGGGGGGCAGAGGTTGCAGTGAGCCGAGATCACGCCACTGAACTCCAGCCTGGGCAACAGAGTAAGACTCTGTCTCAAAAAAAGCAAAATCAGCATCCACTACACATGAAAACGAATCACAGTATTATTTGCACAGGAAGGGTGTAACAAAATATGAATGTATCAAAAAATAGAAATAAAGTCTTTGCAGAAAAAATCTATTTTCTCTGACATGTGTTGAGATTATCTGACAACTCTAAGATTGTACTTAAATTGTCAATAAAAGCATCAAAAGAGCTTCTGAGTCCTGTCTTTGAAGATACATTTTCAATCGATGTGAAGGGTGAGGGAGGGAATCACAGGCTGTGTTACTCCATGGACCTCACTTTTGGCCTGGAGGAATTGCCTGGTGATGAAATGGCCATTTTCATATCGACTCCAAGCTCCAGATTAACCCTGGCTGTTTCTCTTTGTTAAATAAGAGCAGTTTGTTAAGGCTATCAAAATTTGAAAGGGTTCCAAAAAGCAAACCCTCATTCAAAACAATCTAAAAGTTACTTTCTTGCCTTACTAGCCAACTTTTCTAATAAGTGTTTGAAAGAGAAGTTTGCATTGGAAAAGACTATGTATAAGGATTTTTTCTCTCATGTCATTCTAGAAAATTTCAGTGACTTTGCATTGCTGTCAATGTGGGTGTGTTGACGCTGCATTTTCAAAGAACTGCTTTAGTGGCCACGCATTGTTTATTTCACAATGAGCCCATGATGGAGGCATCTGGAATGGAACATTCCGAGAT
'''
        st.code(result)

    with st.expander('Convert gene name to gene ID'):
        st.markdown('Convert gene name to gene ID')
        convert_gene_to_entrez_id = 'NCBIdna.convert_gene_to_entrez_id(gene_name, species)'
        st.code(convert_gene_to_entrez_id)

        st.divider()
        st.markdown('**Parameters**')
        st.markdown('**gene_name** (list):')
        st.markdown('It should be **gene name** or **gene ID**. Support only one by one gene_name.')
        st.markdown('**species** (str):')
        st.markdown(
            'Any species of https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/.')

        st.divider()
        st.markdown('**Return**')
        st.markdown('(str)')
        st.markdown('Gene ID of the gene name for desire species')

        st.divider()
        st.markdown('**Example**')
        exemple = '''gene_list = ['APP', 'PRKN']
species_list = ['Macaca mulatta', 'Human']

results = []
for gene_name in gene_list:
    for species in species_list:
        results.append(NCBIdna.convert_gene_to_entrez_id(gene_name, species))

print(results)'''
        st.code(exemple)
        st.markdown("**Result**")
        result = """'100427716', '351', '710899', '5071'"""
        st.code(result)

    with st.expander('Retrieve gene information'):
        st.markdown(
            'API request for retrieve gene information like accesion number of chromosome, start and stop of the gene')
        get_gene_info = 'NCBIdna.get_gene_info(entrez_id)'
        st.code(get_gene_info)

        st.divider()
        st.markdown('**Parameter**')
        st.markdown('**entrez_id** (str):')
        st.markdown(
            'entrez_id is the ENTREZ_GENE_ID on NCBI. Also it is gene ID with function NCBIdna.convert_gene_to_entrez_id(gene_name, species)')

        st.divider()
        st.markdown('**Return**')
        st.markdown('(dict)')
        st.markdown('Return API request. It is the GeneBank of the gene')

        st.divider()
        st.markdown('**Example**')
        exemple = '''id_list = ['4843', '5071']

results = []
for entrez_id in id_list:
    results.append(NCBIdna.get_gene_info(entrez_id))

print(results)'''
        st.code(exemple)
        st.markdown("**Result**")
        result = """{'uid': '4843', 'name': 'NOS2', 'description': 'nitric oxide synthase 2', 'status': '', 'currentid': '', 'chromosome': '17', 'geneticsource': 'genomic', 'maplocation': '17q11.2', 'otheraliases': 'HEP-NOS, INOS, NOS, NOS2A', 'otherdesignations': 'nitric oxide synthase, inducible|NOS, type II|hepatocyte NOS|inducible NO synthase|inducible NOS|nitric oxide synthase 2, inducible|nitric oxide synthase 2A (inducible, hepatocytes)|nitric oxide synthase, macrophage|peptidyl-cysteine S-nitrosylase NOS2', 'nomenclaturesymbol': 'NOS2', 'nomenclaturename': 'nitric oxide synthase 2', 'nomenclaturestatus': 'Official', 'mim': ['163730'], 'genomicinfo': [{'chrloc': '17', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765, 'exoncount': 27}], 'geneweight': 69810, 'summary': 'Nitric oxide is a reactive free radical which acts as a biologic mediator in several processes, including neurotransmission and antimicrobial and antitumoral activities. This gene encodes a nitric oxide synthase which is expressed in liver and is inducible by a combination of lipopolysaccharide and certain cytokines. Three related pseudogenes are located within the Smith-Magenis syndrome region on chromosome 17. [provided by RefSeq, Jul 2008]', 'chrsort': '17', 'chrstart': 27756765, 'organism': {'scientificname': 'Homo sapiens', 'commonname': 'human', 'taxid': 9606}, 'locationhist': [{'annotationrelease': 'RS_2023_03', 'assemblyaccver': 'GCF_000001405.40', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': 'RS_2023_03', 'assemblyaccver': 'GCF_009914755.1', 'chraccver': 'NC_060941.1', 'chrstart': 28741925, 'chrstop': 28698169}, {'annotationrelease': '110', 'assemblyaccver': 'GCF_000001405.40', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '110', 'assemblyaccver': 'GCF_009914755.1', 'chraccver': 'NC_060941.1', 'chrstart': 28741925, 'chrstop': 28698169}, {'annotationrelease': '109.20211119', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '109.20210514', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '109.20210226', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '109.20201120', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '109.20200815', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '109.20200522', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '109.20200228', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '109.20191205', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '109.20190905', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '109.20190607', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000017.11', 'chrstart': 27800528, 'chrstop': 27756765}, {'annotationrelease': '105.20220307', 'assemblyaccver': 'GCF_000001405.25', 'chraccver': 'NC_000017.10', 'chrstart': 26127554, 'chrstop': 26083791}, {'annotationrelease': '105.20201022', 'assemblyaccver': 'GCF_000001405.25', 'chraccver': 'NC_000017.10', 'chrstart': 26127554, 'chrstop': 26083791}, {'annotationrelease': '105', 'assemblyaccver': 'GCF_000001405.25', 'chraccver': 'NC_000017.10', 'chrstart': 26127554, 'chrstop': 26083791}, {'annotationrelease': '105', 'assemblyaccver': 'GCF_000002125.1', 'chraccver': 'AC_000149.1', 'chrstart': 22334498, 'chrstop': 22290737}, {'annotationrelease': '105', 'assemblyaccver': 'GCF_000306695.2', 'chraccver': 'NC_018928.2', 'chrstart': 26190262, 'chrstop': 26146508}]}, {'uid': '5071', 'name': 'PRKN', 'description': 'parkin RBR E3 ubiquitin protein ligase', 'status': '', 'currentid': '', 'chromosome': '6', 'geneticsource': 'genomic', 'maplocation': '6q26', 'otheraliases': 'AR-JP, LPRS2, PARK2, PDJ', 'otherdesignations': 'E3 ubiquitin-protein ligase parkin|Parkinson disease (autosomal recessive, juvenile) 2, parkin|parkinson juvenile disease protein 2|parkinson protein 2 E3 ubiquitin protein ligase|parkinson protein 2, E3 ubiquitin protein ligase (parkin)', 'nomenclaturesymbol': 'PRKN', 'nomenclaturename': 'parkin RBR E3 ubiquitin protein ligase', 'nomenclaturestatus': 'Official', 'mim': ['602544'], 'genomicinfo': [{'chrloc': '6', 'chraccver': 'NC_000006.12', 'chrstart': 162727765, 'chrstop': 161347416, 'exoncount': 13}], 'geneweight': 77618, 'summary': 'The precise function of this gene is unknown; however, the encoded protein is a component of a multiprotein E3 ubiquitin ligase complex that mediates the targeting of substrate proteins for proteasomal degradation. Mutations in this gene are known to cause Parkinson disease and autosomal recessive juvenile Parkinson disease. Alternative splicing of this gene produces multiple transcript variants encoding distinct isoforms. Additional splice variants of this gene have been described but currently lack transcript support. [provided by RefSeq, Jul 2008]', 'chrsort': '06', 'chrstart': 161347416, 'organism': {'scientificname': 'Homo sapiens', 'commonname': 'human', 'taxid': 9606}, 'locationhist': [{'annotationrelease': 'RS_2023_03', 'assemblyaccver': 'GCF_000001405.40', 'chraccver': 'NC_000006.12', 'chrstart': 162727765, 'chrstop': 161347416}, {'annotationrelease': 'RS_2023_03', 'assemblyaccver': 'GCF_009914755.1', 'chraccver': 'NC_060930.1', 'chrstart': 164090894, 'chrstop': 162705910}, {'annotationrelease': '110', 'assemblyaccver': 'GCF_000001405.40', 'chraccver': 'NC_000006.12', 'chrstart': 162727765, 'chrstop': 161347416}, {'annotationrelease': '110', 'assemblyaccver': 'GCF_009914755.1', 'chraccver': 'NC_060930.1', 'chrstart': 164090894, 'chrstop': 162705910}, {'annotationrelease': '109.20211119', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000006.12', 'chrstart': 162727801, 'chrstop': 161347416}, {'annotationrelease': '109.20210514', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000006.12', 'chrstart': 162727801, 'chrstop': 161347416}, {'annotationrelease': '109.20210226', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000006.12', 'chrstart': 162727801, 'chrstop': 161347416}, {'annotationrelease': '109.20201120', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000006.12', 'chrstart': 162727801, 'chrstop': 161347416}, {'annotationrelease': '109.20200815', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000006.12', 'chrstart': 162727801, 'chrstop': 161347416}, {'annotationrelease': '109.20200522', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000006.12', 'chrstart': 162727801, 'chrstop': 161347416}, {'annotationrelease': '109.20200228', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000006.12', 'chrstart': 162727801, 'chrstop': 161347416}, {'annotationrelease': '109.20191205', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000006.12', 'chrstart': 162727801, 'chrstop': 161347416}, {'annotationrelease': '109.20190905', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000006.12', 'chrstart': 162727801, 'chrstop': 161347416}, {'annotationrelease': '109.20190607', 'assemblyaccver': 'GCF_000001405.39', 'chraccver': 'NC_000006.12', 'chrstart': 162727801, 'chrstop': 161347416}, {'annotationrelease': '105.20220307', 'assemblyaccver': 'GCF_000001405.25', 'chraccver': 'NC_000006.11', 'chrstart': 163148797, 'chrstop': 161768448}, {'annotationrelease': '105.20201022', 'assemblyaccver': 'GCF_000001405.25', 'chraccver': 'NC_000006.11', 'chrstart': 163148797, 'chrstop': 161768448}, {'annotationrelease': '105', 'assemblyaccver': 'GCF_000001405.25', 'chraccver': 'NC_000006.11', 'chrstart': 163148833, 'chrstop': 161768589}, {'annotationrelease': '105', 'assemblyaccver': 'GCF_000002125.1', 'chraccver': 'AC_000138.1', 'chrstart': 160600982, 'chrstop': 159224033}, {'annotationrelease': '105', 'assemblyaccver': 'GCF_000306695.2', 'chraccver': 'NC_018917.2', 'chrstart': 163411070, 'chrstop': 162031090}]}"""
        st.code(result, language="python")

    with st.expander('Extraction of regulatory gene sequences manually'):
        st.markdown('Extraction of regulatory gene sequences manually without gene id. Use chromosome reference and chromosomal coordinates')
        get_dna_sequence = 'NCBIdna.get_dna_sequence(prom_term, upstream, downstream, chraccver, chrstart, chrstop)'
        st.code(get_dna_sequence)

        st.divider()
        st.markdown('**Parameters**')
        st.markdown('**prom_term** ("promoter" or "terminator"):')
        st.markdown(
            'Choose between "promoter" or "terminator". It allows to choose between TSS or gene end')
        st.markdown('**upstream** (integer):')
        st.markdown('Upstream defines the number of base pairs (bp) before the TSS/gene end')
        st.markdown('**downstream** (integer):')
        st.markdown('Downstream defines the number of base pairs (bp) after the TSS/gene end')
        st.markdown('**chraccver** (str):')
        st.markdown(
            'ACCESSION of the chromosome')
        st.markdown('**chrstart** (integer):')
        st.markdown(
            'Coordinate of the TSS on the chromosome')
        st.markdown('**chrstop** (integer):')
        st.markdown(
            'Coordinate of the gene end on the chromosome')

        st.divider()
        st.markdown('**Return**')
        st.markdown('(str)')
        st.markdown('DNA sequence extracted')

        st.divider()
        st.markdown('**Example**')
        exemple = '''prom_term = 'promoter'
upstream = 1000
downstream = 500
chraccver = 'NC_000017.11'
chrstart = 27800528
chrstop = 27756765

results = NCBIdna.get_dna_sequence(prom_term, upstream, downstream, chraccver, chrstart, chrstop)

print(results)'''
        st.code(exemple)
        st.markdown("**Result**")
        result = '''TTTGAGAGGCTGAAGTGGGCAGATCACTTGAGCTTCAGAGTTCGAGACCAGCATGGACAACATGGTGAAACCCAGTCTCTACCAAAAACACAAAAATATTAGCTGGGTGTGGTGGTGCATGCCTGTAGTCCCAGCTACTCAGGAGGCTGAGGTGGGAGGATCGCTTGAGCCTGGGAGGCAGAAGTTGCAATGAGCAGAGATCGTGCCACTCCGCTCCAGTCTTGGTGACAGAATGAGACTCCATCTCAAAAATAAATAAATAAATAAAATAAATGAAATGAAATTATAAGAAATTACCACTTTTTCATGTAAGAAGTGATCATTTCCATTATAAGGGAAGGAATTTAATCCTACCTGCCATTCCACCAAAGCTTACCTAGTGCTAAAGGATGAGGTGTTAGTAAGACCAACATCTCAGAGGCCTCTCTGTGCCAATAGCCTTCCTTCCTTTCCCTTCCAAAAACCTCAAGTGACTAGTTCAGAGGCCTGTCTGGAATAATGGCATCATCTAATATCACTGGCCTTCTGGAACCTGGGCATTTTCCAGTGTGTTCCATACTGTCAATATTCCCCCAGCTTCCTGGACTCCTGTCACAAGCTGGAAAAGTGAGAGGATGGACAGGGATTAACCAGAGAGCTCCCTGCTGAGGAAAAAATCTCCCAGATGCTGAAAGTGAGGCCATGTGGCTTGGCCAAATAAAACCTGGCTCCGTGGTGCCTCTATCTTAGCAGCCACCCTGCTGATGAACTGCCACCTTGGACTTGGGACCAGAAAGAGGTGGGTTGGGTGAAGAGGCACCACACAGAGTGATGTAACAGCAAGATCAGGTCACCCACAGGCCCTGGCAGTCACAGTCATAAATTAGCTAACTGTACACAAGCTGGGGACACTCCCTTTGGAAACCAAAAAAAAAAAAAAAAAAAAGAGACCTTTATGCAAAAACAACTCTCTGGATGGCATGGGGTGAGTATAAATACTTCTTGGCTGCCAGTGTGTTCATAACTTTGTAGCGAGTCGAAAACTGAGGCTCCGGCCGCAGAGAACTCAGCCTCATTCCTGCTTTAAAATCTCTCGGCCACCTTTGATGAGGGGACTGGGCAGTTCTAGACAGTCCCGAAGTTCTCAAGGCACAGGTCTCTTCCTGGTTTGACTGTCCTTACCCCGGGGAGGCAGTGCAGCCAGCTGCAAGGTGAGTTGCCTTCATTTCTGGGGAAGCGGCTGTTTTGAGAGGGTTTTGTTTCTTCCTCTTTGAGAAGGCTCAGAAAATTGTGGGAATTTTCTGCCTACAGAGAGAAGGTGTTGGAAAGTCTAGGTAAAAAAATGCCGACATGAGGTTTGGTCATTTGAACATGGCATCTTGGTCAGATTTCTTTTCTGCAAATATTAGCTGTGGTTGTTACATGACAAGGAAAAACTTTCTGAAAAGTTCAAAATGGGAGTTTGTATCTGCTAAGAATTTTTTTTAACAGAGATTATCTTATCAGTCCTTATAACATGTAA'''
        st.code(result)

    with st.expander('Reverse complement a sequence'):
        st.markdown('As expected...')
        reverse_complement = 'NCBIdna.reverse_complement(dna_sequence)'
        st.code(reverse_complement)

        st.divider()
        st.markdown('**Parameter**')
        st.markdown('**dna_sequence** (str):')
        st.markdown(
            'Use only A, T, G, C')

        st.divider()
        st.markdown('**Return**')
        st.markdown('(str)')
        st.markdown('DNA sequence reverse complemented')

        st.divider()
        st.markdown('**Example**')
        exemple = '''dna_sequence = 'TTTGAGAGGCTGAAGTGGGCAGATCACTTGAGCTTCAGAGTTCGAGACCA'

results = NCBIdna.reverse_complement(dna_sequence)

print(results)'''
        st.code(exemple)
        st.markdown("**Result**")
        result = '''TGGTCTCGAACTCTGAAGCTCAAGTGATCTGCCCACTTCAGCCTCTCAAA'''
        st.code(result)
