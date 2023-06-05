# Responsive Elements Finder 🧬🔍

## About

Python script to quickly extract a promoter region with the NCBI API and search for the presence of transcription factor responsive elements.

## Description

First of all, I think it already exists. But even if I looked hard enough, I couldn't find an application or website that really does it the way I want. Of course, you can do a ctrl+F but it's always the same. You have to look for all the shapes in all the possible ways. Of course, there are applications that do it (SerialCloner), but once again, there's something missing. 

When you have an idea, you want it to happen fast. Searching for a promoter sequence can be tedious, and database websites aren't necessarily designed for novices. And that's where my little script comes in. It extracts the desired gene promoter region. You can choose the distance upstream and downstream. It is capable of knowing the direction of the gene and proceeding to reverse complement.

All you have to do is search for your responsive elements. No need to ctrl+F, it can do it. It also accepts IUPAC code and finds all possible shapes in all directions, reverse, complement, reverse complement. And last but not least, it gives you the coordinates of responsive element from the transcription initiation site.

## OS supported

- Win 7 : ❌ (use [REF StreamLit](https://responsive-elements-finder2.streamlit.app/))
- Win 10/11 : ✅
- Linux/MacOS : Not tested
- WebUI ✨ : [REF StreamLit](https://responsive-elements-finder2.streamlit.app/)

## Functions
### Promoter Finder (requires internet connection)
- Extract mutliple promoter regions using ENTREZ_GENE_ID in FASTA format (NEW ✨: v3.X)

### Responsive Elements Finder (no internet connection required)
- Support multiple promoter regions in FASTA format (NEW ✨: v3.X)
- Find transcription factor responsive elements
- Support IUPAC code for responsive elements
- Calculation of the distance of the found sequence to the transcription start site (TSS)
- Percentage of homology between found sequences and responsive elements (NEW ✨ Fixed: 3.X)
- Find mismatches sequence
- Export results to excel
- Graph of sequence positions found on the promoter (NEW ✨: WebUI only)

## WebUI
Yes, there's a WebUI version too 😊 nothing's too good for you 😊

- Use [REF StreamLit](https://responsive-elements-finder2.streamlit.app/)

Graph of sequence positions found on the promoter
![graph_webui](https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/Graph%20WebUI.png)

## Windows version

- Open latest release page - [Releases](https://github.com/Jumitti/Responsive-Elements-Finder/releases/latest)
- Download ``Responsive.Elements.Finder-v?.exe``
- Run app
- Enjoy 😊

Note: Python packages are not required

## Installation/Requirements for Python version
Made for and on Windows. Maybe works on Linux and MacOS (please install python packages)

- Install ``python-3.10.11`` (or above) https://www.python.org/downloads/
- Install python packages with cmd (Windows: run ``python_packages_(windows).bat``):
    ```shell
    pip install pandas
    pip install pillow
    pip install pyperclip
    pip install openpyxl
    pip install requests
    pip install tabulate
    pip install tk
    ```
- Run ``Responsive Element Finder.v?.py``
- Enjoy ☺

## Promoter Finder

I use the NCBI API (https://www.ncbi.nlm.nih.gov/home/develop/api/). For more information on this part, please refer to ``Promoter_finder_HELP.pdf`` (or in the app top left corner ``How to use``)

## Responsive Elements Finder

### Promoter region

- Put your ``promoter region`` or for **multiple promoter regions** use **FASTA format** like below (Required: all sequences must have the TSS at the same distance, otherwise you assume the inconsistency of the positions of found sequences):
    ```shell
    > Gene Name 1
    ATGCCCGGAGATTTCCGATCGCCGCGAATTTTGGCGCGAGAG
    > Gene Name 2
    TGCCGGTGCTGCCCGTAAATGTAAAATGCGCGATGCGTATGC
    ```

- Responsive Elements (RE) allows [IUPAC nucleotides code](https://www.bioinformatics.org/sms/iupac.html)

- Transcription Start Site (TSS): Set the distance of the TSS from the beginning of the pasted sequence or 'Upstream' used in Promoter Finder. Otherwise leave 0.

- Threshold: excludes sequences with homology < threshold

## Screenshot

![screenshot](https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/Responsive%20Elements%20Finder-v2.png)

## WARNING

``Promoter_finder_HELP.pdf`` and ``REF.png`` must be in the same folder as Responsive-Elements-Finder.py.

## Error

In the PDF there is a small error in relation to NC_XXXXXX.XX. This is the accesion code of the chromosome used. But if you want the chromosome number, it's directly on the gene page.
