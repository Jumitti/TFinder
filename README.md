# Responsive Elements Finder 🧬🔍

## About

Python script to quickly extract a promoter region with the NCBI API and search for the presence of transcription factor responsive elements.

## Description

First of all, I think it already exists. But even if I looked hard enough, I couldn't find an application or website that really does it the way I want. Of course, you can do a ctrl+F but it's always the same. You have to look for all the shapes in all the possible ways. Of course, there are applications that do it (SerialCloner), but once again, there's something missing. 

When you have an idea, you want it to happen fast. Searching for a promoter sequence can be tedious, and database websites aren't necessarily designed for novices. And that's where my little script comes in. It extracts the desired gene promoter region. You can choose the distance upstream and downstream. It is capable of knowing the direction of the gene and proceeding to reverse complement.

All you have to do is search for your responsive elements. No need to ctrl+F, it can do it. It also accepts IUPAC code and finds all possible shapes in all directions, reverse, complement, reverse complement. And last but not least, it gives you the coordinates of responsive element from the transcription initiation site.

## OS supported

- Win 7 : ❌
- Win 10/11 : ✅
- Linux/MacOS : Not tested

## Functions
### Promoter Finder (requires internet connection)
- Extract the promoter region with NCBI API

### Responsive Elements Finder (no internet connection required)
- Find transcription factor responsive elements
- IUPAC code can be used for responsive elements
- Calculation of the distance of the found sequence to the transcription initiation site

#### ✨NEW✨
- Percentage of homology between the sequence found and the responsive elements (partially, see below)
- Find mismatches sequence (Allows only 25% mismatches compared to responsive element)
- Export results to excel

## Windows version

- Open latest release page - [Releases](https://github.com/Jumitti/Responsive-Elements-Finder/releases/latest)
- Download ``Responsive.Elements.Finder-v?.exe``
- Run app

Note: Python packages are not required

## Installation/Requirements for Python version
Made for and on Windows. Maybe works on Linux and MacOS (please install python packages)

Advice: You can use the source code. I recommend the releases (it's user friendly)

- Install ``python-3.10.11`` (or above) https://www.python.org/downloads/
- Install python packages with ``python_packages_(windows).bat``. You can also install with ``cmd.exe``:
    ```shell
    pip install pandas
    pip install pillow
    pip install pyperclip
    pip install openpyxl
    pip install requests
    pip install tabulate
    pip install tk
    ```
- Run ``Responsive-Element-Finder.exe``
- Enjoy ☺

## Promoter Finder

I use the NCBI API (https://www.ncbi.nlm.nih.gov/home/develop/api/). For more information on this part, please refer to ``Promoter_finder_HELP.pdf`` (or in the app top left corner ``How to use``)

## Responsive Elements Finder

All you have to do is to paste  your sequence.

For Responsive Elements (RE), you can use the IUPAC code or just ATGC.

The Transcription Initiation Site (TIS) allows you to calculate the correct coordinates of the REs found. Set the distance of the TIS from the beginning of the pasted sequence or 'Upstream' used in Promoter Finder. Otherwise leave 0.

Threshold excludes sequences with low homology

Note: the homology percentage is calculated by analyzing the difference between the sequence found and the responsive element. However, my script has a small problem with this step. As soon as it finds a sequence and the homology % is higher than the threshold, it won't check whether it's the best homology % with an other responsive element. That's why, if you refine the threshold, you might find the same sequences with a better threshold. Nothing really dramatic

## Enjoy 😊

![universal remotes](https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/Responsive%20Elements%20Finder-v2.png)

## WARNING

``Promoter_finder_HELP.pdf`` and ``REF.png`` must be in the same folder as Responsive-Elements-Finder.py.

## Error

In the PDF there is a small error in relation to NC_XXXXXX.XX. This is the accesion code of the chromosome used. But if you want the chromosome number, it's directly on the gene page.
