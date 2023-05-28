# Responsive Elements Finder

## Topics

Python script to quickly extract a promoter region with the NCBI API and search for the presence of transcription factor responsive elements.

## Description

First of all, I think it already exists. But even if I looked hard enough, I couldn't find an application or website that really does it the way I want. Of course, you can do a ctrl+F but it's always the same. You have to look for all the shapes in all the possible ways. Of course, there are applications that do it (SerialCloner), but once again, there's something missing. 

When you have an idea, you want it to happen fast. Searching for a promoter sequence can be tedious, and database websites aren't necessarily designed for novices. And that's where my little script comes in. It extracts the desired gene promoter region. You can choose the distance upstream and downstream. It is capable of knowing the direction of the gene and proceeding to reverse complement.

All you have to do is search for your responsive elements. No need to ctrl+F, it can do it. It also accepts IUPAC code and finds all possible shapes in all directions, reverse, complement, reverse complement. And last but not least, it gives you the coordinates of responsive element from the transcription initiation site.

## Installation WINDOWS

Advice: if you use the Responsive-Element-Finder.zip release, you'll find everything you need inside. I.e. installer for python 32-bit and 64-bit. The .bat or .exe script for  python packages. The PDF Helper for the Promoter Finder.

- Install ``python-3.10.11`` https://www.python.org/downloads/
- Install python packages with ``python_packages_(windows).EXE`` or ``python_packages_(windows).bat``. You can also install with cmd.exe:
    ```shell
    pip install requests
    pip install pillow
    pip install pyperclip
    pip install tabulate
    ```
- Run ``Responsive-Element-Finder.exe`` or ``Responsive-Element-Finder.py``
- Enjoy :)

## Promoter Finder

I use the NCBI API. For more information on this part, please refer to ``Promoter_finder_HELP.pdf`` (or in the app top left corner ``How to use``)

![universal remotes](https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/HELP_pdf.png)

## Responsive Elements Finder

All you have to do is stick on your sequence.

For Responsive Elements (RE), you can use the IUPAC code or just ATGC.

The Transcription Initiation Site (TIS) allows you to calculate the correct coordinates of the REs found. Set the distance of the TIS from the beginning of the pasted sequence or 'Upstream' used in Promoter Finder. Otherwise leave 0.

## Enjoy 😊

![universal remotes](https://raw.githubusercontent.com/Jumitti/Responsive-Elements-Finder/main/img/Responsive%20Elements%20Finder.png)

## WARNING

``Promoter_finder_HELP.pdf`` and ``REF.png`` must be in the same folder as Responsive-Elements-Finder.py and Responsive-Elements-Finder.exe

## Error

In the PDF there is a small error in relation to NC_XXXXXX.XX. This is the accesion code of the chromosome used. But if you want the chromosome number, it's directly on the gene page.
