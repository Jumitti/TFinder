# TFinder 🧬🔍

## About

Python script to quickly extract a promoter/terminator region with the NCBI API and search for the presence of transcription factor responsive elements from JASPAR and individuals motifs. Comming soon on WebUI ☺

## Description

First of all, I think it already exists. But even if I looked hard enough, I couldn't find an application or website that really does it the way I want. Of course, you can do a ctrl+F but it's always the same. You have to look for all the shapes in all the possible ways. Of course, there are applications that do it (SerialCloner), but once again, there's something missing. 

When you have an idea, you want it to happen fast. Searching for a promoter sequence can be tedious, and database websites aren't necessarily designed for novices. And that's where my little script comes in. It extracts the desired gene promoter/terminator region. You can choose the distance upstream and downstream. It is capable of knowing the direction of the gene and proceeding to reverse complement.

All you have to do is search for your responsive elements. No need to ctrl+F, it can do it. It also accepts IUPAC code and finds all possible shapes in all directions, reverse, complement, reverse complement. You can use also JASPAR_ID of transcription factors and also you can generate a Position Weight MAtrix (PWM) with multiple sequence for a more accuracy. And last but not least, it gives you the coordinates of responsive element from the transcription initiation site.

## Browser compatibility

- Opera GX
- Chrome (also Chromium)
- Safari
- Edge
- Mozilla

## Functions
### Promoter/Terminator Extractor
- Extract mutliple promoter/terminaotr regions using ENTREZ_GENE_ID or NCBI Gene NAme in FASTA format with NCBI API
- Species: Human 🙋🏼‍♂️, Mouse 🖱, Rat 🐀, Drosophila 🦟, Zebrafish 🐟
- Set Upstream and Downstream from Transcription Start Site (TSS) and Gene End

- Mode "Advance": allows to extract for the same gene the promoter and terminator regions for several species

### Binding Site Finder
- Support multiple promoter/terminator regions in FASTA format
- Find transcription factor responsive elements
- Support IUPAC code for responsive elements
- Support PWM transcription factor with JASPAR API
- Generate PWM for personnal responsive elements
- Calculation of the distance of the found sequence to TSS or Gene End
- Relative Score calculation: (score element found - minimum score matrix reference)/(score maximum matrix reference- minimum score matrix reference)
- p-value: 1000000 random sequences of reactive element length are generated based on the proportion of A, T, G, C in the search sequence. p-value=(Nb Rel Score random sequences >= Rel Score of elements founds)/ (Nb random seq generated). p-value is the number of random sequences generated having a relative score greater than or equal to the relative score of the element found divided by the number of random sequences generated
- Export results to excel (.xlsx)
- Export results, sequences, parameters via e-mail
- Graph of sequence positions found on the promoter

![graph_webui](https://raw.githubusercontent.com/Jumitti/TFinder/main/img/promtermoriginal.png)
![graph_webui](https://raw.githubusercontent.com/Jumitti/TFinder/main/img/bsfMS.png)
![graph_webui](https://raw.githubusercontent.com/Jumitti/TFinder/main/img/Graph%20WebUI.png)

