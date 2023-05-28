# Responsive Elements Finder

Python script to quickly extract a promoter region with the NCBI API and search for the presence of transcription factor responsive elements.

First of all, I think it already exists. But even if I looked hard enough, I couldn't find an application or website that really does it the way I want. Of course, you can do a ctrl+F but it's always the same. You have to look for all the shapes in all the possible ways. Of course, there are applications that do it (SerialCloner comes to mind), but once again, there's something missing. 

When you have an idea, you want it to happen fast. Searching for a promoter sequence can be tedious, and database websites aren't necessarily designed for novices. And that's where my little script comes in. It extracts the desired gene promoter region. You can choose the distance upstream and downstream. It is capable of knowing the direction of the gene and proceeding to reverse complement.

All you have to do is search for your responsive elements. No need to ctrl+F, it can do it. It also accepts IUPAC code and finds all possible shapes in all directions, reverse, complement, reverse complement. And last but not least, it gives you the coordinates of responsive element from the transcription initiation site.





