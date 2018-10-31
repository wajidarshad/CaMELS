# Code to run on a local machine
# CaMELS: In silico prediction of calmodulin binding proteins and their binding sites

Due to Ca2+‐dependent binding and the sequence diversity of Calmodulin (CaM) binding proteins, identifying CaM interactions and binding sites in the wet‐lab is tedious and costly. Therefore, computational methods for this purpose are crucial to the design of such wet‐lab experiments. We present an algorithm suite called CaMELS (CalModulin intEraction Learning System) for predicting proteins that interact with CaM as well as their binding sites using sequence information alone. CaMELS offers state of the art accuracy for both CaM interaction and binding site prediction and can aid biologists in studying CaM binding proteins. For CaM interaction prediction, CaMELS uses protein sequence features coupled with a large‐margin classifier. CaMELS models the binding site prediction problem using multiple instance machine learning with a custom optimization algorithm which allows more effective learning over imprecisely annotated CaM‐binding sites during training. CaMELS has been extensively benchmarked using a variety of data sets, mutagenic studies, proteome‐wide Gene Ontology enrichment analyses and protein structures. Our experiments indicate that CaMELS outperforms simple motif‐based search and other existing methods for interaction and binding site prediction. We have also found that the whole sequence of a protein, rather than just its binding site, is important for predicting its interaction with CaM. Using the machine learning model in CaMELS, we have identified important features of protein sequences for CaM interaction prediction as well as characteristic amino acid sub‐sequences and their relative position for identifying CaM binding sites. Python code for training and evaluating CaMELS together with a webserver implementation is available at the URL: http://faculty.pieas.edu.pk/fayyaz/software.html#camels.
# Dependencies
1- Biopython (https://biopython.org/wiki/Download); 
2- Propy-1.0 (https://code.google.com/archive/p/protpy/downloads)

# Usage
run 'main_runner.py' by changing 'example.fasta' at line no. 12 with path of your own fasta file
# Citation details
W. A. Abbasi, A. Asif, S. Andleeb, and F. U. A. A. Minhas, “CaMELS: In silico prediction of calmodulin binding proteins and their binding sites,” Proteins, vol. 85, no. 9, pp. 1724–1740, Sep. 2017.
# Webserver
CaMELS webserver is vailable at https://camels.pythonanywhere.com/.
