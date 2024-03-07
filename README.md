# Motif Mark

## Goals
The goal of this assignment is to use object-oriented programming to visualize motifs on sequences using Python3 compatible code.

Given a fasta file and a motifs file, the script outputs an image in png file format. Additionally, it is capable of handling ambiguous motifs, up to 10 sequences, and up to 5 different motifs. 

## Environment 
The script is able to be run in the following environment:
- conda create -n my_pycairo pycairo
- conda activate my_pycairo

Please ensure Pycairo has been installed on your computer prior to running.

Link to Pycairo website: https://pycairo.readthedocs.io/en/latest/index.html

## File inputs 
The script has the following argparse input options that are required:
- ```-f``` fasta file
- ```-m``` motif file

## Example bash command
```
./motif-mark-oop.py -f Figure_1.fasta -m Fig_1_motif.txt
```