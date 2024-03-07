#!/usr/bin/env python
import cairo
import argparse 
import re

# argparse function #

def get_args(): #set variables in argparse to get input from command line 
    parser = argparse.ArgumentParser(description="A program to assign variables for motif mark assignment")
    parser.add_argument("-f", "--file1", help="File path to fasta file", required=True)
    parser.add_argument("-m", "--file2", help="File path to motif file", required=True)
    return parser.parse_args()

args=get_args()
fasta_file=args.file1
motif_file=args.file2

# Creating dictionary for degenerate base symbols based on IUPAC #
iupac_dict={"A":"A","C":"C","G":"G","T":"T","U":"[TU]",
    "W":"[ATU]","S":"[CG]","M":"[AC]","K":"[GTU]",
    "R":"[AG]","Y":"[CTU]","B":"[CGTU]","D":"[AGTU]",
    "H":"[ACTU]","V":"[ACG]","N":"[ACGTU]"}

# Creating list for colors of motifs # 
orange=[0.933,0.467,0.200]
blue=[0,0.467,0.733]
teal=[0,0.600,0.533]
pink=[0.933,0.200,0.467]
cyan=[0.200,0.733,0.933]
colors_list=[pink,teal,orange,cyan,blue]
motif_color_dictionary={}

# Get the file names for writing the files #
file_name_list=fasta_file.split("/")
file_name=file_name_list[-1]
file_name=file_name.split(".")[0] #get just the prefix for the file name 
#print(file_name)
svg_file_name=file_name+".svg"
png_file_name=file_name+".png"

# Classes #
class Motif: 
    def __init__(self,name:str,start:int,stop:int,length:int,color:list,y_offset:int):
        self.name=name
        self.start=start
        self.stop=stop
        self.length=length
        self.color=color
        self.y_offset=y_offset
    
    def __repr__(self): 
        return f'Motif({self.name},{self.start},{self.stop},{self.color})'
    
    def draw_motif(self,context):
        '''Given a motif object, this function will draw the motif onto a pre-defined context surface.'''
        context.set_source_rgb(self.color[0],self.color[1],self.color[2]) #grab the different color values 
        context.set_line_width(20) 
        context.move_to(left_margin+self.start,self.y_offset) 
        context.line_to((left_margin+self.stop),self.y_offset)
        context.stroke()

class Gene: 
    def __init__(self,name:str,length:int,y_offset:int):
        self.name=name
        self.length=length
        self.y_offset=y_offset
    
    def __repr__(self): 
        return f'Gene({self.name},{self.length})'
        return rep
    
    def draw_gene(self,context):
        '''Given a gene object, this function will draw the gene onto a pre-defined context surface.'''
        context.set_source_rgba(0,0,0,1)
        context.move_to(left_margin,y_offset-10)
        context.set_font_size(12)
        context.show_text(self.name)
        context.set_line_width(3)
        context.move_to(left_margin,self.y_offset)
        context.line_to((left_margin+self.length),self.y_offset)
        context.stroke()

class Exon: 
    def __init__(self,start:int,end:int,y_offset:int):
        self.start=start
        self.end=end
        self.y_offset=y_offset

    def __repr__(self): 
        return f'Exon({self.start},{self.end})'

    def draw_exon(self,context):
        '''Given an exon object, this function will draw the exon onto a pre-defined context surface.'''
        context.set_source_rgba(0,0,0,1)
        context.set_line_width(20)
        context.move_to(left_margin+self.start,self.y_offset)
        context.line_to((left_margin+self.end),self.y_offset)
        context.stroke()

class GeneGroup:
    def __init__(self,gene:Gene,exon:Exon,motif_list:Motif):
        self.motif_list=motif_list
        self.gene=gene
        self.exon=exon

    def __repr__(self): 
        return f'GeneGroup({self.gene},{self.exon},{self.motif_list})'
    
    def draw_genegroup(self,context):
        #context.set_source_rgb(0,0,0)
        #Have the different classes draw themselves
        self.gene.draw_gene(context)
        self.exon.draw_exon(context)
        for motif in self.motif_list:
            motif.draw_motif(context)

# Helper Functions #

def motif_parser(motif_file):
    '''This function will take a motif file and turn it into a list of motifs as strings.'''
    motif_list=[]
    with open(motif_file,"r") as fh1:
        for line in fh1: 
            line=line.strip().upper() #make all of the motifs upper case for ease of use with iupac dict
            motif_list.append(line)
    return motif_list

def regex_it(sequence):
    regex_str=""
    sequence=sequence.upper()
    for base in sequence:
        if base in iupac_dict: 
            regex_str+=iupac_dict[base]
    return regex_str

def find_motif(sequence:str,known_motif_dict:dict,motif_color_dict:dict):
    '''Give a sequence, a dictionary of motifs in the format {motif:regex motif}, and a dictionary of colors for motifs in format {motif, RGB color value} 
    will return a list of the motifs with their sequence, start position, end position, length, and assigned color.'''
    sequence=sequence.upper() #make everything uppercase to make it easier to search
    motif_list=[]
    motif_num=0
    for key,value in known_motif_dict.items():
        motif_seq=key
        motif_color=motif_color_dict[motif_seq] #assign a color to each motif
        found_motif=re.finditer(value,sequence) #find the start and stop for each regex motif in the sequence
        for start_end in found_motif:
            found_motif=start_end.span()
            motif_start=found_motif[0]
            motif_end=found_motif[1]
            motif_length=motif_end-motif_start
            motif_list.append((motif_seq,motif_start,motif_end,motif_length,motif_color)) #add to list so you can do more than one motif at a time
    return(motif_list)

def oneline_fasta(input_file): 
    '''Takes a fasta file and converts the sequence line to one line, if on multiple lines'''
    with open (input_file,"r") as fh1, open("one_line_fasta_output.fa","w") as fh2:
        line=fh1.readline()
        fh2.write(line) #writes the first header line to the file 
        seq=""
        for fasta in fh1:
            if fasta.startswith(">"):
                fh2.write(seq+"\n")
                fh2.write(fasta)
                seq =""
            else:
                seq+=fasta.strip("\n")
        fh2.write(seq+"\n")
    pass 

def fasta_parser(onelinefasta_file:str):
    '''This takes a fasta file where the sequence line has already been consolidated into one line and returns the header, sequence, and length in a list.'''
    header=onelinefasta_file.readline().strip() #grab the header line 
    sequence=onelinefasta_file.readline().strip() #grab the sequence line 
    length=len(sequence) #get the length of the sequence line 
    if header=="":
        return "" #will help break the while loop later 
    return[header,sequence,length]

def find_exon(sequence):
    '''From a sequence this function will identify the exons (assuming they are denoted by uppercase characters)
        and return the start and end positions'''
    exon_sequence = re.finditer('[A-Z]+',sequence) #find the start stop for all uppercase letters 
    for start_end in exon_sequence:
        exon_start_end=start_end.span() #span gets the start end positions from the iterable object 
        return exon_start_end #return the start and end positions in a tuple 

# Main Function # 
known_motif_dictionary={}
#key is ambiguous motifs, value is regex equivalent
motif_object_list=[]

# Get the motifs and store their regex equivalent in a dictionary #
# Create a dictionary so we can assign colors to the motifs #
motif_list=motif_parser(motif_file) #gives us a list of motifs in the file 
i=0
for motif in motif_list:
    motif_regex=regex_it(motif) #translates to regex form 
    known_motif_dictionary[motif]=motif_regex #adds to dictionary for regex names 
    motif_color_dictionary[motif]=colors_list[i] #adds to dictionary for motif colors
    i+=1

oneline_fasta(fasta_file)
# Get longest line for scaling the context #
with open("one_line_fasta_output.fa","r") as oneline_fasta_file:
    longest_line_length=0
    i=0
    while True:
        read=fasta_parser(oneline_fasta_file)
        if read=="":
            break
        i+=1
        sequence=read[1]
        length=read[2]
        if len(sequence) > longest_line_length:
            longest_line_length=len(sequence)
    number_y=i #Get the number of genes for scaling the context

# Set up the surface and context #
left_margin=25 #make it easy to change left margin 
gene_height=80 #to scale the context
surface=cairo.SVGSurface(svg_file_name,(longest_line_length+100),((number_y*gene_height)+gene_height)) #add gene height again to give extra space for legend 
context=cairo.Context(surface)
context.set_source_rgb(1,1,1) #white
context.paint()
context.set_source_rgb(0, 0, 0) 

i=0
with open("one_line_fasta_output.fa","r") as oneline_fasta_file:
    y_offset=100
    while True:
        read=fasta_parser(oneline_fasta_file)
        if read=="":
            break
        i+=1
        # Get info for the gene objects #
        header=read[0] #get the header line
        name=header.split()[0] #split on tabs so we can get gene name
        name=name[1:] #gene name should be first item 
        sequence=read[1] #get the sequence line 
        length=read[2] #get the length 
        gene_object=Gene(name,length,y_offset)

        # Make the exon objects #
        exon=find_exon(sequence)
        exon_start=exon[0]
        exon_end=exon[1]
        exon_object=Exon(exon_start,exon_end,y_offset)

        # Find the motifs #
        motif_list_class=[]
        motif_list_class=find_motif(sequence,known_motif_dictionary,motif_color_dictionary)
        for object in motif_list_class:
            motif_object=Motif(object[0],object[1],object[2],object[3],object[4],y_offset)
            motif_object_list.append(motif_object) #has to be list so that you can draw more than one motif for each gene
        
        #Make the GeneGroup object #
        #gg_name="gene_group"+str(i)
        gg_object=GeneGroup(gene_object,exon_object,motif_object_list)
        gg_object.draw_genegroup(context)
        y_offset+=50

# Make a legend # 
context.set_font_size(12)
context.select_font_face("Arial",
                     cairo.FONT_SLANT_NORMAL,
                     cairo.FONT_WEIGHT_NORMAL)
context.set_source_rgb(0,0,0)
context.move_to(left_margin,y_offset)
context.show_text("Legend:")

for motif in motif_color_dictionary:
    context.move_to(left_margin,y_offset+10) #move down
    context.rel_move_to(0,10) #leave space for color box 
    position=context.get_current_point() #get the current point
    context.rectangle(*position,10,-10) #draw box at current point
    current_color=motif_color_dictionary[str(motif)] #get color for the box
    context.set_source_rgb(current_color[0],current_color[1],current_color[2]) #set color for the box
    context.fill()

    context.move_to(left_margin+15,y_offset+10) #move right to write the name of the motif
    context.rel_move_to(0,10)
    position=context.get_current_point()
    context.set_source_rgb(0,0,0)
    context.show_text(str(motif))

    y_offset+=10 #offset and loop again

context.move_to(left_margin,y_offset+10)
context.rel_move_to(0,10)
position=context.get_current_point()
context.rectangle(*position,10,-10)
context.set_source_rgb(0,0,0) #color black for exons 
context.fill()

context.move_to(left_margin+15,y_offset+10)
context.rel_move_to(0,10)
position=context.get_current_point()
context.set_source_rgb(0,0,0)
context.show_text("Exon")

# Write the header #
context.set_source_rgb(0,0,0)
context.set_font_size(20)
context.move_to(10,40)
context.select_font_face("Arial",
                     cairo.FONT_SLANT_NORMAL,
                     cairo.FONT_WEIGHT_NORMAL)
context.show_text("Motif Mark")

surface.write_to_png(png_file_name)
surface.finish()     














