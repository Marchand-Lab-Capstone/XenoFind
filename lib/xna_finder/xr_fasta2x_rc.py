########################################################################
########################################################################
"""
xm_fasta2x_rc.py 


By: H. Kawabe, N. Kaplan, J. Suamabat, J. A. Marchand

Updated: 5/20/24
"""
########################################################################
########################################################################

import sys
import numpy as np
from Bio.Seq import Seq
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from xf_params import *


#Input fasta 
input_fasta=sys.argv[1]

#Output xfasta
output_fasta=sys.argv[2]

print('XenoFind [STATUS] - Generating forward and reverse xFasta files')
'''
note, consider adding these to xr_tools for reuse
'''
def fetch_xna_pos(xm_header):
    pos=xm_header[xm_header.find('XPOS[')+5:-1].split('-')
    xpos = [x.split(':') for x in pos]
    return xpos
    
def xna_base_rc(xna_base, xna_bp): 
	for x in xna_bp: 
		if xna_base in x:
			if len(x)>1:
			    xx=list(x)
			    xx.remove(xna_base)
			    xrc=xx[0]
			else:
			   xrc=False
	return xrc

#Optional confounding pair override
if len(sys.argv)==5: 
    conf_bases = sys.argv[3]
    xna_bases = sys.argv[4]
    confounding_pairs = [xna_bases[0]+conf_bases[0], xna_bases[1]+conf_bases[1]]
    print('[xFASTA] - Overriding xm_params.py. Using '+confounding_pairs[0]+' and '+confounding_pairs[1])


#Create writable output file
f = open(output_fasta, "w")
detected_xfasta_header = False 
detected_xna = False
require_rc_fasta = True

with open(input_fasta, "r") as fh:
    for line in fh:
    
        #Get header
        if line[0]=='>':
            header = line
            if 'XPOS[' in header: 
                detected_xfasta_header=True
            
        #Get sequence
        if line[0]!='>':
        
            #Make all upper case
            uline = line.upper()
            
            #Look for non-standard base
            diff = list(set(uline.replace('\n',''))-(set(standard_bases)))
            
            #This sequence contains an XNA
            if len(diff)>0:
            
                #Set up a clean fasta sequence 
                uline_gap = uline #Null position
                
                #Get position of each XNA
                for x in range(0,len(diff)):


                    #Location of every modification 
                    x_loc=[i for i in range(len(uline)) if uline.startswith(diff[x], i)]
                    
                    #Subject XNA base detected 
                    xna_base = diff[x]
                    #Get base for segmentation substiution 
                    substitution_base = xna_base_rc(xna_base, confounding_pairs)


                    #If conjugate base requires a different segmentation model
                    if xna_base_rc(xna_base,xna_segmentation_model_sets)==False: 
                        require_rc_fasta=True


                    #Setup output header - 
                    header_x = header.replace('\n','')+'+XPOS['
                    uline_clean = uline #Standard base set


                    #At every modification position 
                    for xi in x_loc:
                        posx=diff[x]+':'+str(xi)+'-'
                        header_x = header_x+posx
                        uline_clean = uline_clean.replace(diff[x],substitution_base)
                        uline_gap = uline_gap.replace(diff[x],'-')
                        detected_xna = True 

                    #Close and write 
                    header_xw=header_x[:-1]+']\n'
                    f.write(header_xw)
                    f.write(uline_clean)

                    if write_gaps==True: 
                        header_gap = header.replace('\n','')+'+-+_GAP[]\n'
                        uline_gap = uline #Standard base set
                        for xi in x_loc:
                            uline_gap = uline_gap.replace(diff[x],'-')
                        f.write(header_gap)
                        f.write(uline_gap.replace('-',''))
            elif len(diff)==0 and write_no_xna_seq==True: 
                f.write(header)
                f.write(uline)
f.close()



#Generate RC xFASTA file
if require_rc_fasta == True: 
    fr = open(output_fasta.replace('.fa','_rc')+'.fa', "w")
    with open(output_fasta, "r") as fo:
        for line in fo: 
            if line[0]=='>' and 'GAP' not in line:
                header = line
                x_pos_base = fetch_xna_pos(header)
                x_pos_to_rc =[]
                for x in x_pos_base: 
                    x_base = x[0]
                    x_pos = x[1]

                    if xna_base_rc(x_base,xna_segmentation_model_sets)!=False: 
                        xpr = [x_base, x_pos.replace(']','')]
                        x_pos_to_rc.append(xpr) 

            if line[0]!='>' and len(x_pos_to_rc)>0:

                #Setup header
                header_rc = header[0:header.find('+XPOS[')]+'+RC+XPOS['
                #Get sequence
                seq = line
                #Take RC
                seq_rc = str(Seq(seq).reverse_complement())
                #Get Length
                seq_len = len(seq_rc) 

                for x in x_pos_to_rc: 
                    x_base_rc=xna_base_rc(x[0],xna_base_pairs)
                    x_base_rc_sub=xna_base_rc(x_base_rc,confounding_pairs)
                    x_base_rc_pos = seq_len-int(x[1])-1
                    seq_rc = seq_rc[0:x_base_rc_pos]+x_base_rc_sub+seq_rc[x_base_rc_pos+1:]
                    header_rc=header_rc+x_base_rc+':'+str(x_base_rc_pos-1)+'-'
                    

                #Close and write
                header_rc=header_rc[:-1]+']'
                fr.write(header_rc)
                fr.write(seq_rc+'\n')
    fr.close()
else: 
    try: 
        rmv = output_fasta.replace('.fa','_rc')+'.fa'
        os.remove(rmv)
    except: 
        require_rc_fasta = False


if detected_xfasta_header == True: 
    print('XenoFind Status - [Error] Fasta input file already in xfasta format')
    fasta_input_error=True 
else: 
    if detected_xna == False: 
        print('XenoFind Status - [Error] No XNAs (BS/PZ/KX/JV/XY) detected in fasta input sequence.')
        fasta_input_error=True 

