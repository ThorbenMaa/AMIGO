import click
import math
import copy
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
random.seed(42)

@click.command()
@click.option(
    "--pdb_file",
    "pdb_file",
    required=True,
    multiple=False,
    type=str,
    default="./structure.pdb",
    help="path to file",
)
@click.option(
    "--leu_scheme",
    "leu_scheme",
    required=True,
    multiple=False,
    type=str,
    default="proS",
    help="proS, proR or both",
)
@click.option(
    "--val_scheme",
    "val_scheme",
    required=True,
    multiple=False,
    type=str,
    default="proS",
    help="proS, proR or both",
)
@click.option(
    "--dist_cut_off",
    "dist_cut_off",
    required=True,
    multiple=False,
    type=float,
    default=8,
    help="distance threshold for NOE",
)
@click.option(
    "--outfile",
    "outfile",
    required=True,
    multiple=False,
    type=str,
    default="test/noe.txt",
    help="outfile path and name",
)
def cli (pdb_file, leu_scheme, val_scheme, dist_cut_off, outfile):

    #print ("hello")
    """
    creation of arrays containing the different cut-off dintancies that will be taken into account
    """
    distArrayMet=[dist_cut_off]
    distArrayIle=[dist_cut_off]
    distArrayLeu=[dist_cut_off]
    distArrayVal=[dist_cut_off]
    distArrayAla=[0]
    distArrayThr=[0]

    d_array= [distArrayMet, distArrayIle, distArrayLeu,distArrayVal, distArrayAla, distArrayThr] #MILVAT cut-off dist
    met=False
    ile=False
    leu=False
    val=False
    ala=False
    thr=False

    if d_array[0]!=[0.0]:
        met=True
    if d_array[1]!=[0.0]:
        ile=True
    if d_array[2]!=[0.0]:
        leu=True
    if d_array[3]!=[0.0]:
        val=True
    if d_array[4]!=[0.0]:
        ala=True
    if d_array[5]!=[0.0]:
        thr=True
    
    # create pdb array
    pdb_array=create_pdb_array( met, ile, leu, val, ala, thr, pdb_file, val_scheme, leu_scheme)

    print ("number of methlys", len(pdb_array))

    # create noe network
    noe_network=calc_distance_between_int_aa_pdb_to_all_aa_pdb (pdb_array, d_array, pdb_array, val_scheme, leu_scheme)

    # create noise by randomly removing noes
    #num_to_remove = len(noe_network)//100*percentage_removed_NOEs
    noe_network_size = 0
    for i in range (0, len(noe_network), 1):
        for j in range(0, len(noe_network[i])-1, 1):
            noe_network_size = noe_network_size+1
    while noe_network_size/len(pdb_array) > 2.1 :
        element_to_remove = random.choice(noe_network)
        noe_network.remove(element_to_remove)
        noe_network_size = 0
        for i in range (0, len(noe_network), 1):
            for j in range(0, len(noe_network[i])-1, 1):
                noe_network_size = noe_network_size+1

    print("noe network size", noe_network_size)

    # write network to txt file
    out_file=open(outfile, "w")
    for i in range (0, len(noe_network), 1):
        for j in range (0, len(noe_network[i])-1, 1):
            out_file.write(str(noe_network[i][j]).replace(",", "").replace("[", ""). replace("]", "").replace("'", "") + "\n")
    out_file.close()


    #print (noe_network[10])








"""
Method that creates array with pdb methyl groups and their positions. Input are booleans which are True if the certain amnio acid type was labelled and False if not. Output is array with amnio acid type, residue number, and xyz positions of each methyl group having an amino acid type that was labelled and therefore should be considered.
"""
def create_pdb_array( met, ile, leu, val, ala, thr, pdb_file, val_scheme, leu_scheme ):
    
    # test whether input is correct
    if leu==True and (leu_scheme not in ["proR", "proS", "both"]):
        raise ValueError('Please provide proper value for leu_scheme. Must be proS, proR or both')
    if val==True and (val_scheme not in ["proR", "proS", "both"]):
        raise ValueError('Please provide proper value for val_scheme. Must be proS, proR or both')

    # gather atom positions               
    pdb_file=open(pdb_file , "r")
    pdb=pdb_file.readlines()
    pdb_file.close()
    aa_pdb_array=[]
    for i in range(0, len(pdb), 1): 
        j=pdb[i].split()
        if j[0]=='ATOM':          
            if j[3]=='MET' and j[2]=='CE'and j[0]=='ATOM' and met==True:
                aa_pdb_array.append([j[3],j[5],j[6],j[7],j[8]])            
            if j[3]=='ALA' and j[2]=='CB'and j[0]=='ATOM' and ala==True:
                aa_pdb_array.append([j[3],j[5],j[6],j[7],j[8]])           
            if j[3]=='VAL' and j[2]=='CG2'and j[0]=='ATOM' and val==True and val_scheme=="proS":
                aa_pdb_array.append([j[3],j[5],j[6],j[7],j[8]])
            if j[3]=='VAL' and j[2]=='CG1'and j[0]=='ATOM' and val==True and val_scheme=="proR":
                aa_pdb_array.append([j[3],j[5],j[6],j[7],j[8]])
            if j[3]=='VAL' and j[2]=='CB'and j[0]=='ATOM' and val==True and val_scheme=="both":
                CG1=pdb[i+1].split()
                CG2=pdb[i+2].split()
                if (CG1[3]=='VAL' and CG1[2]=='CG1'and CG1[0]=='ATOM'
                    and CG2[3]=='VAL' and CG2[2]=='CG2'and CG2[0]=='ATOM'):
                    
                    CG1_pos_temp=[CG1[6],CG1[7],CG1[8]]
                    CG2_pos_temp=[CG2[6],CG2[7],CG2[8]]
                    #pseudo_atom_pos_temp = [(float(CG1_pos_temp[0]) + float(CG2_pos_temp[0]))/2,
                    #                        (float(CG1_pos_temp[1]) + float(CG2_pos_temp[1]))/2,
                    #                        (float(CG1_pos_temp[2]) + float(CG2_pos_temp[2]))/2,
                    #                         ]
                    #aa_pdb_array.append([j[3], j[5], pseudo_atom_pos_temp[0], pseudo_atom_pos_temp[1], pseudo_atom_pos_temp[2]])
                    aa_pdb_array.append([j[3], j[5], CG1_pos_temp[0], CG1_pos_temp[1], CG1_pos_temp[2], CG2_pos_temp[0], CG2_pos_temp[1], CG2_pos_temp[2]])
                else:
                    raise ValueError('no correct pdb format. VAL atoms CB, CG1, and CG2 must be in subsequent rows')
                
            if j[3]=='LEU' and j[2]=='CD2'and j[0]=='ATOM' and leu==True and leu_scheme=="proS":
                aa_pdb_array.append([j[3],j[5],j[6],j[7],j[8]]) 
            if j[3]=='LEU' and j[2]=='CD1'and j[0]=='ATOM' and leu==True and leu_scheme=="proR":
                aa_pdb_array.append([j[3],j[5],j[6],j[7],j[8]]) 
            if j[3]=='LEU' and j[2]=='CG'and j[0]=='ATOM' and leu==True and leu_scheme=="both":
                CD1=pdb[i+1].split()
                CD2=pdb[i+2].split()
                if (CD1[3]=='LEU' and CD1[2]=='CD1'and CD1[0]=='ATOM'
                    and CD2[3]=='LEU' and CD2[2]=='CD2'and CD2[0]=='ATOM'):
                    
                    CD1_pos_temp=[CD1[6],CD1[7],CD1[8]]
                    CD2_pos_temp=[CD2[6],CD2[7],CD2[8]]
                    #pseudo_atom_pos_temp = [(float(CD1_pos_temp[0]) + float(CD2_pos_temp[0]))/2,
                    #                        (float(CD1_pos_temp[1]) + float(CD2_pos_temp[1]))/2,
                    #                        (float(CD1_pos_temp[2]) + float(CD2_pos_temp[2]))/2,
                    #                         ]
                    #aa_pdb_array.append([j[3], j[5], pseudo_atom_pos_temp[0], pseudo_atom_pos_temp[1], pseudo_atom_pos_temp[2]])
                    aa_pdb_array.append([j[3], j[5], CD1_pos_temp[0], CD1_pos_temp[1], CD1_pos_temp[2], CD2_pos_temp[0], CD2_pos_temp[1], CD2_pos_temp[2]])

            if j[3]=='ILE' and j[2]=='CD1'and j[0]=='ATOM' and ile==True:
                aa_pdb_array.append([j[3],j[5],j[6],j[7],j[8]])
            if j[3]=='THR' and j[2]=='CG2'and j[0]=='ATOM' and thr==True:
                aa_pdb_array.append([j[3],j[5],j[6],j[7],j[8]])
            
    return aa_pdb_array


"""
Method to create pdb based building blocks. Inputs are the array with all pdb methyl groups and their positions (pdb_array), the array with cut-off disatancies of different amino acid types (d_array; consists of 6 arrays with all cut-off distances for MILVAT, respectively)
and the pdb methyl groups of interest that corresponding building blocks should be created (int_aa_pdb). Output is an array of pdb based building blocks of the pdb methyl groups of interest.

Note: in case of multimers, the same methyl group can pop up multiple times here.
"""
def calc_distance_between_int_aa_pdb_to_all_aa_pdb (pdb_array, d_array, int_aa_pdb, val_scheme, leu_scheme):
    if int_aa_pdb==[]:
        return []
    
    # test whether input is correct
    if leu_scheme not in ["proR", "proS", "both"]:
        raise ValueError('Please provide proper value for leu_scheme. Must be proS, proR or both')
    if val_scheme not in ["proR", "proS", "both"]:
        raise ValueError('Please provide proper value for val_scheme. Must be proS, proR or both')
    
    max_dist_array=[]
    
    #figure out amino acid type of pdb methyl group of interest. Note: Even though int_aa_pdb can contain several methyl groups of interest, they all have always the same amino acid type. That's why you can determine the amino acid type by only taking the first methyl group into account.
    aa_type=int_aa_pdb[0][0]
    if aa_type=='MET':
         max_dist_array=d_array[0]
    if aa_type=='ILE':
         max_dist_array=d_array[1]
    if aa_type=='LEU':
         max_dist_array=d_array[2]
    if aa_type=='VAL':
         max_dist_array=d_array[3]
    if aa_type=='ALA':
         max_dist_array=d_array[4]
    if aa_type=='THR':
         max_dist_array=d_array[5]
    
    int_aa=[]   
    
    
    #calculate distancies to surrounding pdb methyl groups based on the provided cut-off distancies and creation of array with building blocks corresponding to methyl groups of interest.
    for i in range (0, len(int_aa_pdb), 1):
        for k in range(0, len( max_dist_array), 1):
            int_aa.append([])
            for j in range (0, len(pdb_array), 1):
                if pdb_array[j][1]!=int_aa_pdb[i][1]: #calculates proR dists in case one or both are val or leu
                    x=float(pdb_array[j][2])-float(int_aa_pdb[i][2])
                    y=float(pdb_array[j][3])-float(int_aa_pdb[i][3])
                    z=float(pdb_array[j][4])-float(int_aa_pdb[i][4])
                    dist=math.sqrt(x**2+y**2+z**2)

                    # calculate additional distances in case of VAL and LEU and "both" as labeling scheme
                    if (aa_type in ["LEU", "VAL"] and pdb_array[j][0] in ["LEU", "VAL"]) and leu_scheme == "both" and val_scheme == "both": #if both are leu or val
                        # dist int_aa proR and pdb_aa proS
                        x_alt1=float(pdb_array[j][5])-float(int_aa_pdb[i][2])
                        y_alt1=float(pdb_array[j][6])-float(int_aa_pdb[i][3])
                        z_alt1=float(pdb_array[j][7])-float(int_aa_pdb[i][4])
                        dist_alt1=math.sqrt(x_alt1**2+y_alt1**2+z_alt1**2)

                        # dist int_aa proS and pdb_aa proS
                        x_alt2=float(pdb_array[j][5])-float(int_aa_pdb[i][5])
                        y_alt2=float(pdb_array[j][6])-float(int_aa_pdb[i][6])
                        z_alt2=float(pdb_array[j][7])-float(int_aa_pdb[i][7])
                        dist_alt2=math.sqrt(x_alt2**2+y_alt2**2+z_alt2**2)

                        # dist int_aa proS and pdb_aa proR
                        x_alt3=float(pdb_array[j][2])-float(int_aa_pdb[i][5])
                        y_alt3=float(pdb_array[j][3])-float(int_aa_pdb[i][6])
                        z_alt3=float(pdb_array[j][4])-float(int_aa_pdb[i][7])
                        dist_alt3=math.sqrt(x_alt3**2+y_alt3**2+z_alt3**2)

                        dists=np.array([dist, dist_alt1, dist_alt2, dist_alt3])
                        dist=np.min(dists)

                    elif (aa_type in ["LEU", "VAL"] and pdb_array[j][0] not in ["LEU", "VAL"]) and leu_scheme == "both" and val_scheme == "both": #if int_aa is leu or val
                        # dist int_aa proS and pdb_aa normal
                        x_alt3=float(pdb_array[j][2])-float(int_aa_pdb[i][5])
                        y_alt3=float(pdb_array[j][3])-float(int_aa_pdb[i][6])
                        z_alt3=float(pdb_array[j][4])-float(int_aa_pdb[i][7])
                        dist_alt3=math.sqrt(x_alt3**2+y_alt3**2+z_alt3**2)

                        dists=np.array([dist, dist_alt3])
                        dist=np.min(dists)

                    elif (aa_type not in ["LEU", "VAL"] and pdb_array[j][0] in ["LEU", "VAL"]) and leu_scheme == "both" and val_scheme == "both": #if pdb methyl is leu or val
                        # dist int_aa normal and pdb_aa proS
                        x_alt3=float(pdb_array[j][5])-float(int_aa_pdb[i][2])
                        y_alt3=float(pdb_array[j][6])-float(int_aa_pdb[i][3])
                        z_alt3=float(pdb_array[j][7])-float(int_aa_pdb[i][4])
                        dist_alt3=math.sqrt(x_alt3**2+y_alt3**2+z_alt3**2)
                        
                        dists=np.array([dist, dist_alt3])
                        dist=np.min(dists)

                    if dist<= max_dist_array[k]:
                        int_aa[-1].append([int_aa_pdb[i][0], int_aa_pdb[i][1], pdb_array[j][0], pdb_array[j][1]])
            int_aa[-1].append(max_dist_array[k])
    return int_aa
                
if __name__ == "__main__":
    cli()