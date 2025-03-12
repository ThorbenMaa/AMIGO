"""
Code written by Thorben Maass, Insitute of Chemistry and Metabolomics, University of Luebeck, Germany.
"""

"""
Load python necessary python libraries.
"""
import click
import math
import copy
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

"""
define input parameters
"""
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
    "--noe_file",
    "noe_file",
    required=True,
    multiple=False,
    type=str,
    default="./noe.txt",
    help="path to file",
)
@click.option(
    "--additional_measure_file",
    "additional_measure_file",
    required=True,
    multiple=False,
    type=str,
    default="./additional_measure.txt",
    help="path to file",
)
@click.option(
    "--additional_xtal_file",
    "additional_xtal_file",
    required=True,
    multiple=False,
    type=str,
    default="./additional_xtal.txt",
    help="path to file",
)
@click.option(
    "--preassignments_file",
    "preassignments_file",
    required=True,
    multiple=False,
    type=str,
    default="./preassignments.txt",
    help="path to file",
)
@click.option(
    "--min_met",
    "min_met",
    required=True,
    multiple=False,
    type=float,
    default=3,
    help="minimal cut-off distance methionine",
)
@click.option(
    "--max_met",
    "max_met",
    required=True,
    multiple=False,
    type=float,
    default=9,
    help="maximal cut-off distance methionine",
)
@click.option(
    "--step_met",
    "step_met",
    required=True,
    multiple=False,
    type=float,
    default=0.2,
    help="steps from minimal to maximal cut-off distance methionine",
)
@click.option(
    "--min_ile",
    "min_ile",
    required=True,
    multiple=False,
    type=float,
    default=3,
    help="minimal cut-off distance ile",
)
@click.option(
    "--max_ile",
    "max_ile",
    required=True,
    multiple=False,
    type=float,
    default=9,
    help="maximal cut-off distance ile",
)
@click.option(
    "--step_ile",
    "step_ile",
    required=True,
    multiple=False,
    type=float,
    default=0.2,
    help="steps from minimal to maximal cut-off distance ile",
)
@click.option(
    "--min_leu",
    "min_leu",
    required=True,
    multiple=False,
    type=float,
    default=3,
    help="minimal cut-off distance leu",
)
@click.option(
    "--max_leu",
    "max_leu",
    required=True,
    multiple=False,
    type=float,
    default=9,
    help="maximal cut-off distance leu",
)
@click.option(
    "--step_leu",
    "step_leu",
    required=True,
    multiple=False,
    type=float,
    default=0.2,
    help="steps from minimal to maximal cut-off distance leu",
)
@click.option(
    "--min_val",
    "min_val",
    required=True,
    multiple=False,
    type=float,
    default=3,
    help="minimal cut-off distance val",
)
@click.option(
    "--max_val",
    "max_val",
    required=True,
    multiple=False,
    type=float,
    default=9,
    help="maximal cut-off distance val",
)
@click.option(
    "--step_val",
    "step_val",
    required=True,
    multiple=False,
    type=float,
    default=0.2,
    help="steps from minimal to maximal cut-off distance val",
)
@click.option(
    "--min_ala",
    "min_ala",
    required=True,
    multiple=False,
    type=float,
    default=3,
    help="minimal cut-off distance ala",
)
@click.option(
    "--max_ala",
    "max_ala",
    required=True,
    multiple=False,
    type=float,
    default=9,
    help="maximal cut-off distance ala",
)
@click.option(
    "--step_ala",
    "step_ala",
    required=True,
    multiple=False,
    type=float,
    default=0.2,
    help="steps from minimal to maximal cut-off distance ala",
)
@click.option(
    "--min_thr",
    "min_thr",
    required=True,
    multiple=False,
    type=float,
    default=3,
    help="minimal cut-off distance thr",
)
@click.option(
    "--max_thr",
    "max_thr",
    required=True,
    multiple=False,
    type=float,
    default=9,
    help="maximal cut-off distance thr",
)
@click.option(
    "--step_thr",
    "step_thr",
    required=True,
    multiple=False,
    type=float,
    default=0.2,
    help="steps from minimal to maximal cut-off distance thr",
)
@click.option(
    "--weight_noe",
    "weight_noe",
    required=True,
    multiple=False,
    type=float,
    default=1,
    help="weight factor noes",
)
@click.option(
    "--weight_additional_resraint1",
    "weight_additional_resraint1",
    required=True,
    multiple=False,
    type=float,
    default=0,
    help="weight factor weight_additional_resraint1",
)
@click.option(
    "--weight_additional_resraint2",
    "weight_additional_resraint2",
    required=True,
    multiple=False,
    type=float,
    default=0,
    help="weight factor weight_additional_resraint2",
)
@click.option(
    "--weight_additional_resraint3",
    "weight_additional_resraint3",
    required=True,
    multiple=False,
    type=float,
    default=0,
    help="weight factor weight_additional_resraint3",
)
@click.option(
    "--weight_additional_resraint4",
    "weight_additional_resraint4",
    required=True,
    multiple=False,
    type=float,
    default=0,
    help="weight factor weight_additional_resraint4",
)
@click.option(
    "--percentage_starting_points",
    "percentage_starting_points",
    required=True,
    multiple=False,
    type=float,
    default=1,
    help="weight factor weight_additional_resraint4",
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
def cli (pdb_file, noe_file, additional_measure_file, additional_xtal_file, preassignments_file,
        min_met, max_met, step_met,
        min_ile, max_ile, step_ile,
        min_leu, max_leu, step_leu,
        min_val, max_val, step_val,
        min_ala, max_ala, step_ala,
        min_thr, max_thr, step_thr,
        weight_noe, weight_additional_resraint1, weight_additional_resraint2, weight_additional_resraint3, weight_additional_resraint4,
        percentage_starting_points,
        leu_scheme, val_scheme):
    """
    This code block is responsible for reading input files.

    noe.txt contains the NOE list. The text file must contain out of four columns with the first one beeing the amino acid type of the donor methyl group resonance,
    the second column the peak ID of the donor methyl group resonance, the third column the amino acid type of the acceptor methyl group resonance and the fourth column the peak ID
    of the acceptor methyl group resonance. Make sure the amino acid types are written in the pdb format (-> MET, ILE, VAL, LEU, ALA, THR) and the peak IDs are integers. 

    IMPORTANT: the NOE list must be sorted by the peak ID of the donor methyl group resonance! Only then the algorithm will be able to disassemble the NMR based graph in its building blocks (see Paper)

    additional_xtal.txt  contains additional parameters such as chemical shifts (not recommended), PCSs, PREs or RDCs that are predicted from a given structural model of your protein of 
    interest. It must contain out of six columns (amino acid type in pdb format, residue number in pdb file, first additional parameter, second additional parameter, third additional parameter, fourth additional parameter).
    In case you have less then four additional parameters for a certain methyl group, you must give it the value 999 instead for the algorithm will recocnize only then that this value should not
    taken into account.

    additional_measure.txt contains additional parameters such as chemical shifts (not recommended), PCSs, PREs or RDCs that were measured. It must contain out of six columns (amino acid type, peak ID, first additional parameter, second additional parameter, third additional parameter, fourth additional parameter).
    In case you have less then four additional parameters for a certain methyl group, you must give it the value 999 instead for the algorithm will recocnize only then that this value should not
    taken into account.

    preassignments.txt contains a list with preassigned methyl group resonances. It must contain four columns (Aminoacid type of methyl group resonance, peak ID, amino acid type of methyl group that the methyl group resonance was assigned to, residue number in pdb file that this methyl group resonance was assigned to).
    Note: yes it is somehow weird to give the methyl group types two times. Just do it ;-).
    """
    gemessen=open(noe_file , "r") 
    arrayGemessen=gemessen.readlines()
    gemessen.close()

    additional_xtal=open(additional_xtal_file, "r")
    chemShift_PRE_PCS_RDC_xtal=additional_xtal.readlines()
    additional_xtal.close()
    chemShift_PRE_PCS_RDC_array_xtal=create_chemShift_PRE_PCS_RDC_xtal(chemShift_PRE_PCS_RDC_xtal)

    additional_measure=open(additional_measure_file, "r")
    chemShift_PRE_PCS_RDC_measure=additional_measure.readlines()
    additional_measure.close()
    chemShift_PRE_PCS_RDC_array_measure=create_chemShift_PRE_PCS_RDC_measure(chemShift_PRE_PCS_RDC_measure)

    pre_ass=open(preassignments_file, "r")
    pre_assignments=pre_ass.readlines()
    pre_ass.close()

    """
    This block of code provides the interactive determination of different parameters.
    """

    weightNOE=float(weight_noe) #(input("tpe weight factor for noes:    "))
    weightChemShift=float(weight_additional_resraint1)#(input("type weight factor for additional restraint1:   "))
    weightPRE=float(weight_additional_resraint2)#(input("type weight factor for additional restraint2:   "))
    weightPCS=float(weight_additional_resraint3)#(input("type weight factor for additional restraint3:   "))
    weightRDC=float(weight_additional_resraint4)#(input("type weight factor for additional restraint4:   "))
    weightArray=[weightNOE, weightChemShift, weightPRE, weightPCS, weightRDC]

    #print("you will be asked for cut-off distances for the different labeled aminoacid types. If you have not labeled a certain amino acid type, use 0 for minimum and maximum and 1 for step.")

    rangeMin=float(min_val)#input("Type minimal cut-off distance in Ångström for Val:   ")
    rangeMax=float(max_val)#input("Type maximum cut-off distance in Ångström for Val:   ")
    step=float(step_val)#input("Type distance step range, with that you want to go from the minimal cut-off distance to the maximal cut-off distance for Val in Ångström. Use only one decimal place:   ")
    distArrayVal=[]  
    rangeMin=float(rangeMin)
    rangeMax=float(rangeMax)+float(step)
    step=int(float(step)*10)

    for i in range (int(rangeMin*10), int(rangeMax*10), step):
        distArrayVal.append(i*0.1)
        
    rangeMin=float(min_ile)#input("Type minimal cut-off distance in Ångström for Ile:   ")
    rangeMax=float(max_ile)#=input("Type maximum cut-off distance in Ångström for Ile:   ")
    step=float(step_ile)#input("Type distance step range, with that you want to go from the minimal cut-off distance to the maximal cut-off distance for Ile in Ångström. Use only one decimal place:   ")
    distArrayIle=[]  
    rangeMin=float(rangeMin)
    rangeMax=float(rangeMax)+float(step)
    step=int(float(step)*10)

    for i in range (int(rangeMin*10), int(rangeMax*10), step):
        distArrayIle.append(i*0.1)
        
    rangeMin=float(min_met)#input("Type minimal cut-off distance in Ångström for Met:   ")
    rangeMax=float(max_met)#input("Type maximum cut-off distance in Ångström for Met:   ")
    step=float(step_met)#input("Type distance step range, with that you want to go from the minimal cut-off distance to the maximal cut-off distance for Met in Ångström. Use only one decimal place:   ")
    distArrayMet=[]  
    rangeMin=float(rangeMin)
    rangeMax=float(rangeMax)+float(step)
    step=int(float(step)*10)

    for i in range (int(rangeMin*10), int(rangeMax*10), step):
        distArrayMet.append(i*0.1)
        
    rangeMin=float(min_leu)#input("Type minimal cut-off distance in Ångström for Leu:   ")
    rangeMax=float(max_leu)#input("Type maximum cut-off distance in Ångström for Leu:   ")
    step=float(step_leu)#input("Type distance step range, with that you want to go from the minimal cut-off distance to the maximal cut-off distance for Leu in Ångström. Use only one decimal place:   ")
    distArrayLeu=[]  
    rangeMin=float(rangeMin)
    rangeMax=float(rangeMax)+float(step)
    step=int(float(step)*10)

    for i in range (int(rangeMin*10), int(rangeMax*10), step):
        distArrayLeu.append(i*0.1)
        
    rangeMin=float(min_ala)#input("Type minimal cut-off distance in Ångström for Ala:   ")
    rangeMax=float(max_ala)#input("Type maximum cut-off distance in Ångström for Ala:   ")
    step=float(step_ala)#input("Type distance step range, with that you want to go from the minimal cut-off distance to the maximal cut-off distance for Ala in Ångström. Use only one decimal place:   ")
    distArrayAla=[]  
    rangeMin=float(rangeMin)
    rangeMax=float(rangeMax)+float(step)
    step=int(float(step)*10)

    for i in range (int(rangeMin*10), int(rangeMax*10), step):
        distArrayAla.append(i*0.1)
        
    rangeMin=float(min_thr)#input("Type minimal cut-off distance in Ångström for Thr:   ")
    rangeMax=float(max_thr)#input("Type maximum cut-off distance in Ångström for Thr:   ")
    step=float(step_thr)#input("Type distance step range, with that you want to go from the minimal cut-off distance to the maximal cut-off distance for Thr in Ångström. Use only one decimal place:   ")
    distArrayThr=[]  
    rangeMin=float(rangeMin)
    rangeMax=float(rangeMax)+float(step)
    step=int(float(step)*10)

    for i in range (int(rangeMin*10), int(rangeMax*10), step):
        distArrayThr.append(i*0.1)
        
    """
    creation of arrays containing the different cut-off dintancies that will be taken into account
    """
    arrayDistMet=distArrayMet
    arrayDistIle=distArrayIle
    arrayDistVal=distArrayVal
    arrayDistAla=distArrayAla
    arrayDistLeu=distArrayLeu
    arrayDistThr=distArrayThr

    d_array= [distArrayMet, distArrayIle, distArrayLeu,distArrayVal, distArrayAla,distArrayThr] #MILVAT cut-off dist
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




    """Calculating expected penalty"""

    ad1_NMR=[]
    ad2_NMR=[]
    ad3_NMR=[]
    ad4_NMR=[]

    ad1_xtal=[]
    ad2_xtal=[]
    ad3_xtal=[]
    ad4_xtal=[]


    for i in range (0, len(chemShift_PRE_PCS_RDC_array_measure), 1):
        if chemShift_PRE_PCS_RDC_array_measure[i][2]!=999:
            ad1_NMR.append(abs(chemShift_PRE_PCS_RDC_array_measure[i][2]))
        if chemShift_PRE_PCS_RDC_array_measure[i][3]!=999:
            ad2_NMR.append(abs(chemShift_PRE_PCS_RDC_array_measure[i][3]))
        if chemShift_PRE_PCS_RDC_array_measure[i][4]!=999:
            ad3_NMR.append(abs(chemShift_PRE_PCS_RDC_array_measure[i][4]))
        if chemShift_PRE_PCS_RDC_array_measure[i][5]!=999:
            ad4_NMR.append(abs(chemShift_PRE_PCS_RDC_array_measure[i][5]))

    for i in range (0, len(chemShift_PRE_PCS_RDC_array_xtal),1): 
        if chemShift_PRE_PCS_RDC_array_xtal[i][2]!=999:
            ad1_xtal.append(abs(chemShift_PRE_PCS_RDC_array_xtal[i][2]))
        if chemShift_PRE_PCS_RDC_array_xtal[i][3]!=999:
            ad2_xtal.append(abs(chemShift_PRE_PCS_RDC_array_xtal[i][3]))
        if chemShift_PRE_PCS_RDC_array_xtal[i][4]!=999:
            ad3_xtal.append(abs(chemShift_PRE_PCS_RDC_array_xtal[i][4]))
        if chemShift_PRE_PCS_RDC_array_xtal[i][5]!=999:
            ad4_xtal.append(abs(chemShift_PRE_PCS_RDC_array_xtal[i][5]))

    if len(ad1_NMR)>0 and len(ad1_xtal)>0:          
        ad1_NMR_np=np.array(ad1_NMR)
        ad1_xtal_np=np.array(ad1_xtal) 
    else:
        ad1_NMR_np=0
        ad1_xtal_np=0   
    if len(ad2_NMR)>0 and len(ad2_xtal)>0:
        ad2_NMR_np=np.array(ad2_NMR) 
        ad2_xtal_np=np.array(ad2_xtal)
    else:
        ad2_NMR_np=0
        ad2_xtal_np=0  
    if len(ad3_NMR)>0 and len(ad3_xtal)>0:
        ad3_NMR_np=np.array(ad3_NMR) 
        ad3_xtal_np=np.array(ad3_xtal) 
    else:
        ad3_NMR_np=0
        ad3_xtal_np=0    
    if len(ad4_NMR)>0 and len(ad4_xtal)>0:
        ad4_NMR_np=np.array(ad4_NMR) 
        ad4_xtal_np=np.array(ad4_xtal)
    else:
        ad4_NMR_np=0
        ad4_xtal_np=0

    

    penalty_array=[10*abs(np.mean(ad1_NMR_np)-np.mean(ad1_xtal_np)), 10*abs(np.mean(ad2_NMR_np)-np.mean(ad2_xtal_np)), 10*abs(np.mean(ad3_NMR_np)-np.mean(ad3_xtal_np)), 10*abs(np.mean(ad4_NMR_np)-np.mean(ad4_xtal_np))]


    """
    creation of lists for output files
    """
    assignment_pair_array_all_solutions=[]
    assignment_score_array=[]
    methyl_walk_array_all_solutions=[]

    cluster_number=0

    """
    creation of array with most unique building blocks of NMR based graph and array with pdb methyl groups and their positions.
    """
    g_nmr=create_G_nmr(arrayGemessen)
    most_unique_aa_array=find_most_unique_aa(g_nmr, percentage_starting_points)
    pdb_array=create_pdb_array( met, ile, leu, val, ala, thr, pdb_file, val_scheme, leu_scheme)


    """
    Assignment is started. In this block of code, the assignment procedure with the (remaining) NMR building blocks as starting points for the methyls walk is done.
    """
    for l in range (0, len(g_nmr), 1):
        print ('progress: ', l//len(g_nmr), ' %')
        cluster_number=cluster_number+1       
        best_assignment_pair=[]
        best_methyl_walk_array=[]
        best_assignment_score=0
        
        #consider each NMR building block part of the most unique NMR building blocks array as potential starting point.
        #Note: the percentage of most unique NMR building blocks that should be considered as potential starting point is interactively provided by the user.
        for g in range (0, len(most_unique_aa_array), 1):
            if most_unique_aa_array[g][1]>1:
                break
                
            assignment_score=0
            assignment_pair=[]
            methyl_walk_array=[]
            
            #loads prevoius assignments and partially reconstructed pdb and NMR based graphs
            if len(assignment_pair_array_all_solutions)>0:
                for h in range (0, len(assignment_pair_array_all_solutions[-1]), 1):   
                    assignment_pair.append(assignment_pair_array_all_solutions[-1][h])
                    methyl_walk_array.append(methyl_walk_array_all_solutions[-1][h])
            
            #looks for pre-assignemnt provided by the user in the array of NMR based most unqiue buliding blocks and in the array of pdb metyhl goups. Results are array with assigned pdb methyl groups (it can have more than one element if the protein of interest is a multimer) and building block of assigned methyl group resonance. 
            most_unique_aa=most_unique_aa_array[g]
            assigned_NMR='not found'
            assigned_pdb='not found'
            for j in range (0, len(pre_assignments), 1):
                if most_unique_aa[0][0][1]==pre_assignments[j].split()[1]:
                    assigned_NMR=pre_assignments[j]
            if assigned_NMR!='not found':
                for j in range (0, len(pdb_array), 1):
                    if pdb_array[j][1]==assigned_NMR.split()[3]:
                        assigned_pdb= [pdb_array[j]]  
             
            #checks if the certain NMR building block that is now regarded as starting point is already assigned. If yes, it won't be considered
            #most_unique_aa=most_unique_aa_array[g]
            for o in range (0, len(assignment_pair), 1):
                if assignment_pair[o][0][0][0][1]==most_unique_aa[0][0][1]:
                    most_unique_aa=[]
                    break
             
            # if most uniq NMR builiding block is pre-assigned
            if assigned_NMR!='not found' and assigned_pdb!='not found':
                
                #creates array with pdb building blocks of assigned pdb methyl group indentified in the previous step.
                int_aa=calc_distance_between_int_aa_pdb_to_all_aa_pdb (pdb_array, d_array, assigned_pdb, val_scheme, leu_scheme)
            
            #if not preassigned
            else:
                
                #identify all pdb methyl groups with the same amino acid type as the node corresponding to the NMR building block that is now regarded as starting point. Result is an array with the corresponding pdb methyl groups and their positions. Alredy assigned ones are excluded.
                int_aa_pdb=find_aa_with_aa_type_of_most_unique_aa_in_pdb (most_unique_aa, pdb_array, assignment_pair)
                
                #creates array with pdb building blocks of pdb methyl groups indentified in the previous step.
                int_aa=calc_distance_between_int_aa_pdb_to_all_aa_pdb (pdb_array, d_array, int_aa_pdb, val_scheme, leu_scheme)
            
            #Tests which building block of the pdb based building blocks matches the NMR based building block regarded as starting point best and assignes the corresponding nodes to each other. Starts/continues reconstruction of NMR and pdb based graphes with the identified building blocks. Starts new methyl walk with the identified building blocks. Calculates comparison score of these building blocks and adds it to the overall assignment score. If no assignment pair can be identified the best_fit method returns 'stop'. Then the second most unique building block will be tested
            if best_fit(int_aa, most_unique_aa,chemShift_PRE_PCS_RDC_array_measure, chemShift_PRE_PCS_RDC_array_xtal, weightArray, penalty_array)!='stop':
                assignment_pair.append(best_fit(int_aa, most_unique_aa,chemShift_PRE_PCS_RDC_array_measure, chemShift_PRE_PCS_RDC_array_xtal, weightArray, penalty_array)[0:2])
                methyl_walk_array.append('start')
                assignment_score=assignment_score+best_fit(int_aa, most_unique_aa,chemShift_PRE_PCS_RDC_array_measure, chemShift_PRE_PCS_RDC_array_xtal, weightArray, penalty_array)[2]
                
                #continue methyl walk until no more assignments can be identified
                p=0
                while len(assignment_pair)>p:
                    p=len(assignment_pair)
                    
                    #creates array with all NMR based building blocks corresponding to nodes in the partially reconstructed NMR based graph that are not already assigned
                    int_aa_umgebung_nmr=[]
                    for j in range (0, len(assignment_pair), 1):
                        int_aa_umgebung_nmr=int_aa_umgebung_nmr+find_neighbors_of_assigned_nmr(assignment_pair, g_nmr, assignment_pair[j])
                    
                    #sorts array from previous step by their uniqueness
                    most_unique_aa_umgebung_array=find_most_unique_aa_in_neighbors (int_aa_umgebung_nmr)
                    
                    #starts with the most unique NMR based building block.
                    for k in range (0, len(most_unique_aa_umgebung_array), 1):
                        
                        #looks for a node within this NMR building block that is already assigned. Memorizes pdb building block corresponding to the node in the partially reconstructed pdb based graph the NMR node is assigned to. Important for the methyl walk
                        pos='not found'
                        most_unique_aa_umgebung=most_unique_aa_umgebung_array[k]
                        for h in range (0, len (assignment_pair), 1):
                            for r in range (0, len(assignment_pair[h][0][0]), 1):                            
                                if assignment_pair[h][0][0][r][3]==most_unique_aa_umgebung_array[k][0][0][1]:
                                    pos=h
                                    break

                        #looks for pre-assignemnt provided by the user in the array of NMR based most unqiue buliding blocks and in the array of pdb metyhl goups. Results are array with assigned pdb methyl groups (it can have more than one element if the protein of interest is a multimer) and building block of assigned methyl group resonance. 
                        #most_unique_aa=most_unique_aa_array[g]
                        assigned_NMR='not found'
                        assigned_pdb='not found'
                        for j in range (0, len(pre_assignments), 1):
                            if most_unique_aa_umgebung[0][0][1]==pre_assignments[j].split()[1]:
                                assigned_NMR=pre_assignments[j]
                        if assigned_NMR!='not found':
                            for j in range (0, len(pdb_array), 1):
                                if pdb_array[j][1]==assigned_NMR.split()[3]:
                                    assigned_pdb= [pdb_array[j]]  
                        
                         # if most uniq NMR builiding block is pre-assigned
                        if assigned_NMR!='not found' and assigned_pdb!='not found':
                            
                            #creates array with pdb building blocks of pdb methyl groups indentified in the previous step.
                            int_aa_umgebung_crystal=calc_distance_between_int_aa_pdb_to_all_aa_pdb (pdb_array, d_array, assigned_pdb, val_scheme, leu_scheme)
                        
                        #if not preassigned
                        else:
                            #creates an array with pdb methyl groups and their positions. The included methyl groups correspond to nodes of the in the previous step identified pdb building block that are not already assigned. Important for the methyl walk.
                            int_aa_umgebung_crystal_pdb=find_neighbors_of_assigned_aa_crystal (assignment_pair, pdb_array, assignment_pair[pos], most_unique_aa_umgebung)
                            
                            #creates array with pdb based building blocks corresponding to the in the previous step idetified pdb methyl groups.
                            int_aa_umgebung_crystal=calc_distance_between_int_aa_pdb_to_all_aa_pdb (pdb_array, d_array, int_aa_umgebung_crystal_pdb, val_scheme, leu_scheme)
                        
                        #Tests which building block of the pdb based building blocks matches the most unique NMR based building block best and assignes the corresponding nodes to each other. Continues reconstruction of NMR and pdb based graphes with the identified building blocks. Continues methyl walk with the identified building blocks. Calculates comparison score of these building blocks and adds it to the overall assignment score. If no assignment pair can be identified the best_fit method returns 'stop'. Then the second most unique building block will be tested
                        if best_fit(int_aa_umgebung_crystal, most_unique_aa_umgebung,chemShift_PRE_PCS_RDC_array_measure, chemShift_PRE_PCS_RDC_array_xtal, weightArray, penalty_array)!='stop':
                            assignment_pair.append(best_fit(int_aa_umgebung_crystal, most_unique_aa_umgebung,chemShift_PRE_PCS_RDC_array_measure, chemShift_PRE_PCS_RDC_array_xtal, weightArray, penalty_array)[0:2])
                            methyl_walk_array.append(methyl_walk_array[pos] + ' => ' + assignment_pair[pos][0][0][0][0]+' '+assignment_pair[pos][0][0][0][1] + ' // ' + assignment_pair[pos][1][0][0]+' '+assignment_pair[pos][1][0][1])
                            assignment_score=assignment_score+best_fit(int_aa_umgebung_crystal, most_unique_aa_umgebung,chemShift_PRE_PCS_RDC_array_measure, chemShift_PRE_PCS_RDC_array_xtal, weightArray, penalty_array)[2]
                            break
            
            #compares overall comparison score of assignments done with the selected starting point to overall comparison score of the starting point with the best overall comparison score. Memorizes the better overall assignment score , the corresponding assignments and partially reconstructed pdb and NMR based graphs and the corresponding methyl walks.
            if assignment_score>best_assignment_score:
                best_assignment_score=assignment_score
                best_assignment_pair=copy.deepcopy(assignment_pair)
                best_methyl_walk_array=copy.deepcopy(methyl_walk_array)
        
        #adds assignments, the corresponding partially reconstructed pdb and NMR based graphs and methyl walk to the final solution. Repeat until no more starting points are identified with that assignments can be achieved.
        methyl_walk_array_all_solutions.append(best_methyl_walk_array)
        assignment_pair_array_all_solutions.append(best_assignment_pair)
        assignment_score_array.append(best_assignment_score)
        if best_assignment_pair==[]:
            print ('progress: ', 100, ' %')
            break

    """
    Here, the final soltion is evaluated. This is important for parameter optimization.
    """
    start=0
    for i in range(0, len(assignment_pair_array_all_solutions)-1, 1):
        evaluate_assignment_pair(assignment_pair_array_all_solutions[i], start)
        start=len(assignment_pair_array_all_solutions[i])
        
    """
    This block of code is dedicated for creation of output files. 
    """
    #Creation of .txt file containing assignments, information about the assignments and methyl walks.
    outfile=open("result.txt", 'w')
    i=0
    j=0
    array_additional_NMR_output=[]
    array_additional_xtal_output=[]
    outfile.write('distance ranges for [M] [I] [L] [V] [A] [T] methyl groups: '+ str(d_array))
    outfile.write('\n')
    outfile.write('\n')
    outfile.write('weighting factors [w_NOE, w_1, w_2, w_3, w_4]            : '+str(weightArray))
    outfile.write('\n')
    outfile.write('\n')
    d_array_met=[]
    d_array_ile=[]
    d_array_leu=[]
    d_array_val=[]
    d_array_ala=[]
    d_array_thr=[]
    if assignment_pair_array_all_solutions[0]==[] and len(assignment_pair_array_all_solutions)==1:
        print ("No assignments identified (that's the reason of the error message)- sorry!")

    total_amount_ass=0
    total_amount_ass_m=0
    total_amount_ass_i=0
    total_amount_ass_l=0
    total_amount_ass_v=0
    total_amount_ass_a=0
    total_amount_ass_t=0

    total_amount_perfect=0
    total_amount_perfect_m=0
    total_amount_perfect_i=0
    total_amount_perfect_l=0
    total_amount_perfect_v=0
    total_amount_perfect_a=0
    total_amount_perfect_t=0

    while i< len(assignment_pair_array_all_solutions)-1:
        
        amount_zug_cluster=0
        amount_plus=0
        amount_minus=0
        amount_o=0
        outfile.write('\ncluster       :      '+ str(i))
        #print('Cluster Score  :      ', assignment_score_array[i])
        max_pos_cluster_score=0
        while j<len(assignment_pair_array_all_solutions[i]):
            max_pos_cluster_score=max_pos_cluster_score+len(assignment_pair_array_all_solutions[i][j][0][0])  
            amount_zug_cluster=amount_zug_cluster+1
            total_amount_ass=total_amount_ass+1
            
            #count perfectly matching building blocks
            if assignment_pair_array_all_solutions[i][j][0][2]=='+':
                amount_plus=amount_plus+1
                total_amount_perfect=total_amount_perfect+1
            elif assignment_pair_array_all_solutions[i][j][0][2]=='-':
                amount_minus=amount_minus+1
            elif assignment_pair_array_all_solutions[i][j][0][2]=='o':
                amount_o=amount_o+1
                
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='MET':
                total_amount_ass_m=total_amount_ass_m+1
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='ILE':
                total_amount_ass_i=total_amount_ass_i+1
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='LEU':
                total_amount_ass_l=total_amount_ass_l+1
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='VAL':
                total_amount_ass_v=total_amount_ass_v+1
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='ALA':
                total_amount_ass_a=total_amount_ass_a+1
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='THR':
                total_amount_ass_t=total_amount_ass_t+1
            
            if assignment_pair_array_all_solutions[i][j][0][2]=='+' and assignment_pair_array_all_solutions[i][j][1][0][0]=='MET':
                total_amount_perfect_m=total_amount_perfect_m+1
            if assignment_pair_array_all_solutions[i][j][0][2]=='+' and assignment_pair_array_all_solutions[i][j][1][0][0]=='ILE':
                total_amount_perfect_i=total_amount_perfect_i+1
            if assignment_pair_array_all_solutions[i][j][0][2]=='+' and assignment_pair_array_all_solutions[i][j][1][0][0]=='LEU':
                total_amount_perfect_l=total_amount_perfect_l+1
            if assignment_pair_array_all_solutions[i][j][0][2]=='+' and assignment_pair_array_all_solutions[i][j][1][0][0]=='VAL':
                total_amount_perfect_v=total_amount_perfect_v+1
            if assignment_pair_array_all_solutions[i][j][0][2]=='+' and assignment_pair_array_all_solutions[i][j][1][0][0]=='ALA':
                total_amount_perfect_a=total_amount_perfect_a+1
            if assignment_pair_array_all_solutions[i][j][0][2]=='+' and assignment_pair_array_all_solutions[i][j][1][0][0]=='THR':
                total_amount_perfect_t=total_amount_perfect_t+1
            
            outfile.write('\n')
            outfile.write('\n')
            outfile.write('\nmethyl walk (NMR//pdb)                   :     ' + str(methyl_walk_array_all_solutions[i][j])+ '\n')
            
            outfile.write('was extened with methyl group resonance  :     ' +str(assignment_pair_array_all_solutions[i][j][0][0][0][0])+' '+str(assignment_pair_array_all_solutions[i][j][0][0][0][1]) + '  with NOEs to resonances ')
    
            
            for k in range (0, len(assignment_pair_array_all_solutions[i][j][0][0]), 1):    
                outfile.write(str(assignment_pair_array_all_solutions[i][j][0][0][k][2])+' ' +str(assignment_pair_array_all_solutions[i][j][0][0][k][3])+' / ')
            outfile.write('\nand thereby assigned to pdb methyl group :     '+str(assignment_pair_array_all_solutions[i][j][1][0][0])+' '+str(assignment_pair_array_all_solutions[i][j][1][0][1]) + ' with expected NOEs to   ')
            
            for k in range (0, len(assignment_pair_array_all_solutions[i][j][1])-1, 1):
                outfile.write( str(assignment_pair_array_all_solutions[i][j][1][k][2]) + ' ' +str(assignment_pair_array_all_solutions[i][j][1][k][3]) + ' / ')
            outfile.write(' applying a cut-off distance of ' +str(assignment_pair_array_all_solutions[i][j][1][-1]) + ' Amgstroem  \n')
            
            #save used cut-off diatnces for different amino acid types for creation of histograms
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='MET':
                d_array_met.append(assignment_pair_array_all_solutions[i][j][1][-1])
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='ILE':
                d_array_ile.append(assignment_pair_array_all_solutions[i][j][1][-1])
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='LEU':
                d_array_leu.append(assignment_pair_array_all_solutions[i][j][1][-1])
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='VAL':
                d_array_val.append(assignment_pair_array_all_solutions[i][j][1][-1])
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='ALA':
                d_array_ala.append(assignment_pair_array_all_solutions[i][j][1][-1])
            if assignment_pair_array_all_solutions[i][j][1][0][0]=='THR':
                d_array_thr.append(assignment_pair_array_all_solutions[i][j][1][-1])
            
            
                
            pos_in_additional_xtal='not found'
            pos_in_additional_measure='not found'
                    
            for f in range (0, len(chemShift_PRE_PCS_RDC_array_xtal),1):
                    if chemShift_PRE_PCS_RDC_array_xtal[f][1]==assignment_pair_array_all_solutions[i][j][1][0][1]:#chemShift_PRE_PCS_RDC_array_xtal[f][0]==assignment_pair_array_all_solutions[i][j][1][0][0] and chemShift_PRE_PCS_RDC_array_xtal[f][1]==assignment_pair_array_all_solutions[i][j][1][0][1]:
                        pos_in_additional_xtal=f
                        
            for f in range (0, len(chemShift_PRE_PCS_RDC_array_measure),1):  		
                    if chemShift_PRE_PCS_RDC_array_measure[f][1]==assignment_pair_array_all_solutions[i][j][0][0][0][1]:#chemShift_PRE_PCS_RDC_array_measure[f][0]==assignment_pair_array_all_solutions[i][j][0][0][0][0] and chemShift_PRE_PCS_RDC_array_measure[f][1]==assignment_pair_array_all_solutions[i][j][0][0][0][1]:
                        pos_in_additional_measure=f
                
            if pos_in_additional_xtal!='not found' and pos_in_additional_measure !='not found':
                array_additional_NMR_output.append(chemShift_PRE_PCS_RDC_array_measure[pos_in_additional_measure])
                array_additional_xtal_output.append(chemShift_PRE_PCS_RDC_array_xtal[pos_in_additional_xtal])             
                outfile.write("additional experimental restraints are [additional restraint 1, 2, 3 and 4]    :"+  str(chemShift_PRE_PCS_RDC_array_measure[pos_in_additional_measure][2:])+ "\n")
                outfile.write("additional theoretical restraints are [additional restraint 1, 2, 3 and 4]     :"+  str(chemShift_PRE_PCS_RDC_array_xtal[pos_in_additional_xtal][2:])+ "\n")      
            
            j=j+1
            
        outfile.write('\namount of assignments in cluster                  : '+str( amount_zug_cluster))
        outfile.write('\namount of perfectly matching building block pairs : '+str( amount_plus))
        outfile.write('\nCluster Score  :      '+ str(assignment_score_array[i]) +' of '+str(weightNOE*max_pos_cluster_score*2))
        outfile.write('\n')
        i=i+1

    outfile.write('\ntotal amount of assignments                             : '+str( total_amount_ass))
    outfile.write('\ntotal amount of perfectly matching building block pairs: '+str( total_amount_perfect))
    outfile.write('\ntotal ratio of perfectly matching building block pairs to toal assignments for MET: '+str( total_amount_perfect_m) + " / "+ str(total_amount_ass_m))
    outfile.write('\ntotal ratio of perfectly matching building block pairs to toal assignments for ILE: '+str( total_amount_perfect_i) + " / "+ str(total_amount_ass_i))
    outfile.write('\ntotal ratio of perfectly matching building block pairs to toal assignments for LEU: '+str( total_amount_perfect_l) + " / "+ str(total_amount_ass_l))
    outfile.write('\ntotal ratio of perfectly matching building block pairs to toal assignments for VAL: '+str( total_amount_perfect_v) + " / "+ str(total_amount_ass_v))
    outfile.write('\ntotal ratio of perfectly matching building block pairs to toal assignments for ALA: '+str( total_amount_perfect_a) + " / "+ str(total_amount_ass_a))
    outfile.write('\ntotal ratio of perfectly matching building block pairs to toal assignments for THR: '+str( total_amount_perfect_t) + " / "+ str(total_amount_ass_t))
    outfile.close()

    #creation of interactive surfaces showing the final pdb and NMR based graphs
    von_gem=[]
    nach_gem=[]
    for i in range(0, len(assignment_pair_array_all_solutions[-2]), 1):
        for j in range(0, len(assignment_pair_array_all_solutions[-2][i][0][0]), 1):
            von_gem.append(str(assignment_pair_array_all_solutions[-2][i][0][0][j][0])+ ' ' +str(assignment_pair_array_all_solutions[-2][i][0][0][j][1]))
            nach_gem.append(str(assignment_pair_array_all_solutions[-2][i][0][0][j][2])+ ' ' +str(assignment_pair_array_all_solutions[-2][i][0][0][j][3]))

    df = pd.DataFrame({ 'from':von_gem, 'to':nach_gem})
    df

    G = nx.from_pandas_edgelist(df, source='from', target='to', create_using=nx.DiGraph())
    H = G.to_directed()

    fig=plt.figure(1, figsize=(40, 40))
    fig.suptitle('G_NMR')

    color_map=[]
    for node in H:
        if node.split()[0]=='MET':
            color_map.append('red')
        elif node.split()[0]=='ILE':
            color_map.append('orange')
        elif node.split()[0]=='LEU':
            color_map.append('yellow')
        elif node.split()[0]=='VAL':
            color_map.append('blue')
        elif node.split()[0]=='ALA':
            color_map.append('cyan')
        elif node.split()[0]=='THR':
            color_map.append('green')
        else:
            color_map.append('grey')


    nx.draw(H, pos=nx.spring_layout(H),  with_labels=True, font_size=8 , node_color=color_map, node_size=800, width=0.2)
            
    von_pdb=[]
    nach_pdb=[]
    for i in range(0, len(assignment_pair_array_all_solutions[-2]), 1):
        for j in range(0, len(assignment_pair_array_all_solutions[-2][i][1])-1, 1):
            von_pdb.append(str(assignment_pair_array_all_solutions[-2][i][1][j][0])+ ' ' +str(assignment_pair_array_all_solutions[-2][i][1][j][1]))
            nach_pdb.append(str(assignment_pair_array_all_solutions[-2][i][1][j][2])+ ' ' +str(assignment_pair_array_all_solutions[-2][i][1][j][3]))

    df = pd.DataFrame({ 'from':von_pdb, 'to':nach_pdb})
    df

    G = nx.from_pandas_edgelist(df, source='from', target='to', create_using=nx.DiGraph())
    H = G.to_directed()
    fig.savefig('G_NMR.svg')

    fig=plt.figure(2, figsize=(40, 40))
    fig.suptitle('G_pdb')

    color_map=[]
    for node in H:
        if node.split()[0]=='MET':
            color_map.append('red')
        elif node.split()[0]=='ILE':
            color_map.append('orange')
        elif node.split()[0]=='LEU':
            color_map.append('yellow')
        elif node.split()[0]=='VAL':
            color_map.append('blue')
        elif node.split()[0]=='ALA':
            color_map.append('cyan')
        elif node.split()[0]=='THR':
            color_map.append('green')
        else:
            color_map.append('grey')

    nx.draw(H, pos=nx.spring_layout(H),  with_labels=True, font_size=8 , node_color=color_map, node_size=800, width=0.2)
    fig.savefig('G_pdb.svg')   

    #creation of interactive surfaces showing predicted and measured additional parameters.
    fig=plt.figure(3, figsize=(8, 8))
    fig.suptitle('additional restraints')
    cs_nmr=[]
    cs_xtal=[]
    pre_nmr=[]
    pre_xtal=[]
    pcs_nmr=[]
    pcs_xtal=[]
    rdc_nmr=[]
    rdc_xtal=[]

    for i in range (0, len(array_additional_xtal_output),1):
        cs_nmr.append(array_additional_NMR_output[i][2])
        pre_nmr.append(array_additional_NMR_output[i][3])
        pcs_nmr.append(array_additional_NMR_output[i][4])
        rdc_nmr.append(array_additional_NMR_output[i][5])
        cs_xtal.append(array_additional_xtal_output[i][2])
        pre_xtal.append(array_additional_xtal_output[i][3])
        pcs_xtal.append(array_additional_xtal_output[i][4])
        rdc_xtal.append(array_additional_xtal_output[i][5])

    i=0
    while i<len (cs_nmr):
        if cs_nmr[i]==999 or cs_xtal[i]==999:
            del cs_nmr[i]
            del cs_xtal[i]
        else:
            i=i+1
    i=0
    while i<len(pre_nmr):
        if pre_nmr[i]==999 or pre_xtal[i]==999:
            del pre_nmr[i]
            del pre_xtal[i]
        else:
            i=i+1
    i=0
    while i<len(pcs_nmr):
        if pcs_nmr[i]==999 or pcs_xtal[i]==999:
            del pcs_nmr[i]
            del pcs_xtal[i]
        else:
            i=i+1
    i=0
    while i<len(rdc_nmr):
        if rdc_nmr[i]==999 or rdc_xtal[i]==999:
            del rdc_nmr[i]
            del rdc_xtal[i]
        else:
            i=i+1

    plt.subplot(2,2,1)
    plt.scatter(cs_nmr, cs_xtal, marker='o', label="additional restraint1")
    plt.xlabel("experimental values")
    plt.ylabel("theoretical values")
    plt.legend()
    plt.subplot(2,2,2)
    plt.scatter(pre_nmr, pre_xtal, marker='o',label="additional restraint2")
    plt.xlabel("experimental values")
    plt.ylabel("theoretical values")
    plt.legend()
    plt.subplot(2,2,3)
    plt.scatter(pcs_nmr, pcs_xtal,marker='o',label="additional restraint3")
    plt.xlabel("experimental values")
    plt.ylabel("theoretical values")
    plt.legend()
    plt.subplot(2,2,4)
    plt.scatter(rdc_nmr, rdc_xtal, marker='o',label="additional restraint4")
    plt.xlabel("experimental values")
    plt.ylabel("theoretical values")
    plt.legend()
    plt.tight_layout()
    fig.savefig('additional_restraints.svg')   

    fig=plt.figure(4, figsize=(10, 10))
    d_met=np.array(d_array_met)
    d_ile=np.array(d_array_ile)
    d_leu=np.array(d_array_leu)
    d_val=np.array(d_array_val)
    d_ala=np.array(d_array_ala)
    d_thr=np.array(d_array_thr)

    fig.suptitle('cut-off diatnces')
    plt.subplot(3,2,1)
    plt.hist(d_array_met, bins=d_array[0], label='MET, \u03bc '+str(round(np.mean(d_met),1))+', standard deviation '+str(round(np.std(d_met),1)))
    plt.legend()
    plt.subplot(3,2,2)
    plt.hist(d_array_ile, bins=d_array[1], label='ILE, \u03bc '+str(round(np.mean(d_ile),1))+', standard deviation '+str(round(np.std(d_ile),1)))
    plt.legend()
    plt.subplot(3,2,3)
    plt.hist(d_array_leu, bins=d_array[2], label='LEU, \u03bc '+str(round(np.mean(d_leu),1))+', standard deviation '+str(round(np.std(d_leu),1)))
    plt.legend()
    plt.subplot(3,2,4)
    plt.hist(d_array_val, bins=d_array[3], label='VAL, \u03bc '+str(round(np.mean(d_val),1))+', standard deviation '+str(round(np.std(d_val),1)))
    plt.legend()
    plt.subplot(3,2,5)
    plt.hist(d_array_ala, bins=d_array[4], label='ALA, \u03bc '+str(round(np.mean(d_ala),1))+', standard deviation '+str(round(np.std(d_ala),1)))
    plt.legend()
    plt.subplot(3,2,6)
    plt.hist(d_array_thr, bins=d_array[5], label='THR, \u03bc '+str(round(np.mean(d_thr),1))+', standard deviation '+str(round(np.std(d_thr),1)))
    plt.legend()
    plt.tight_layout()
    fig.savefig('cut-off-distances.svg')

    plt.show()





"""
Method to sort NMR based building blocks concerning their uniqueness score. Input is array with NMR based building blocks. Output is array with NMR building blocks sorted by their uniqueness score. 
"""
def find_most_unique_aa (g_nmr, percentage_starting_points): 
    
    #counting amount of certain amino acid types occuring as neighbors of active nodes within their building block.
    amount_met=0
    amount_ile=0
    amount_leu=0
    amount_val=0
    amount_ala=0
    amount_thr=0
    amount_total=0
    
    #counting active nodes having a certain amino acid type.
    aa_met=0
    aa_ile=0
    aa_leu=0
    aa_val=0
    aa_ala=0
    aa_thr=0
    amount_aa_total=0
    
    #counts event of an active node having a certain amico acid type as a neighbor.
    for i in range (0, len(g_nmr), 1):
        for j in range (0, len(g_nmr[i]), 1):
            if g_nmr[i][j][2]=='MET':
                amount_met=amount_met+1
                amount_total=amount_total+1
            if g_nmr[i][j][2]=='ILE':
                amount_ile=amount_ile+1
                amount_total=amount_total+1
            if g_nmr[i][j][2]=='LEU':
                amount_leu=amount_leu+1
                amount_total=amount_total+1
            if g_nmr[i][j][2]=='VAL':
                amount_val=amount_val+1
                amount_total=amount_total+1
            if g_nmr[i][j][2]=='ALA':
                amount_ala=amount_ala+1
                amount_total=amount_total+1
            if g_nmr[i][j][2]=='THR':
                amount_thr=amount_thr+1
                amount_total=amount_total+1
                
        #counting amount of active nodes having acertain amnio acid type.  
        if g_nmr[i][0][0]=='MET':
            aa_met=aa_met+1
            amount_aa_total=amount_aa_total+1
        if g_nmr[i][0][0]=='ILE':
            aa_ile=aa_ile+1
            amount_aa_total=amount_aa_total+1
        if g_nmr[i][0][0]=='LEU':
            aa_leu=aa_leu+1
            amount_aa_total=amount_aa_total+1
        if g_nmr[i][0][0]=='VAL':
            aa_val=aa_val+1
            amount_aa_total=amount_aa_total+1
        if g_nmr[i][0][0]=='ALA':
            aa_ala=aa_ala+1
            amount_aa_total=amount_aa_total+1
        if g_nmr[i][0][0]=='THR':
            aa_thr=amount_thr+1
            amount_aa_total=amount_aa_total+1
     
    #probability of an active node having a certain amino acid type.
    p_aa_met=aa_met/amount_aa_total
    p_aa_ile=aa_ile/amount_aa_total
    p_aa_leu=aa_leu/amount_aa_total
    p_aa_val=aa_val/amount_aa_total
    p_aa_ala=aa_ala/amount_aa_total
    p_aa_thr=aa_thr/amount_aa_total    
                          
    #probability of an active node having a certain amino acid type as a neighbor.
    p_met=amount_met/amount_total
    p_ile=amount_ile/amount_total
    p_leu=amount_leu/amount_total
    p_val=amount_val/amount_total
    p_ala=amount_ala/amount_total
    p_thr=amount_thr/amount_total
    
    #calculation of uniqueness score for building blocks. Result is p_array. The i-th element in p_array corresponds to the i-th building block in g_nmr (NMR based building block array).
    p_array=[]
    for i in range (0, len(g_nmr), 1):
        p_array.append(1)
        
    for i in range (0, len(g_nmr), 1):      
        if g_nmr[i][0][0]=='MET':
            p_array[i]=p_array[i]*(p_aa_met**4)
        if g_nmr[i][0][0]=='ILE':
            p_array[i]=p_array[i]*(p_aa_ile**4)
        if g_nmr[i][0][0]=='LEU':
            p_array[i]=p_array[i]*(p_aa_leu**4)
        if g_nmr[i][0][0]=='VAL':
            p_array[i]=p_array[i]*(p_aa_val**4)
        if g_nmr[i][0][0]=='ALA':
            p_array[i]=p_array[i]*(p_aa_ala**4)
        if g_nmr[i][0][0]=='THR':
            p_array[i]=p_array[i]*(p_aa_thr**4)
               
        for j in range (0, len(g_nmr[i]), 1):
            if g_nmr[i][j][2]=='MET':
                p_array[i]=p_array[i]*p_met
            if g_nmr[i][j][2]=='ILE':
                p_array[i]=p_array[i]*p_ile
            if g_nmr[i][j][2]=='LEU':
                p_array[i]=p_array[i]*p_leu
            if g_nmr[i][j][2]=='VAL':
                p_array[i]=p_array[i]*p_val
            if g_nmr[i][j][2]=='ALA':
                p_array[i]=p_array[i]*p_ala
            if g_nmr[i][j][2]=='THR':
                p_array[i]=p_array[i]*p_thr
    
    # interactive determination of which percentage of building block should be considered as potential starting points sorted by their uniqueness score. Result is an array with most unique NMR based building blocks.   
    percentage_of_considered_methyl_groups_as_start_points=float(percentage_starting_points)#input("Which percentage of total Methylgroups do you want to use as potential cluster starting points? The more you use, the more assignments will be possible but the longer it will take to execute the skript. Type value between 0 (0%) and 1 (100%):   ")
    #print('Please wait. Depending on your computional power and data, running the skript can take a few minutes to 24h \n ...\n ...')
    amount_of_unique_as=int(len(g_nmr)*float(percentage_of_considered_methyl_groups_as_start_points))
    most_unique_aa_array=[]
    for i in range (0, amount_of_unique_as, 1):
        p_most_unique=2
        most_unique_aa_pos='no one found'
        for j in range (0, len (p_array), 1):
             if 0<p_array[j]<p_most_unique:
                 p_most_unique=p_array[j]
                 most_unique_aa_pos=j
        most_unique_aa_array.append([g_nmr[most_unique_aa_pos], p_most_unique])
        p_array[most_unique_aa_pos]=2
    
    if most_unique_aa_pos=='no one found':
        return []  
          
    return most_unique_aa_array

"""
Method to create NMR based building blocks. Input is NOE list. Output is array with NMR based building blocks.
"""             
def create_G_nmr (a): 
    result_array=[]
    i=0
    while i<len(a):
        b=i
        loc=[]
        while b<len(a) and a[b].split()[1]==a[i].split()[1]:#es darf kein \n am ende des textdokumentes sein, sonst gibt es  "list index out of range" fehler hier
                b=b+1
        for f in range ( i, b, 1):
           loc.append([a[f].split()[0], a[f].split()[1], a[f].split()[2], a[f].split()[3]])
        result_array.append(loc)
        i=b-i+i
    return result_array

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
Method that creates array with additional predicted parameters. Input is array with additional parameters provided by the user. Output is formated array with additional parameters.
"""
def create_chemShift_PRE_PCS_RDC_xtal (chemShift_PRE_PCS_RDC_xtal):
	chemShift_PRE_PCS_RDC_array_xtal=[]
	for i in range (0, len(chemShift_PRE_PCS_RDC_xtal), 1):
		chemShift_PRE_PCS_RDC_array_xtal.append([chemShift_PRE_PCS_RDC_xtal[i].split()[0], chemShift_PRE_PCS_RDC_xtal[i].split()[1], float(chemShift_PRE_PCS_RDC_xtal[i].split()[2]), float(chemShift_PRE_PCS_RDC_xtal[i].split()[3]), float(chemShift_PRE_PCS_RDC_xtal[i].split()[4]), float(chemShift_PRE_PCS_RDC_xtal[i].split()[5])]) #as nr, as type, chemshift, pre,pcs,rdc
	return chemShift_PRE_PCS_RDC_array_xtal

"""
Method that creates array with additional measured parameters. Input is array with additional parameters provided by the user. Output is formated array with additional parameters.
"""	
def create_chemShift_PRE_PCS_RDC_measure (chemShift_PRE_PCS_RDC_measure):
	chemShift_PRE_PCS_RDC_array_measure=[]
	for i in range (0, len(chemShift_PRE_PCS_RDC_measure), 1):
		chemShift_PRE_PCS_RDC_array_measure.append([chemShift_PRE_PCS_RDC_measure[i].split()[0], chemShift_PRE_PCS_RDC_measure[i].split()[1], float(chemShift_PRE_PCS_RDC_measure[i].split()[2]), float(chemShift_PRE_PCS_RDC_measure[i].split()[3]), float(chemShift_PRE_PCS_RDC_measure[i].split()[4]), float(chemShift_PRE_PCS_RDC_measure[i].split()[5])]) #as nr, as type, chemshift, pre,pcs,rdc
	return chemShift_PRE_PCS_RDC_array_measure

"""
Method that identifies pdb methyl groups with same amnico acid type as a active node of a certain NMR building block. It excludes already assigned pdb methyl groups.
Input is an NMR building block (most_unique_aa), an array of pdb methyl groups and their position (pdb_array) and the partially reconstructed pdb and NMR based graphs.
Output is an array of pdb methyl groups and their position. Already assigned pdb methyl groups are excluded.
"""
def find_aa_with_aa_type_of_most_unique_aa_in_pdb (most_unique_aa, pdb_array, assignment_pair):
    if most_unique_aa==[] or pdb_array==[]:
        return []
    
    #looks for pdb methyls with same amino acid type as active node of NMR building block.
    most_unique_aa_type=most_unique_aa[0][0][0]
    aa_with_aa_type_of_most_unique_aa_array=[]
    for i in range(0, len(pdb_array), 1):
        if most_unique_aa_type==pdb_array[i][0]:
           aa_with_aa_type_of_most_unique_aa_array.append(pdb_array[i])
    
    #sorts out already assigned pdb methyls
    for i in range (0, len(assignment_pair), 1): #entfernen der interessanten as die schon zug sind
        for j in range (0, len( aa_with_aa_type_of_most_unique_aa_array), 1):
            if aa_with_aa_type_of_most_unique_aa_array[j]!=['']:
                if assignment_pair[i][1][0][1]== aa_with_aa_type_of_most_unique_aa_array[j][1]:
                    aa_with_aa_type_of_most_unique_aa_array[j]=['']
    aa_with_aa_type_of_most_unique_aa_array_neu = [feld for feld in  aa_with_aa_type_of_most_unique_aa_array if feld != ['']] 
    
    return aa_with_aa_type_of_most_unique_aa_array_neu

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
                
                    
  
"""
Determines assignment pair by comparison of NMR based building block to elements of array of pdb based building blocks. Calculates comparison score of assignment pair that was identified.
Input is NMR based building block (most_unique_aa), array of pdb based building blocks (int_aa), arrays with predicted and measured additional parameters (chemShift_PRE_PCS_RDC_array_measure,chemShift_PRE_PCS_RDC_array_xtal)
and an array with weight factors interactively provided by the user (weightArray). The Output is either an array with the NMR based building block and the best matching pdb based building block and their comparison score or just 'stop' as a sign that no assignment pair was identified.
"""
def best_fit (int_aa, most_unique_aa, chemShift_PRE_PCS_RDC_array_measure,chemShift_PRE_PCS_RDC_array_xtal , weightArray, penalty_array):
    #read weight factors
    weight_noe=weightArray[0]
    weight_chemShift=weightArray[1]
    weight_pre=weightArray[2]
    weight_pcs=weightArray[3]
    weight_rdc=weightArray[4]
    
    
    if most_unique_aa == [] or int_aa==[]:
        return 'stop'
    #compares all pdb based building blocks against the NMR based building block and looks for the best matching pair and the second best matching pair.
    max_score=-100000
    pos_best_fit='no one found'
    best_score=-100000   
    max_score_zwei=-100000
    pos_best_fit_zwei='no one found'
    for i in range (0, len(int_aa), 1):
        
        chemShift_mostUnique=0
        chemShift_intAA=0
        pre_mostUnique=0
        pre_intAA=0
        pcs_mostUnique=0
        pcs_intAA=0
        rdc_mostUnique=0
        rdc_intAA=0
        
        score=0
        score_tot=-100000000
        int_aa_loc=copy.deepcopy(int_aa)
        most_unique_aa_loc=copy.deepcopy(most_unique_aa[0])
        
        # compares amino acid types of neighboring amino acids of active nodes in the corresponding pdb and NMR building blocks.
        if len(int_aa[i])>1:
            for j in range (0, len(most_unique_aa[0]), 1):
                for k in range (0, len(int_aa[i])-1, 1): #-1 weil ja die letzte Zahl die dist ist
                    if int_aa_loc[i][k][2]== most_unique_aa_loc[j][2]:
                        score=score+1
                        most_unique_aa_loc[j][2]='---'
                        int_aa_loc[i][k][2]='+++'
        
        #compares additional parameters.
        intAAfoundInList=False
        mostUniqueAAfoundInList=False
        chemShift_intAA=False
        pre_intAA=False
        pcs_intAA=False
        rdc_intAA=False
        chemShift_mostUnique=False
        pre_mostUnique=False
        pcs_mostUnique=False
        rdc_mostUnique=False
  
        if len(int_aa[i])>1:         
            for j in range (0, len(chemShift_PRE_PCS_RDC_array_xtal),1):
                if chemShift_PRE_PCS_RDC_array_xtal[j][1]==int_aa[i][0][1]: #chemShift_PRE_PCS_RDC_array_xtal[j][0]==int_aa[i][0][0] and chemShift_PRE_PCS_RDC_array_xtal[j][1]==int_aa[i][0][1]:
                    if chemShift_PRE_PCS_RDC_array_xtal[j][2]!=999:
                        chemShift_intAA=chemShift_PRE_PCS_RDC_array_xtal[j][2]
                    if chemShift_PRE_PCS_RDC_array_xtal[j][3]!=999:
                        pre_intAA=chemShift_PRE_PCS_RDC_array_xtal[j][3]
                    if chemShift_PRE_PCS_RDC_array_xtal[j][4]!=999:
                        pcs_intAA=chemShift_PRE_PCS_RDC_array_xtal[j][4]
                    if chemShift_PRE_PCS_RDC_array_xtal[j][5]!=999:
                        rdc_intAA=chemShift_PRE_PCS_RDC_array_xtal[j][5]
                    intAAfoundInList=True
            for j in range (0, len(chemShift_PRE_PCS_RDC_array_measure),1):	
                if chemShift_PRE_PCS_RDC_array_measure[j][1]==most_unique_aa[0][0][1]: # chemShift_PRE_PCS_RDC_array_measure[j][0]==most_unique_aa[0][0][0] and chemShift_PRE_PCS_RDC_array_measure[j][1]==most_unique_aa[0][0][1]:
                    if chemShift_PRE_PCS_RDC_array_measure[j][2]!=999:
                        chemShift_mostUnique=chemShift_PRE_PCS_RDC_array_measure[j][2]
                    if chemShift_PRE_PCS_RDC_array_measure[j][3]!=999:
                        pre_mostUnique=chemShift_PRE_PCS_RDC_array_measure[j][3]
                    if chemShift_PRE_PCS_RDC_array_measure[j][4]!=999:
                        pcs_mostUnique=chemShift_PRE_PCS_RDC_array_measure[j][4]
                    if chemShift_PRE_PCS_RDC_array_measure[j][5]!=999:
                        rdc_mostUnique=chemShift_PRE_PCS_RDC_array_measure[j][5]
                    mostUniqueAAfoundInList=True
              
        if chemShift_intAA!=False and chemShift_mostUnique!=False:
            chemShift=abs(chemShift_intAA-chemShift_mostUnique)
        else:
            chemShift=penalty_array[0]
        if pre_intAA!=False and pre_mostUnique!=False:
            pre=abs(pre_intAA-pre_mostUnique)
        else:
            pre=penalty_array[1]
        if pcs_intAA!=False and pcs_mostUnique!=False:
            pcs=abs(pcs_intAA-pcs_mostUnique)
        else:
            pcs=penalty_array[2]
        if rdc_intAA!=False and rdc_mostUnique!=False:
            rdc=abs(rdc_intAA-rdc_mostUnique)
        else:
            rdc=penalty_array[3]
        
        #calculates comparison score based on neighboring amino acid types and additional parameters provided by the user.
        if mostUniqueAAfoundInList==True and intAAfoundInList==True:
            score_tot=weight_noe*(2*score-abs(score-(len(int_aa[i])-1))-abs(score-len(most_unique_aa[0]))) - weight_chemShift*chemShift- weight_pre*pre-weight_pcs*pcs-weight_rdc*rdc
        else:
            score_tot=weight_noe*(2*score-abs(score-(len(int_aa[i])-1))-abs(score-len(most_unique_aa[0])))-(weight_chemShift*penalty_array[0]+weight_pre*penalty_array[1]+weight_pcs*penalty_array[2]+weight_rdc*penalty_array[3])
            
        #memorizes the best fitting pdb based building block and second best fitting building block
            
        if score_tot>max_score and len(int_aa[i])>1:
            max_score_zwei=max_score
            max_score=score_tot
            pos_best_fit_zwei=pos_best_fit
            pos_best_fit=i
            best_score=score
        elif score_tot>max_score_zwei and len(int_aa[i])>1:
            max_score_zwei=score_tot
            pos_best_fit_zwei=i
    
    #If either the NMR based bulding block or the pdb based building block have more than 5 neighbors more as both of them have in common, the assignment will be rejected. If both together have 5 neighbors more than they have in common, the assignment will be rejected as well. If there is a second building block having another active node, the score has to be significantly bigger for the best matching one
    if pos_best_fit!='no one found':
        return [most_unique_aa, int_aa[pos_best_fit], max_score] 
    else:
        return 'stop'

"""
Method that creates array with pdb methyl groups and their positions corresponding to nodes that are neighbors of a certain already assigned node in the partially reconstructed pdb based graph.
Inputs are partially reconstructed pdb based and NMR based graphs (assignment_pair), array with all methyl groups and their positions (pdb_array), a specific building block of the partially reconstructed pdb based graph and the corresponding building block in the partially reconstructed NMR based graph (assignment_pair_spez),
and the most unique building block corresponding to a node of the partially reconstructed NMR based graph. 
Output is an array with pdb methyl groups and their positions that ar neighbored to the specific already assigned node of the partially reconstructed pdb based graph and are not already assigned.
This backtracking is imprortant for the methyl walk.
"""
def find_neighbors_of_assigned_aa_crystal (assignment_pair, pdb_array, assignment_pair_spez, most_unique_aa_umgebung_nmr):    
    if most_unique_aa_umgebung_nmr==[]:
        return []
    
    #identify amino acid type of active node (this one is the already assigned one) in NMR building block.
    as_type_most_unique_nmr=most_unique_aa_umgebung_nmr[0][0][0]
    
    #looks for neighbors of the corresponding specific already assigned node of the partially reconstructed pdb based graph.
    int_aa_pdb_umgebung_crystal=[]
    assigned_crystal=assignment_pair_spez[1][0:-1]
    for i in range (0, len(assigned_crystal), 1):
        for j in range (0, len(pdb_array), 1):
            if pdb_array[j][1]==assigned_crystal[i][3] and as_type_most_unique_nmr==pdb_array[j][0]:
                int_aa_pdb_umgebung_crystal.append(pdb_array[j])
    
    #sorts out already assignedpdb methyl groups.
    for i in range (0, len(assignment_pair), 1):
        for j in range (0, len(int_aa_pdb_umgebung_crystal), 1):
            if int_aa_pdb_umgebung_crystal[j]!=['']:
                if assignment_pair[i][1][0][1]==int_aa_pdb_umgebung_crystal[j][1]:
                    int_aa_pdb_umgebung_crystal[j]=['']
    int_aa_pdb_umgebung_crystal_neu = [feld for feld in int_aa_pdb_umgebung_crystal if feld != ['']] 
    
    return int_aa_pdb_umgebung_crystal_neu

"""
Method that creates array with NMR based building blocks corresponding to nodes beeing neighbors of an active node in a specific NMR based bulding block. Input is array with NMR based building blocks (g_nmr), partially reconstructed NMR based graph and pdb based graph (assignment_pair)
and a specific NMR based building block (assignment_pair_spez). Output is array with bulding block corresponding to nodes beeing part of the NMR based building block (except for the specific active node).
"""
def find_neighbors_of_assigned_nmr(assignment_pair, g_nmr, assignment_pair_spez):
    #all building blocks corresponding to nodes beeing neighbors of the active node in the specific NMR based buidling block are identified in the array of all NMR based building blocks and memorized.
    int_aa_umgebung_nmr=[]
    assigned_nmr=assignment_pair_spez[0][0]
    for i in range (0, len(assigned_nmr), 1):
        for j in range (0, len(g_nmr), 1):
            if g_nmr[j][0][1]==assigned_nmr[i][3]:
                int_aa_umgebung_nmr.append(g_nmr[j])
    #all building blocks identified in the previous step that correspond to alraedy assigned nodes are sorted out
    for i in range (0, len(assignment_pair), 1):
        for g in range (0, len(int_aa_umgebung_nmr), 1):
            if int_aa_umgebung_nmr[g]!=['']:
                if assignment_pair[i][0][0][0][1]==int_aa_umgebung_nmr[g][0][1]:
                    int_aa_umgebung_nmr[g]=['']
    int_aa_umgebung_nmr_neu = [feld for feld in int_aa_umgebung_nmr if feld != ['']]  
    return int_aa_umgebung_nmr_neu

"""
Method that sorts array of NMR building blocks by their uniqueness score. Input is an array of NMR building blocks. Output is array of NMR building blocks sorted by their uniqueness score.
"""
def find_most_unique_aa_in_neighbors (g_nmr):
    if g_nmr==[]:
        return []
    
    #amount of active nodes having a specific amino acid type as a neighbor
    amount_met=0
    amount_ile=0
    amount_leu=0
    amount_val=0
    amount_ala=0
    amount_thr=0
    amount_total=0
    
    #amount of active nodes having a specific amino acid type
    aa_met=0
    aa_ile=0
    aa_leu=0
    aa_val=0
    aa_ala=0
    aa_thr=0
    amount_aa_total=0
    
    for i in range (0, len(g_nmr), 1):
        
        #counts event of active nodes having a certain amico acid type as a neighbor.
        for j in range (0, len(g_nmr[i]), 1):
            if g_nmr[i][j][2]=='MET':
                amount_met=amount_met+1
                amount_total=amount_total+1
            if g_nmr[i][j][2]=='ILE':
                amount_ile=amount_ile+1
                amount_total=amount_total+1
            if g_nmr[i][j][2]=='LEU':
                amount_leu=amount_leu+1
                amount_total=amount_total+1
            if g_nmr[i][j][2]=='VAL':
                amount_val=amount_val+1
                amount_total=amount_total+1
            if g_nmr[i][j][2]=='ALA':
                amount_ala=amount_ala+1
                amount_total=amount_total+1
            if g_nmr[i][j][2]=='THR':
                amount_thr=amount_thr+1
                amount_total=amount_total+1

        #counts active nodes having a certain amino acid type.
        if g_nmr[i][0][0]=='MET':
            aa_met=aa_met+1
            amount_aa_total=amount_aa_total+1
        if g_nmr[i][0][0]=='ILE':
            aa_ile=aa_ile+1
            amount_aa_total=amount_aa_total+1
        if g_nmr[i][0][0]=='LEU':
            aa_leu=aa_leu+1
            amount_aa_total=amount_aa_total+1
        if g_nmr[i][0][0]=='VAL':
            aa_val=aa_val+1
            amount_aa_total=amount_aa_total+1
        if g_nmr[i][0][0]=='ALA':
            aa_ala=aa_ala+1
            amount_aa_total=amount_aa_total+1
        if g_nmr[i][0][0]=='THR':
            aa_thr=aa_thr+1
            amount_aa_total=amount_aa_total+1
    
    #probability of an active node having a certain amino acid type.
    p_aa_met=aa_met/amount_aa_total
    p_aa_ile=aa_ile/amount_aa_total
    p_aa_leu=aa_leu/amount_aa_total
    p_aa_val=aa_val/amount_aa_total
    p_aa_ala=aa_ala/amount_aa_total
    p_aa_thr=aa_thr/amount_aa_total
    
    #probability of an active node having a certain amino acid type as a neighbor
    p_met=amount_met/amount_total
    p_ile=amount_ile/amount_total
    p_leu=amount_leu/amount_total
    p_val=amount_val/amount_total
    p_ala=amount_ala/amount_total
    p_thr=amount_thr/amount_total
    
    #calculation of uniqueness score for active nodes. Result is p_array. The i-th element in p_array corresponds to the i-th building block in g_nmr (NMR based building block array).
    p_array=[]
    for i in range (0, len(g_nmr), 1):
        p_array.append(1)

    for i in range (0, len(g_nmr), 1):   
        if g_nmr[i][0][0]=='MET':
            p_array[i]=p_array[i]*(p_aa_met**4)
        if g_nmr[i][0][0]=='ILE':
            p_array[i]=p_array[i]*(p_aa_ile**4)
        if g_nmr[i][0][0]=='LEU':
            p_array[i]=p_array[i]*(p_aa_leu**4)
        if g_nmr[i][0][0]=='VAL':
            p_array[i]=p_array[i]*(p_aa_val**4)
        if g_nmr[i][0][0]=='ALA':
            p_array[i]=p_array[i]*(p_aa_ala**4)
        if g_nmr[i][0][0]=='THR':
            p_array[i]=p_array[i]*(p_aa_thr**4)
       
        for j in range (0, len(g_nmr[i]), 1):
            if g_nmr[i][j][2]=='MET':
                p_array[i]=p_array[i]*p_met
            if g_nmr[i][j][2]=='ILE':
                p_array[i]=p_array[i]*p_ile
            if g_nmr[i][j][2]=='LEU':
                p_array[i]=p_array[i]*p_leu
            if g_nmr[i][j][2]=='VAL':
                p_array[i]=p_array[i]*p_val
            if g_nmr[i][j][2]=='ALA':
                p_array[i]=p_array[i]*p_ala
            if g_nmr[i][j][2]=='THR':
                p_array[i]=p_array[i]*p_thr
    
    #sorting NMR building blocks for their uniquenes scores
    amount_of_unique_as=int(len(g_nmr))
    most_unique_aa_array=[]
    for i in range (0, amount_of_unique_as, 1):
        p_most_unique=2
        most_unique_aa_pos='no one found'
        for j in range (0, len (p_array), 1):
             if 0<p_array[j]<p_most_unique:
                 p_most_unique=p_array[j]
                 most_unique_aa_pos=j
        most_unique_aa_array.append([g_nmr[most_unique_aa_pos], p_most_unique])
        p_array[most_unique_aa_pos]=2
    
    if most_unique_aa_pos=='no one found':
        return []  
          
    return most_unique_aa_array

"""
Method to evaluate assignments. This is important for parameter optimization. Inputs are the final pdb and NMR based graphs with the assignments (best_assignment_pair).
Outputs are the same as the inputs but with additional notes (+,-,o) that will appear in the output file (result.txt). "+" does mean the edge pattern of two certain building blocks match very well.
'o' means that the edge pattern doesn't match 100 percent but still not too bad. "-" means that the edge pattern doesn't match very well.
"""
def evaluate_assignment_pair(best_assignment_pair,start):
    
    #compare edge patterns of pdb and NMR building blocks that are assigned to each other and classify the matching.     
    for i in range (0, len (best_assignment_pair),1 ):
        score=0
        assignment_pair_loc=copy.deepcopy(best_assignment_pair[i])     
        for j in range (0, len (assignment_pair_loc[0][0]),1):
            for k in range (0, len(assignment_pair_loc[1])-1, 1): 
                if assignment_pair_loc[0][0][j][2]== assignment_pair_loc[1][k][2]:
                    score=score+1
                    assignment_pair_loc[0][0][j][2]='---'
                    assignment_pair_loc[1][k][2]='+++'
        if abs(score-(len(assignment_pair_loc[1])-1))<1 and abs(score-len (assignment_pair_loc[0][0]))<1 and abs(2*score-((len(assignment_pair_loc[1])-1)+len (assignment_pair_loc[0][0])))<2:
            best_assignment_pair[i][0].append('+')
        elif abs(score-(len(assignment_pair_loc[1])-1))<2 and abs(score-len (assignment_pair_loc[0][0]))<2 and abs(2*score-((len(assignment_pair_loc[1])-1)+len (assignment_pair_loc[0][0])))<2:
            best_assignment_pair[i][0].append('o')
        else:
            best_assignment_pair[i][0].append('-')
            

if __name__ == "__main__":
    cli()

