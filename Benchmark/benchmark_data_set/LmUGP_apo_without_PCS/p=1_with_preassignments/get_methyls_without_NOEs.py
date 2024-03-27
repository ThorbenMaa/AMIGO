import click
import pandas as pd
"""
define input parameters
"""
@click.command()
@click.option(
    "--shared_ass_file",
    "shared_ass_file",
    required=True,
    multiple=False,
    type=str,
    default="shared_min_cluster_1.tsv",
    help="name of outfile",
)
@click.option(
    "--amb_ass_file",
    "amb_ass_file",
    required=True,
    multiple=False,
    type=str,
    default="ambiguous_shared_min_cluster_1.tsv",
    help="name of outfile",
)
@click.option(
    "--not_ass_file",
    "not_ass_file",
    required=True,
    multiple=False,
    type=str,
    default="unassigned_shared_min_cluster_1.tsv",
    help="name of outfile",
)
@click.option(
    "--pdb_file",
    "pdb_file",
    required=True,
    multiple=False,
    type=str,
    default="shared_min_cluster_1.tsv",
    help="name of outfile",
)
@click.option(
    "--labeling_scheme",
    "labeling_scheme",
    required=True,
    multiple=True,
    type=str,
    default=["M", "I", "L", "V", "A"],#, "T"],
    help="name of outfile",
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

def cli (shared_ass_file, amb_ass_file, not_ass_file, pdb_file, labeling_scheme, val_scheme, leu_scheme):
    
    # load AMIGO results
    shared_ass_df=pd.read_csv(shared_ass_file, sep="\t")
    amb_ass_df=pd.read_csv(amb_ass_file, sep="\t")
    not_ass_df=pd.read_csv(not_ass_file, sep="\t")

    number_of_aa_with_noes = shared_ass_df.shape[0] + amb_ass_df.shape[0] + not_ass_df.shape[0]
    print("number of methyls with NOEs: ", number_of_aa_with_noes)


    # create pdb list
    met=False
    ile=False
    leu=False
    val=False
    ala=False
    thr=False
    for labeled in labeling_scheme:
        if labeled=="M":
            met=True
        if labeled=="I":
            ile=True
        if labeled=="L":
            leu=True
        if labeled=="V":
            val=True
        if labeled=="A":
            ala=True
        if labeled=="T":
            thr=True
    #print(met, ile, leu, val, ala, thr, pdb_file, val_scheme, leu_scheme)
    pdb_list = create_pdb_array( met, ile, leu, val, ala, thr, pdb_file, val_scheme, leu_scheme )
    print ("number of labelled methyls in pdb", len(pdb_list))

    #print (shared_ass_df["pdb_id"].values)
    #print(amb_ass_df.iloc[:,4].values)
    #print(not_ass_df["nmr_id"].values)

    # identify methyl groups without NOEs; i.e., unconsidered by AMIGO
    not_found=[]
    for i in range (0, len(pdb_list), 1):    
        pdb_id = pdb_list[i][1]
        aa_type = pdb_list[i][0]
        #print(pdb_list[i])
        #print(pdb_id)
        found = False

        if int(pdb_id) in shared_ass_df["nmr_id"].values :
            found = True
        if int(pdb_id) in amb_ass_df["nmr_id"].values :
            found = True
        if int(pdb_id) in not_ass_df["nmr_id"].values :
            found = True
        
        # save pdb_id if not found
        if found == False:
            not_found.append(str(pdb_id)+" "+str(aa_type))

    not_found = list(set(not_found))
    print("number of methyls without NOE", len(not_found))

    outfile = open("no_NOEs_or_no_reference_assignment.txt", "w")
    for i in range (0, len(not_found), 1):
        outfile.write(not_found[i].replace("[", "").replace("]", "").replace(",", "")+"\n")

        
        
        

    
    



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






if __name__ == "__main__":
    cli()