import click
import pandas as pd
"""
define input parameters
"""
@click.command()
@click.option(
    "--outfile_name",
    "outfile_name",
    required=True,
    multiple=False,
    type=str,
    default="final_assignment.tsv",
    help="name of outfile",
)
@click.option(
    "--results_files",
    "results_files",
    required=True,
    multiple=True,
    type=str,
    default=["./10/result.txt", "./11/result.txt"],
    help="path to file",
)
@click.option(
    "--noe_file",
    "noe_file",
    required=True,
    multiple=False,
    type=str,
    default="./7/noe.txt",
    help="path to noe file to get number of assignable resonances",
)
@click.option(
    "--min_cluster_size",
    "min_cluster_size",
    required=True,
    multiple=False,
    type=int,
    default=5,
    help="minimal cluster size to consider a particular cluster for assigning methyl groups",
)
def cli (results_files, outfile_name, min_cluster_size, noe_file):
    if len(results_files)!=3:
        raise ValueError('please provide exactly three assignment files')
    
    all_assignments=[]
    # iterate over result files
    for i in range (0, len(results_files), 1):
        # assignment list
        assignment_list=[]

        result_file=open(results_files[i], "r")
        result_file_lines=result_file.readlines()
        result_file.close()
        # iterate over lines
        for j in range (0, len(result_file_lines), 1):
            if result_file_lines[j].split(":")[0]=="was extened with methyl group resonance  ":
                assignment_list.append([result_file_lines[j].split(":")[1].split()[0], result_file_lines[j].split(":")[1].split()[1],
                                       result_file_lines[j+1].split(":")[1].split()[0], result_file_lines[j+1].split(":")[1].split()[1]])
            
            # remove cluster assignments if cluster is smaller than 5
            if (result_file_lines[j].split(":")[0]=="amount of assignments in cluster                  "
                 and int(result_file_lines[j].split(":")[1]) <= min_cluster_size):
                cluster_size = int(result_file_lines[j].split(":")[1])
                del assignment_list[-cluster_size:]

        
        
        columns=["nmr_aa_type", "nmr_id", "pdb_aa_type", "pdb_id"]
        df_temp=pd.DataFrame(assignment_list, columns=columns)
        all_assignments.append(df_temp)
        #print (df_temp)

        #print(assignment_list)
    # generate df with unambiguous assignments
    df_final_assignment=all_assignments[0]
    print ("unambiguous assignments using 1 file")
    print (df_final_assignment.shape[0])
    for i in range (1, len(all_assignments), 1):
        df_final_assignment=pd.merge(df_final_assignment, all_assignments[i])
        print ("unambiguous assignments using "+ str(i+1) + " file")
        print (df_final_assignment.shape[0])

    # generate df with all assignments
    df_ambiguous_assignments=all_assignments[0]
    df_ambiguous_assignments = df_ambiguous_assignments.rename(columns={'pdb_aa_type': 'pdb_aa_type_'+str(results_files[0]), 'pdb_id': 'pdb_id_'+str(results_files[0])})
    for i in range (1, len(all_assignments), 1):
        all_assignments[i] = all_assignments[i].rename(columns={'pdb_aa_type': 'pdb_aa_type_'+str(results_files[i]), 'pdb_id': 'pdb_id_'+str(results_files[i])})
        df_ambiguous_assignments=pd.merge(df_ambiguous_assignments, all_assignments[i], on=["nmr_id", "nmr_aa_type"], how="outer")
        #print ("assignments using "+ str(i+1) + " file")
        #print (df_ambiguous_assignments.shape[0])
    
    
    # process
    df_ambiguous_assignments=df_ambiguous_assignments[(df_ambiguous_assignments["pdb_id_"+str(results_files[0])] != df_ambiguous_assignments["pdb_id_"+str(results_files[1])]) 
                                                      | (df_ambiguous_assignments["pdb_id_"+str(results_files[0])] != df_ambiguous_assignments["pdb_id_"+str(results_files[2])])
                                                      | (df_ambiguous_assignments["pdb_id_"+str(results_files[1])] != df_ambiguous_assignments["pdb_id_"+str(results_files[2])])]
    print ("number of ambiguous assignments :")
    print (df_ambiguous_assignments.shape[0])
    

    # generate list of unassigned methyls
    noe_list_file=open(noe_file, "r")
    noe_list_lines=noe_list_file.readlines()
    noe_list_file.close()
    noe_list=[]

    # generate df with all methyls
    ## iterate over lines
    for i in range (0, len(noe_list_lines), 1):
        noe_list.append(noe_list_lines[i].split())
    
    
    ## generate df
    columns=["nmr_aa_type", "nmr_id", "nmr_aa_type_acceptor", "nmr_id_acceptor"]
    df_noes=pd.DataFrame(noe_list, columns=columns)
    df_noes=df_noes.drop_duplicates(subset=["nmr_aa_type", "nmr_id"]) #noe donors
    df_noes = df_noes.drop(columns=["nmr_aa_type_acceptor", "nmr_id_acceptor"], axis=1)
    print("assignabel methyls :")
    print ( df_noes.shape[0])
    

    ## only keep unassigned methyls
    df_all_assignments=pd.concat(all_assignments, axis=0)
    df_all_assignments=df_all_assignments.drop_duplicates(subset=["nmr_aa_type", "nmr_id"])
    df_all_assignments = df_all_assignments.drop(columns=["pdb_aa_type", "pdb_id"], axis=1)
    #print(df_all_assignments)

    common=df_noes.merge(df_all_assignments, on=["nmr_aa_type", "nmr_id"])
    df_unanssigned=df_noes[~df_noes.nmr_id.isin(common.nmr_id)]
    #print (common)

    print ("number of unassigned methyls: ")
    print (df_unanssigned.shape[0])


    # safe files
    df_final_assignment.to_csv(outfile_name, sep="\t") 
    df_ambiguous_assignments.to_csv("ambiguous_"+outfile_name, sep="\t") 

    df_unanssigned.to_csv("unassigned_"+outfile_name, sep="\t") 



        


    #for 







if __name__ == "__main__":
    cli()