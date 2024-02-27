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
def cli (results_files, outfile_name):
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
                 and int(result_file_lines[j].split(":")[1]) <= 5):
                cluster_size = int(result_file_lines[j].split(":")[1])
                del assignment_list[-cluster_size:]

        
        
        columns=["nmr_aa_type", "nmr_id", "pdb_aa_type", "pdb_id"]
        df_temp=pd.DataFrame(assignment_list, columns=columns)
        all_assignments.append(df_temp)
        #print (df_temp)

        #print(assignment_list)
    print (all_assignments)
    df_final_assignment=all_assignments[0]
    print ("assignments using 1 file")
    print (df_final_assignment.shape[0])
    for i in range (1, len(all_assignments), 1):
        df_final_assignment=pd.merge(df_final_assignment, all_assignments[i])
        print ("assignments using "+ str(i+1) + " file")
        print (df_final_assignment.shape[0])
    
    # safe file
    df_final_assignment.to_csv(outfile_name, sep="\t") 

        


    #for 







if __name__ == "__main__":
    cli()