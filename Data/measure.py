import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import re
from collections import deque
import sys


def comp_within(mutant_df, control_df, dist:int=10) -> pd.DataFrame:
    """
    Compare positions between mutant and control dataframes and remove matched rows from mutant dataframe.
    This function iterates through the unique chromosomes in the mutant dataframe and compares the positions
    of the rows in the mutant dataframe with the positions in the control dataframe within a specified distance.
    If a match is found within the distance, the row is removed from the mutant dataframe.
    Parameters:
        mutant_df (pd.DataFrame): DataFrame containing mutant data with at least 'chrom' and position columns.
        control_df (pd.DataFrame): DataFrame containing control data with at least 'chrom' and position columns.
        dist (int, optional): Distance within which to consider positions as matching. Default is 10.
    Returns:
    pd.DataFrame: The mutant dataframe with matched rows removed.
    """
    indices_to_drop = []
    
    # Cycle through each chromosome in the mutant dataframe
    for chrom in mutant_df['#CHROM'].unique():
        control_df_chrom = control_df[control_df['#CHROM'] == chrom]
        
        # Cycle through the rows in the mutant dataframe
        for i in range(0, mutant_df.shape[0], dist):
            pos = mutant_df.iloc[i, 1]
            
            matched = False
            j = 0
            # Cycle through the rows in the control dataframe
            while (j < control_df_chrom.shape[0]) & (matched == False):
                pos_2 = control_df_chrom.iloc[j, 1]
                # Check if the positions match within the specified distance, and if so, add the index to the list
                if pos - dist <= pos_2 <= pos + dist:
                    indices_to_drop.append(mutant_df.index[i])
                    matched = True
                    
                j += 1
                
    return mutant_df

def compare_between(all_dict, distance:int=10, agreement:int = 2) -> pd.DataFrame:
        """
        Compare positions between different methods and return a combined dataframe.
        This function iterates through the methods in the dictionary and compares the positions
        of the rows in the dataframes within a specified distance. If a match is found within the distance
        in at least the specified number of methods, the row is included in the combined dataframe.
        Parameters:
            all_dict (dict): Dictionary containing dataframes for different methods.
            distance (int, optional): Distance within which to consider positions as matching. Default is 10.
            agreement (int, optional): Minimum number of methods that should agree on a position. Default is 2.
        Returns:
        pd.DataFrame: The combined dataframe with matched rows.
        """
        consensus = []
       
        # making a dictionary to remember which column each pos is stored in.
        pos_col_dict = {}
        i = 1
        for key in all_dict.keys():
            pos_col_dict[key] = i
            i += 1
        
        # Cycle every method
        for curr_method in all_dict:
            df_keys = list(all_dict.keys())
            # df_keys.remove(curr_method)
            leader = all_dict[curr_method]
            
            # Cycle every row in the leader
            for i in range(0, leader.shape[0]):
                chrom = leader.iloc[i]['#CHROM']
                pos = leader.iloc[i]['POS']
                match_row = [-1 for x in range(len(all_dict.keys()) + 2)]
                count = 0
                average = 0
                
                # Cycle every method
                for alt_method in df_keys:
                    matches = all_dict[alt_method]
                    
                    # Filter for chromosome matches 
                    matches = matches[matches['#CHROM'] == chrom]
                    matches = matches[(matches['POS'] >= pos - distance) & (matches['POS'] <= pos + distance)]
                    matches['POS'] = matches['POS'].astype(int)
                
                    # If there are matches, store the position and update the running sums
                    if not matches.empty:
                        storage_pos = pos_col_dict[alt_method]
                        # store the position
                        match_row[storage_pos] = tuple(matches['POS'].tolist())
                        
                        # update the running sums
                        average += sum(matches['POS'].tolist())
                        count += 1
                
                # If the count is greater than the agreement, add the row to the consensus
                if count >= agreement:
                    # remember the chromosome
                    match_row[0] = chrom
                    consensus.append(match_row)
                    average = average / count
                    # remember the average
                    match_row[-1] = int(average)
                    
        
        columns = ['CHROM'] + list(pos_col_dict.keys()) + ['AVERAGE']
        
        consensus_df = pd.DataFrame(consensus, columns = columns)
        consensus_df = consensus_df.drop_duplicates()
        
        consensus_df = consensus_df.sort_values(by='CHROM', ascending=True)
        
        return consensus_df
                    

class file:
    """
    A class to represent a file with specific attributes and methods.

    Attributes:
    -----------
        path : str
            The file path.
        method : str
            The method associated with the file.
        condition : str
            The condition associated with the file.
        df : DataFrame
            The DataFrame associated with the file.

    Methods:
    --------
        __init__(self, path, method, condition, df)
            Initializes the file with the given path, method, condition, and DataFrame.
    """
    def __init__(self, path, method, condition, df):
        self.method = method
        self.condition = condition
        self.path = path
        self.df = df
        
        
def get_data(directory):
    """
    Reads CSV files from a specified directory, extracts method and condition from filenames,
    and returns a list of file objects containing the file path, method, condition, and DataFrame.
    Args:
        directory (str): The path to the directory containing the CSV files.
    Returns:
        list: A list of file objects, each containing:
            - file_path (str): The path to the CSV file.
            - method (str): The method extracted from the filename.
            - condition (str): The condition extracted from the filename.
            - df (pd.DataFrame): The DataFrame containing the data from the CSV file.
    """
    files = []
    
    # Iterate over the files in the directory
    for filename in os.listdir(directory):
        # Check if the file is a CSV file
        if filename.endswith('.csv'):
            # Extract the method and condition from the filename
            method = re.search(r'(?<=_).+(?=.csv)', filename).group(0)
            condition = re.search(r'.+(?=_)', filename).group(0)
            
            # print(f'File Name: {filename}, Method: {method}, Condition: {condition}')

            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path)
            
            files.append(file(file_path, method, condition, df))
            
    return files
            
def filter(files):
    """
    Filters a list of file objects based on specific criteria.
    This function iterates over a list of file objects and filters out those
    whose 'method' attribute is 'seekSV'. For the remaining files, it further
    filters the DataFrame (df attribute) to include only rows where the 'FILTER'
    column has the value 'PASS'. The filtered file objects are then returned.
    Args:
        files (list): A list of file objects. Each file object is expected to have
                      a 'method' attribute and a 'df' attribute, where 'df' is a
                      pandas DataFrame.
    Returns:
        list: A list of filtered file objects.
    """
    
    filtered_files = []
    
    for file in files:
        # For now, we are only interested in files that are not seekSV
        if file.method != 'seekSV':
            df = file.df
            # Filter out rows with quality flags
            df = df[(df['FILTER'] == 'PASS') | (df['FILTER'] == '.')]
            file.df = df
            filtered_files.append(file)
            if file.df.empty == False:
                printdf = file.df[['#CHROM', 'POS']]
                print(f'File: {file.path} ########## \n {printdf}')
    return filtered_files


if __name__ == "__main__":
    directory = './Data'
    if len(sys.argv) == 3:
        distance = int(sys.argv[1])
        agreement = int(sys.argv[2])
    else:
        distance = int(input("Enter the distance: "))
        agreement = int(input("Enter the agreement: "))
        
        if distance < 1 or agreement < 1:
            print("Distance and agreement must be positive integers.")
            sys.exit(1)
    
    
    # distance = 30
    # agreement = 3
    
    files = get_data(directory)
    files = filter(files)
    
    filtered_within ={}
    
    methods = set(file.method for file in files)
    # Iterate over the methods to grab the N2 and GE2722 dataframes
    for method in methods:
        n2_df, ge_df = None, None  
        for file in files:
            if (file.method == method) & (file.condition == 'N2'):
                N2_df = file.df
            elif (file.method == method) & (file.condition == 'GE2722'):
                GE_df = file.df
        if N2_df is None or GE_df is None:
            print(f'No matching condition found for{method}')
            continue
        # Compare GE with N2 and add filtered file to the dictionary
        filtered_within[method] = comp_within(GE_df, N2_df, distance)
        
    
    final_consensus = compare_between(filtered_within, distance, agreement)
    
    print(final_consensus)
    output_path = 'final_consensus.csv'
    final_consensus.to_csv(output_path, index=False)


    
        