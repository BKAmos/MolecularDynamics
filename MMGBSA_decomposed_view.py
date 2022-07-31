# -*- coding: utf-8 -*-

###############################################################################################
# Authors: B Kirtley Amos and Owen Yao
# Email ID: k.amos@atombioworks.com
# Usage: python script (The upper directory variable needs to be changed for each new run)
# Outputs: Produces heatmaps of models as well as the averaged heat map for models and outputs them in a HeatMap directory
# Output format: Heatmaps in .png for each model and an average heatmap of all models
# Caveats: In its current state this only works for averages when the time frames are specific (45, 95, 145, 195)
###############################################################################################

import matplotlib.pyplot as plt
import pandas as pd
pd.options.mode.chained_assignment = None
import seaborn as sns
import glob
import os


def heatmapper(data, directory='N/A', typ='N/A'):
    # Formats data
    data = pd.concat(data)
    data['Time'] = [int(x) for x in data['Time']]
    data['TotalE'] = [float(x) for x in data['TotalE']]
    top5 = data[data['Ranking'] <= 5]
    
    # Create heatmaps folder
    try:
        os.mkdir(directory+'/'+'heatmaps')
    except FileExistsError:
        pass
    
    # Formats heatmap
    pivotTabel = top5.pivot('Concat', 'Time', 'TotalE')
    pivotTabel.reindex(sorted(pivotTabel.columns), axis=1)
    sns.heatmap(pivotTabel, annot=True, linewidths=.5, cbar_kws={'label': 'kcal/mol'}) \
    .set(ylabel='Residue pairs', xlabel='Time (ns)')
    
    # Display/save heatmap
    # plt.savefig(f'{directory}/heatmaps/{directory.split("/")[-1]}_{typ}', bbox_inches='tight')
    plt.show()
    
    
def averaging(data):
    ulty = []
    
    # Put all trials together, calculate mean per concat
    for time, files in data.items():
        for x in range(len(files)):
            try:
                df_m = pd.merge(df_m, files[x].iloc[:, 2:4], on="Concat", how='outer').fillna(0)
            except UnboundLocalError:
                df_m = files[0].iloc[:, 2:4]

        df_m['Resnum']=df_m['Concat'].str.split(':').str[1]
        df_m.sort_values(by=['Resnum'], inplace=True, ascending=False)
        df_m.iloc[:,1:len(df_m.columns)-1] = df_m.iloc[:,1:len(df_m.columns)-1].astype(float)
        df_m['Mean']=df_m.iloc[:,1:len(df_m.columns)-1].mean(axis=1)
        df_m['Time']=time
        ulty.append(df_m)
        del df_m

    # Merge concat results accordingly
    for dfx in ulty:
        try:
            df_means = pd.merge(df_means, dfx.loc[:, ['Concat', 'Mean', 'Time']], 
                                on='Concat', how='outer').fillna(0)
        except UnboundLocalError:
            df_means = dfx.loc[:, ['Concat', 'Mean', 'Time']]

    df_means['Resnum']=df_means['Concat'].str.split(':').str[1]
    df_means.sort_values(by=['Resnum'], inplace=True, ascending=False)
    return df_means


# Retrieve model csv file
upperDirectory = '/home/abwer/owen/mmgbsa_pipeline_refinement/rbd_apt25_model1/rbd_apt25_model1'
text_files = glob.glob(upperDirectory + '/rbd*/*.dat', recursive = True)
timecheck = {}
dataframeList = []


# Loop through model csv file
for i in [x for x in text_files if 'top20' not in x]:
    # Open file
    try:
        df = pd.read_csv(i, header=None)
        basename = os.path.basename(i)

    except pd.errors.ParserError:
        with open(i) as file:
            lines = file.readlines()
            start = lines.index('Complex:\n')
            end = lines[start:].index('Sidechain Energy Decomposition:\n')
            fix = lambda x: [z.rstrip('\n') for z in x]
            lines = [fix(x.split(',')) for x in lines[start+4:][:end-5]]
            df = pd.DataFrame(lines)
            basename = i.split('/')[-2]

    # Retrieve first 20 sorted rows with r:a: in 0 and l:b: in 1
    samesRemoved = df[(df[0] != df[1])]
    droppingOpposites = samesRemoved.loc[samesRemoved[0].str.contains('R:A:') 
                                     == samesRemoved[1].str.contains('L:B:')]
    droppingLB = droppingOpposites[droppingOpposites[0].str.contains('L:B:') == False]
    droppingLB[list(range(2, 20))] = droppingLB[list(range(2, 20))].astype(float)
    sorted_df = droppingLB.sort_values(by = [17], ascending=True).head(20)

    # Get the model name
    stripped = basename.split('.', 1)[0]
    stripped2 = stripped.split('-', 1)[0]
    ranking = []

    # Rank the things
    for j in range(len(sorted_df)):
        ranking.append(j + 1)

    # Create new df with labels
    sorted_df['Ranking'] = ranking
    sorted_df['Model'] = stripped2    
    df2 = pd.DataFrame().assign(Res1=sorted_df[0], Res2=sorted_df[1], 
                                Concat=[x.replace('R:A:', '') + ":" + y.replace('L:B:', '')
                                        for x,y in zip(sorted_df[0], sorted_df[1])],
                                TotalE=sorted_df[17], Ranking=sorted_df['Ranking'], 
                                Model=sorted_df['Model'], 
                                Time=[x.split('_')[-1] for x in sorted_df['Model']])
    dataframeList.append(df2)

    # Gather data for average
    try:
        timecheck[stripped2.split('_')[-1]].append(df2)
    except KeyError:
        timecheck[stripped2.split('_')[-1]] = [df2]

# Generates heatmaps
first = []
second = []
third = []
fdf = pd.DataFrame(columns=['Concat', 'TotalE', 'Time'])
data = [first.append(x) if x.iloc[0]['Model'].split('_')[-2] == '1' 
        else second.append(x) if x.iloc[0]['Model'].split('_')[-2] == '2'
        else third.append(x) for x in dataframeList]
heatmapper(first, upperDirectory, 1)
heatmapper(second, upperDirectory, 2)
heatmapper(third, upperDirectory, 3)
averages = averaging(timecheck)

# Reorder df into rows with concat, mean, and time
for i in range(len(averages)):
    st = i
    i *= 4
    for x in [1, 3, 5, 7]:
        fdf.loc[i] = [averages.iloc[st, 0], averages.iloc[st, x], averages.iloc[st, x+1]]
        i += 1

# Redistribute df into time frame dfs
fdf = fdf[fdf["Time"] != 0]
t45, t95, t145, t190 = fdf[fdf["Time"] == '45'], fdf[fdf["Time"] == '95'], \
                        fdf[fdf["Time"] == '145'], fdf[fdf["Time"] == '190']
times = [t45, t95, t145, t190]

for time in times:
    time.sort_values(by=['TotalE'], inplace=True)
    ranking = [j+1 for j in range(len(time))]
    time['Ranking'] = ranking

heatmapper(times, upperDirectory, 'average')

# Save results
result = pd.concat(dataframeList)
# result.to_csv('output1.csv')
