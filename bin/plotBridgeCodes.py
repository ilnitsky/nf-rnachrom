import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import sys

def plot_bridge_codes(statistics_file, mode, input_path, output_path):

    SE_code = {
        "Code": 
            ["0","F-00","F-0R","F-D0","F-DR","R-00","R-0R","R-D0","R-DR","FR"],
        "Description": 
            ["bridge not found",
            "F-bridge found, RNA and DNA parts did not pass the length filter",
            "F-bridge found, only the RNA part passed the length filter",
            "F-bridge found, only the DNA part passed the length filter",
            "F-bridge found, RNA and DNA parts were filtered for length",
            "R-bridge found, RNA and DNA parts did not pass the length filter",
            "R-bridge found, only the RNA part passed the length filter",
            "R-bridge found, only the DNA part passed the length filter",
            "R-bridge found, RNA and DNA parts were filtered for length",
            "both forward and reverse bridges found"]
    }

    PE_code = {
        "Code": ["R-0-DR","0-R-DR","F-0-DR","0-F-DR","R-0-0R","R-0-D0","R-0-00","0-R-0R","0-R-D0","0-R-00","F-0-0R","F-0-D0","F-0-00","0-F-0R","0-F-D0","0-F-00","0-0","0-FR","FR-0","R-R","F-F","R-F","F-R","R-FR","F-FR","FR-R","FR-F","FR-FR"],
        "Description": 
            ["reverse (R) bridge found in R1, RNA and DNA parts were filtered for length",
            "R-bridge found in R2, RNA and DNA parts were filtered for length",
            "forward (F) bridge found in R1, RNA and DNA parts were filtered for length",
            "F-bridge found in R2, RNA and DNA parts were filtered for length",
            "R-bridge found in R1, only the RNA part passed the length filter",
            "R-bridge found in R1, only the DNA part passed the length filter",
            "R-bridge found in R1, RNA and DNA parts did not pass the length filter",
            "R-bridge found in R2, only the RNA part passed the length filter",
            "R-bridge found in R2, only the RNA part passed the length filter",
            "R-bridge found in R2, RNA and DNA parts did not pass the length filter",
            "F-bridge found in R1, only the RNA part passed the length filter",
            "F-bridge found in R1, only the DNA part passed the length filter",
            "F-bridge found in R1, RNA and DNA parts did not pass the length filter",
            "F-bridge found in R2, only the RNA part passed the length filter",
            "F-bridge found in R2, only the DNA part passed the length filter",
            "F-bridge found in R2, RNA and DNA parts did not pass the length filter",
            "bridge not found",
            "two bridges forward and reverse found in R2",
            "two bridges forward and reverse found in R1",
            "found reverse bridge in R1 and in R2",
            "found forward bridge in R1 and in R2",
            "found reverse bridge in R1 and forward bridge in R2",
            "found forward bridge in R1 and reverse bridge in R2",
            "found reverse bridge in R1, and also forward and reverse bridge in R2",
            "found forward bridge in R1, and also forward and reverse bridge in R2",
            "forward and reverse bridge found in R1, and reverse bridge in R2",
            "forward and reverse bridge found in R1, and forward bridge in R2",
            "forward and reverse bridges are found in both R1 and R2"] 
        }

    statistics = pd.read_csv(input_path + statistics_file, header=None, sep='\t', dtype = 'str')
    statistics = statistics.rename({0: 'Code', 1: 'N-reads'}, axis='columns')
    statistics['N-reads'] = statistics['N-reads'].apply(lambda x: int(x))
    
    if mode == "SE":
        code_df = pd.DataFrame.from_dict(SE_code)
    else:
        code_df = pd.DataFrame.from_dict(PE_code)
    result = pd.merge(statistics, code_df, on="Code", how="left")
    
    total = result['N-reads'].sum()
    
    plt.figure()
    ylim_max = result['N-reads'].max()*2
    
    g = sns.barplot(result, x="Code", y="N-reads", estimator="sum", errorbar=None)
    g.set(title='{N} reads'.format(N=total))

    g.bar_label(g.containers[0], fontsize=10)
    g.set(ylim=(0.9, ylim_max))
    plt.yscale('log')
    plt.xticks(rotation=90)
    
    text = ""
    for k in str(result[['Code','Description']].to_dict('records'))[2:-2].split('}, {'):
        text = text + k.replace("'Code': ","").replace(", 'Description'","").replace("'","") + "\n"
    plt.figtext(0, 1, text, ha="left", fontsize=12) #bbox={"facecolor":"orange", "alpha":0.5, "pad":5} 
    # plt.show()
    plt.savefig('{output_path}{statistics_file}_{mode}.png'.format(output_path=output_path, statistics_file=statistics_file, mode=mode), dpi=360, bbox_inches="tight")

plot_bridge_codes(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])