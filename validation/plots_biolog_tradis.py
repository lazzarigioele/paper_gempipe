import numpy as np
np.random.seed(0)
import pandas as pnd
import seaborn as sbn
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Patch






def plot_predictions(match_df, tool_colors, focus_on='substrate', exclude=[], alltools=False, legend=False):
    match_df = match_df.copy()
    
    
    # exclude tools:
    match_df = match_df[match_df['tool'].isin(exclude)==False]
    
    match_dfs = {
        'match': match_df,
        'match_nounk': match_df[match_df['call'] != 'unknown'].copy()  # remove unknowns:
    }
        
        
    # define main theme
    sbn.set_theme(style="whitegrid")
        
    
    # create subplots according to the datasets
    datasets = sorted(match_df['dataset'].unique())
    fig, axes = plt.subplots(1, len(datasets), figsize=(3.5 * len(datasets), 5), dpi=300, sharey=False)
    if len(datasets) == 1: axes = [axes]  # ensure iterable if only one dataset / subplot


    for ax, dataset in zip(axes, datasets):
        
        for label, df in match_dfs.items():

            df_sub = df[df["dataset"] == dataset]
            
            
            count = df_sub.groupby(["strain", "tool"])[focus_on].nunique().reset_index(name="count")
            if alltools: 
                # eg show space for "pan-" tools also when "pan-" tools where not used in the dataset
                count['new_index'] = count['strain'] + count['tool'] 
                count = count.set_index('new_index', drop=True)
                strains = df_sub['strain'].unique() 
                tools = match_df['tool'].unique()
                missing_rows = []
                for strain in strains:
                    for tool in tools:
                        if strain + tool not in count.index:
                            missing_rows.append({
                                'future_index': strain + tool, 
                                'strain': strain,
                                'tool': tool,
                                'count': 0,
                            })
                missing_rows = pnd.DataFrame.from_records(missing_rows)
                if len(missing_rows) > 0:
                    missing_rows = missing_rows.set_index('future_index', drop=True, verify_integrity=True)
                    count = pnd.concat([count, missing_rows])
                count = count.reset_index(drop=True)
            summary = count.groupby(["tool"])["count"].agg(["mean", "std"]).reset_index()


            # ensure consistent tool order
            summary = summary.set_index("tool").reindex(tool_colors.keys()).dropna().reset_index()

            
            
        # ---- PLOT ----
            # bar plot (means with error bars)
            bars = ax.bar(
                summary["tool"],
                summary["mean"],
                #yerr = summary_sub["std"] if label == 'match' else None,
                color = ['lightgrey' for t in summary["tool"]] if label == 'match' else [tool_colors[t] for t in summary["tool"]],
                capsize = 5
            )  

            # overlay individual strain points with jitter
            for i, tool in enumerate(summary["tool"]):

                # Avoid '0' points consequences of 'alltools'==True:
                if tool not in match_df[match_df['dataset']==dataset]['tool'].unique():
                    continue
                    
                    
                # get y values:
                tool_points = count[count["tool"] == tool]["count"].values

                # get x values (jitter). 
                # loc= center of the distribution.
                # scale= std deviation of the normal distribution
                # size= how many random points to generate
                jitter_x = np.random.normal(loc=i, scale=0.1, size=len(tool_points)) 
                ax.scatter(
                    jitter_x,
                    tool_points,
                    color = "grey" if label == 'match' else 'black',
                    s = 10,  # circle dimension
                    alpha = 0.7,
                    zorder = 1 if label == 'match' else 2  # layer
                )

        ax.set_title(dataset, fontsize = 12)
        ax.set_xticks(range(len(summary["tool"])))
        ax.set_xticklabels(summary["tool"], rotation=45, ha="right")
        #ax.spines["top"].set_visible(False)
        #ax.spines["right"].set_visible(False)
        
        
    # add the legend
    if legend:
        patches = []
        tools_shown = [t for t in tool_colors.keys() if t in match_df['tool'].unique()]  # sort tools
        for tool in tools_shown:
            patches.append(Patch(color=tool_colors[tool], label=tool))    
        legend = plt.legend(
            handles=patches,
            bbox_to_anchor=(1.05, 0.5), loc='center left', frameon=False, title='tool')
    

    fig.tight_layout()
    return fig
    
    

def plot_predictions_match(match_df, tool_colors, focus_on='substrate', perc=False, common=False, std=True, exclude=[], alltools=False):
    match_df = match_df.copy()  # do not alter the original dataframe
    
    
    # exclude tools:
    match_df = match_df[match_df['tool'].isin(exclude)==False]
    # remove unknown:
    match_df = match_df[match_df['call'] != 'unknown']
    
    
    # define main theme
    sbn.set_theme(style="whitegrid")
        
            
    # define colors
    match_colors = {"TP": "yellowgreen", "TN": "forestgreen", "FP": "violet", "FN": "blueviolet", 'infeasible': 'grey'}

    
    # create subplots according to the datasets
    datasets = sorted(match_df['dataset'].unique())
    fig, axes = plt.subplots(1, len(datasets), figsize=(3.5 * len(datasets), 5), dpi=300, sharey=True if perc else False)
    if len(datasets) == 1: axes = [axes]  # ensure iterable if only one dataset / subplot

    
    # iterate subplots / datasets
    for ax, dataset in zip(axes, datasets):
        match_df_sub = match_df[match_df["dataset"] == dataset]


        # work with gids detected by all tools:
        if common:
            tools = match_df_sub["tool"].unique()
            # get the set of gids per tool
            gid_sets = {tool: set(match_df_sub.loc[match_df_sub["tool"] == tool, focus_on]) for tool in tools}
            # take intersection across all sets
            common_gids = set.intersection(*gid_sets.values())
            match_df_sub = match_df_sub[match_df_sub[focus_on].isin(common_gids)]
            

        # get count and summary tables
        count = match_df_sub.groupby(["strain", "tool"])[focus_on].nunique().reset_index(name="count")
        if alltools:  
            # eg show space for "pan-" tools also when "pan-" tools where not used in the dataset
            count['new_index'] = count['strain'] + count['tool'] 
            count = count.set_index('new_index', drop=True)
            strains = match_df_sub['strain'].unique() 
            tools = match_df['tool'].unique()
            missing_rows = []
            for strain in strains:
                for tool in tools:
                    if strain + tool not in count.index:
                        missing_rows.append({
                            'future_index': strain + tool, 
                            'strain': strain,
                            'tool': tool,
                            'count': 0,
                        })
            missing_rows = pnd.DataFrame.from_records(missing_rows)
            if len(missing_rows) > 0:
                missing_rows = missing_rows.set_index('future_index', drop=True, verify_integrity=True)
                count = pnd.concat([count, missing_rows])
            count = count.reset_index(drop=True)
        # finally create 'summary'
        summary = count.groupby(["tool"])["count"].agg(["mean", "std"]).reset_index()
        
        
        # get count and summary tables (call)
        count_call = match_df_sub.groupby(["strain", "tool", "call"])[focus_on].nunique().reset_index(name="count")
        # some strain in some tool might not have eg FN. If not corrected, this can alter mean and std comps.
        # add missing rows (0 rows)
        count_call['new_index'] = count_call['strain'] + count_call['tool'] + count_call['call']
        count_call = count_call.set_index('new_index', drop=True)
        strains = count_call['strain'].unique() 
        tools = count_call['tool'].unique() 
        if alltools:  # eg show space for "pan-" tools also when "pan-" tools where not used in the dataset
            tools = match_df['tool'].unique()
        calls = count_call['call'].unique() 
        missing_rows = []
        for strain in strains:
            for tool in tools:
                for call in calls:
                    if strain + tool + call not in count_call.index:
                        missing_rows.append({
                            'future_index': strain + tool + call, 
                            'strain': strain,
                            'tool': tool,
                            'call': call,
                            'count': 0,
                        })
        missing_rows = pnd.DataFrame.from_records(missing_rows)
        if len(missing_rows) > 0:
            missing_rows = missing_rows.set_index('future_index', drop=True, verify_integrity=True)
            count_call = pnd.concat([count_call, missing_rows])
        count_call = count_call.reset_index(drop=True)
        # finally create 'summary_call'
        summary_call = count_call.groupby(["tool", "call"])["count"].agg(["mean", "std"]).reset_index()
        

        # express 'summary_call' in percentages if requested
        if perc:
            for index, row in summary_call.iterrows(): 
                summary_sub = summary[summary['tool']==row['tool']].iloc[0]
                
                if summary_call.loc[index, 'mean'] != 0:
                    summary_call.loc[index, 'mean'] = summary_call.loc[index, 'mean'] / summary_sub['mean'] * 100
                    summary_call.loc[index, 'std'] = summary_call.loc[index, 'std'] / summary_sub['mean'] * 100
                        
        
        # get pivot
        pivots = {}
        for metric in ['mean', 'std']:
            pivots[metric] = summary_call.pivot(index="tool", columns="call", values=metric)
            pivots[metric] = pivots[metric].loc[[t for t in tool_colors.keys() if t in summary_call['tool'].unique()], ]  # correct row order
            pivots[metric] = pivots[metric][[m for m in match_colors.keys() if m in summary_call['call'].unique()]]  # correct col order
            pivots[metric] = pivots[metric].fillna(0)
                        
        
        
        # ---- PLOT ----
        bottom = 0
        for call in pivots['mean'].columns:

            # draw stacked bars
            ax.bar(
                pivots['mean'].index,
                pivots['mean'][call],
                bottom=bottom,
                label=call,
                color=match_colors[call]
            )
            
            # draw errors:
            if std:
                ax.errorbar(
                    pivots['mean'].index,
                    bottom + pivots['mean'][call],
                    yerr = pivots['std'][call],
                    fmt = "none",
                    ecolor = "black",
                    capsize = 3
                )
            
            # update bottom for stacking
            bottom = bottom + pivots['mean'][call]
            
    
        ax.set_title(dataset, fontsize = 12)
        ax.set_xticks(range(len(pivots['mean'].index)))
        ax.set_xticklabels(pivots['mean'].index, rotation=45, ha="right")
        #ax.spines["top"].set_visible(False)
        #ax.spines["right"].set_visible(False)
        if perc and dataset[:2]=='01': ax.set_ylabel("percentage")
    
    
    # add the legend
    patches = []
    match_shown = match_df['call'].unique()
    match_shown = reversed([m for m in match_colors.keys() if m in match_shown])  # sort tools
    for match in match_shown:
        patches.append(Patch(color=match_colors[match], label=match))    
    legend = plt.legend(
        handles=patches,
        bbox_to_anchor=(1.05, 0.5), loc='center left', frameon=False, title='metric')
    
    
    fig.tight_layout()
    return fig
    
    
    
def plot_predictions_metrics(match_df, tool_colors, focus_on='substrate', common=False, std=True, exclude=[], alltools=False):
    match_df = match_df.copy()  # do not alter the original dataframe
    
    
    # exclude tools:
    match_df = match_df[match_df['tool'].isin(exclude)==False]
    # remove unknown:
    match_df = match_df[match_df['call'] != 'unknown']
    
    
    # define main theme
    sbn.set_theme(style="whitegrid") 
    
    
    # create subplots according to the datasets
    datasets = sorted(match_df['dataset'].unique())
    fig, axes = plt.subplots(1, len(datasets), figsize=(5 * len(datasets), 4), dpi=300, sharey=True)
    if len(datasets) == 1: axes = [axes]  # ensure iterable if only one dataset / subplot
    

    # iterate subplots / datasets
    for ax, dataset in zip(axes, datasets):
        match_df_sub = match_df[match_df["dataset"] == dataset]
        
        
        # work with gids detected by all tools:
        if common:
            tools = match_df_sub["tool"].unique()
            # get the set of gids per tool
            gid_sets = {tool: set(match_df_sub.loc[match_df_sub["tool"] == tool, focus_on]) for tool in tools}
            # take intersection across all sets
            common_gids = set.intersection(*gid_sets.values())
            match_df_sub = match_df_sub[match_df_sub[focus_on].isin(common_gids)]
            
            
        # get count and summary tables (call)
        count_call = match_df_sub.groupby(["strain", "tool", "call"])[focus_on].nunique().reset_index(name="count")
        # some strain in some tool might not have eg FN. If not corrected, this can alter mean and std comps.
        # add missing rows (0 rows)
        count_call['new_index'] = count_call['strain'] + count_call['tool'] + count_call['call']
        count_call = count_call.set_index('new_index', drop=True)
        strains = count_call['strain'].unique() 
        tools = count_call['tool'].unique() 
        if alltools:  # eg show space for "pan-" tools also when "pan-" tools where not used in the dataset
            tools = match_df['tool'].unique()
        calls = count_call['call'].unique() 
        missing_rows = []
        for strain in strains:
            for tool in tools:
                for call in calls:
                    if strain + tool + call not in count_call.index:
                        missing_rows.append({
                            'future_index': strain + tool + call, 
                            'strain': strain,
                            'tool': tool,
                            'call': call,
                            'count': 0,
                        })
        missing_rows = pnd.DataFrame.from_records(missing_rows)
        if len(missing_rows) > 0:
            missing_rows = missing_rows.set_index('future_index', drop=True, verify_integrity=True)
            count_call = pnd.concat([count_call, missing_rows])
        count_call = count_call.reset_index(drop=True)

        
        # 'pivoted' is a dataframe
        # rows : strain-tool PAIRS    
        # columns : strain ,tool, FN, FP, TN, TP
        pivoted = count_call.pivot_table(
            index=["strain", "tool"],
            columns="call",
            values="count"
        ).reset_index()


        # compute metrics
        TP, TN, FP, FN = pivoted["TP"], pivoted["TN"], pivoted["FP"], pivoted["FN"]
        infeasible = pivoted["infeasible"] if 'infeasible' in pivoted.columns else 0
        factor = (TP + FP + TN + FN ) / (TP + FP + TN + FN + infeasible)
        pivoted["precision"] = TP / (TP + FP) * 100 * factor
        pivoted["recall"] = TP / (TP + FN) * 100 * factor
        pivoted["specificity"] = TN / (TN + FP) * 100 * factor
        pivoted["accuracy"] = (TP + TN) / (TP + TN + FP + FN) * 100 * factor
        
        # when means.loc[tool].to_numpy() produces all-NaN values, the bar has no valid height to plot.
        for m in ["precision", "recall", "specificity", "accuracy"]:
            pivoted[m] = pivoted[m].replace([np.inf, -np.inf], np.nan).fillna(0)
        

        # group by tool and compute mean + std of metrics across strains
        metrics = ["precision", "recall", "specificity", "accuracy"]
        means = pivoted.groupby("tool")[metrics].mean()
        stds = pivoted.groupby("tool")[metrics].std().fillna(0)

        # order rows (tools) according to tool_colors
        tools_ordered = [k for k in tool_colors.keys() if k in means.index]
        means = means.loc[tools_ordered]
        stds = stds.loc[tools_ordered]

        
        
        # ---- PLOT ----
        x = np.arange(len(metrics))  # one group per metric
        width = 0.8 / len(tools_ordered)  # divide bar width among tools

        for i, tool in enumerate(tools_ordered):

            ax.bar(
                x + i * width,
                means.loc[tool].to_numpy(),
                width,
                yerr = stds.loc[tool].to_numpy() if std else None,
                label = tool,
                color = tool_colors[tool],
                capsize = 0,
            )

        ax.set_xticks(x + width * (len(tools_ordered) - 1) / 2)
        ax.set_xticklabels(metrics, rotation=0, ha="center")
        ax.set_title(dataset)
        if dataset[:2]=='01': ax.set_ylabel("percentage")
        
        
    # add the legend  
    patches = []
    tools_shown = [t for t in tool_colors.keys() if t in match_df['tool'].unique()]  # sort tools
    for tool in tools_shown:
        patches.append(Patch(color=tool_colors[tool], label=tool))    
    legend = plt.legend(
        handles=patches,
        bbox_to_anchor=(1.05, 0.5), loc='center left', frameon=False, title='tool')
        
    
    fig.tight_layout()
    return fig
    
    
    
def plot_predictions_match2(match_df, tool_colors, focus_on='substrate', common=False, exclude=[], alltools=False, legend=False):
    match_df = match_df.copy()  # do not alter the original dataframe
    
    
    # exclude tools:
    match_df = match_df[match_df['tool'].isin(exclude)==False]
    # remove unknown:
    match_df = match_df[match_df['call'] != 'unknown']
    
    
    # define main theme
    sbn.set_theme(style="whitegrid")
        
            
    # define colors
    colors = {"TP": "yellowgreen", "TN": "forestgreen", "FP": "violet", "FN": "blueviolet", 'infeasible': 'grey'}

    
    # create subplots according to the datasets
    datasets = sorted(match_df['dataset'].unique())
    calls = [c for c in colors.keys() if c in match_df['call'].unique()]
    fig, axes = plt.subplots(len(calls), len(datasets), figsize=(4 * len(datasets), 3 * len(calls)), dpi=300)
    

    
    # iterate subplots / datasets
    for c, dataset in enumerate(datasets):
        match_df_sub = match_df[match_df["dataset"] == dataset]


        # work with gids detected by all tools:
        if common:
            tools = match_df_sub["tool"].unique()
            # get the set of gids per tool
            gid_sets = {tool: set(match_df_sub.loc[match_df_sub["tool"] == tool, focus_on]) for tool in tools}
            # take intersection across all sets
            common_gids = set.intersection(*gid_sets.values())
            match_df_sub = match_df_sub[match_df_sub[focus_on].isin(common_gids)]
            
    
        
        # get count and summary tables (call)
        count_call = match_df_sub.groupby(["strain", "tool", "call"])[focus_on].nunique().reset_index(name="count")
        # some strain in some tool might not have eg FN. If not corrected, this can alter mean and std comps.
        # add missing rows (0 rows)
        count_call['new_index'] = count_call['strain'] + count_call['tool'] + count_call['call']
        count_call = count_call.set_index('new_index', drop=True)
        strains = count_call['strain'].unique() 
        tools = count_call['tool'].unique() 
        if alltools:  # eg show space for "pan-" tools also when "pan-" tools where not used in the dataset
            tools = match_df['tool'].unique()
        calls = count_call['call'].unique() 
        missing_rows = []
        for strain in strains:
            for tool in tools:
                for call in calls:
                    if strain + tool + call not in count_call.index:
                        missing_rows.append({
                            'future_index': strain + tool + call, 
                            'strain': strain,
                            'tool': tool,
                            'call': call,
                            'count': 0,
                        })
        missing_rows = pnd.DataFrame.from_records(missing_rows)
        if len(missing_rows) > 0:
            missing_rows = missing_rows.set_index('future_index', drop=True, verify_integrity=True)
            count_call = pnd.concat([count_call, missing_rows])
        count_call = count_call.reset_index(drop=True)
        

        
        # ---- PLOT ----
        tool_order = [t for t in tool_colors.keys() if t in count_call['tool'].unique()]
        calls = [c for c in colors.keys() if c in match_df['call'].unique()]
        for r, call in enumerate(calls):
            count_call_sub = count_call[count_call["call"] == call]
            
            
            # draw boxplots
            sbn.boxplot(
                data=count_call_sub,
                x="tool",
                y="count",
                hue="tool",
                order=tool_order,
                palette=tool_colors,
                dodge=False,
                ax=axes[r, c],
                fliersize=0,  # size of the outlier markers (“fliers”) in the boxplot (0 because of the swarmplot, see below)
                legend=False
            )
            
            
            # Overlay swarmplot to show each strain as a point.
            # Avoid '0' points consequences of 'alltools'==True
            original_tools = match_df[(match_df["dataset"] == dataset) & (match_df["call"] == call)]['tool'].unique()
            count_call_sub = count_call_sub[count_call_sub['tool'].isin(original_tools)]
            sbn.stripplot(
                data=count_call_sub,
                x="tool",
                y="count",
                ax=axes[r, c],
                color="white",
                alpha=0.5,
                linewidth=0.5,
                edgecolor='black',
                size=3,
                jitter=0.3,  # increase spread (0.4 ≈ 40% of category width)
            )
            
            
            if r==0: axes[r, c].set_title(dataset, fontsize = 12)
            axes[r, c].set_xlabel('')
            if c==0: axes[r, c].set_ylabel(call, fontsize=12)
            else:
                axes[r, c].set_ylabel('')
            if r == len(calls)-1:
                axes[r, c].set_xticks(range(len(tool_order)))
                axes[r, c].set_xticklabels(tool_order, rotation=45, ha="right")
            else:
                axes[r, c].set_xticks(range(len(tool_order)))
                axes[r, c].set_xticklabels([])
                
            # ensure y-axis shows only integers
            axes[r, c].yaxis.set_major_locator(MaxNLocator(integer=True))
            
            #ax.spines["top"].set_visible(False)
            #ax.spines["right"].set_visible(False)
            
            
    
    fig.tight_layout()
    return fig
  