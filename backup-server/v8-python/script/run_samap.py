#!/usr/bin/env python
# coding: utf-8

# In[27]:


from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os


# In[ ]:


##variable###


# In[26]:


datasetlist=['M-MG','R-MG','S-MG','R-AG','R-CG','S-AG']
processed_pairs=set()
catgory="newcelltype"
epi='-Epi'


# In[ ]:


#####SAMAP#############


# In[ ]:


def do_samap(sp1,sp2,gd1,gd2,fn1,fn2):
    sam1=SAM()
    sam1.load_data(fn1)
    sam2=SAM()
    sam2.load_data(fn2) 
    sams = {sp1:sam1,sp2:sam2}
    sm = SAMAP(
            sams,
            f_maps = '/data01/sunxuebo/project/scrnaseq/v8-python/samap/re-maps/' 
        )
    sm.run(pairwise=False)
    samap = sm.samap
    
    keys = {sp1:catgory,sp2:catgory}
    D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
    D.head()
    D.to_csv(f"{sp1}{gd1}-{sp2}{gd2}-{catgory}{epi}-D.csv")
    
    MappingTable.head()
    MappingTable.to_csv(f"{sp1}{gd1}-{sp2}{gd2}-{catgory}{epi}-MappingTable.csv")
    
    #MappingTable.set_index(MappingTable.columns[0], inplace=True)
    
    # Filter the columns and rows based on the sp1 and sp2 prefixes
    filtered_columns = [col for col in MappingTable.columns if col.startswith(sp1)]
    filtered_rows = [row for row in MappingTable.index if row.startswith(sp2)]
    
    # Filter the MappingTable based on these selected columns and rows
    filtered_data = MappingTable.loc[filtered_rows, filtered_columns]
    
    # Set up figure dimensions
    fig = plt.figure(figsize=(12, 8))
    
    # Create a grid to allocate space for the colorbar on the right
    grid = plt.GridSpec(1, 2, width_ratios=[0.95, 0.05], wspace=0.3)
    
    # Add axes for the heatmap
    ax = fig.add_subplot(grid[0, 0])
    
    # Add axes for the colorbar on the right
    cbar_ax = fig.add_subplot(grid[0, 1])
    
    # Plot the heatmap
    sns.heatmap(
        filtered_data,
        cmap="viridis",         # Use viridis colormap
        annot=False,            # Disable cell annotations
        linewidths=0.5,         # Add grid lines for better separation
        linecolor="gray",       # Grid line color
        cbar=True,              # Enable color bar
        cbar_ax=cbar_ax,        # Place the color bar in the defined space
        cbar_kws={"label": "Cell Similarity Scores"},  # Label for the color bar
        ax=ax                   # Assign heatmap to the main axes
    )
    
    # Adjust colorbar orientation
    cbar_ax.yaxis.set_label_position("left")
    cbar_ax.yaxis.tick_left()
    
    # Rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    
    # Add title and axis labels
    ax.set_title(f"{sp1}{gd1}-{sp2}{gd2} Cell Similarity Scores Based on SAMap{epi}", fontsize=16)
    ax.set_xlabel(f"{catgory}", fontsize=12)
    ax.set_ylabel(f"{catgory}", fontsize=12)
    
    # Save the figure
    plt.savefig(f'{sp1}{gd1}-{sp2}{gd2}-{catgory}{epi}-samap-heatmap.png', dpi=300, bbox_inches='tight')
    
    # Show the plot
    plt.show()


# In[6]:


base_dir = f"./{catgory}"
if not os.path.exists(base_dir):
    os.makedirs(base_dir)
for dataset1 in datasetlist:
    for dataset2 in datasetlist:  
        if dataset1:
            ls=dataset1.split(sep='-')
        if ls[0]=='M':
            sp1='ms'
        elif ls[0]=='R':
            sp1='rb'
        elif ls[0]=='S':
            sp1='sg'
        gd1=ls[1]
        fn1 = f'./{dataset1}_counts_pr.h5ad'
        if dataset2:
            ls=dataset2.split(sep='-')
        if ls[0]=='M':
            sp2='ms'
        elif ls[0]=='R':
            sp2='rb'
        elif ls[0]=='S':
            sp2='sg'
        gd2=ls[1]
        fn2 = f'./{dataset2}_counts_pr.h5ad'  
        if sp1 == sp2:
            continue
        pair = tuple(sorted([dataset1, dataset2]))
        if pair in processed_pairs:
            continue
        do_samap(sp1,sp2,gd1,gd2,fn1,fn2)
        processed_pairs.add(pair)
        print(pair)
#MappingTable=pd.read_csv("D:/111/rbMG-sgMG-MappingTable.csv",index_col = 0)
#MappingTable


# In[ ]:


#MappingTable=pd.read_csv("D:/111/rbMG-sgMG-MappingTable.csv",index_col = 0)


# In[2]:


import pandas as pd
import numpy as np
import holoviews as hv
from holoviews import opts

def sankey_plot(M,species_order=None,align_thr=0.1,save_to='sankeyplot.png',title='sankeyplot',**params):
    """Generate a sankey plot
    
    Parameters
    ----------
    M: pandas.DataFrame
        Mapping table output from `get_mapping_scores` (second output).

    align_thr: float, optional, default 0.1
        The alignment score threshold below which to remove cell type mappings.
    
    species_order: list, optional, default None
        Specify the order of species (left-to-right) in the sankey plot.
        For example, `species_order=['hu','le','ms']`.

    Keyword arguments
    -----------------
    Keyword arguments will be passed to `sankey.opts`.
    """    
    def q(x):
        return np.array(list(x))
    if species_order is not None:
        ids = np.array(species_order)
    else:
        ids = np.unique([x.split('_')[0] for x in M.index])

    if len(ids)>2:
        d = M.values.copy()
        d[d<align_thr]=0
        x,y = d.nonzero()
        x,y = np.unique(np.sort(np.vstack((x,y)).T,axis=1),axis=0).T
        values = d[x,y]
        nodes = q(M.index)

        node_pairs = nodes[np.vstack((x,y)).T]
        sn1 = q([xi.split('_')[0] for xi in node_pairs[:,0]])
        sn2 = q([xi.split('_')[0] for xi in node_pairs[:,1]])
        filt = np.logical_or(
            np.logical_or(np.logical_and(sn1==ids[0],sn2==ids[1]),np.logical_and(sn1==ids[1],sn2==ids[0])),
            np.logical_or(np.logical_and(sn1==ids[1],sn2==ids[2]),np.logical_and(sn1==ids[2],sn2==ids[1]))
        )
        x,y,values=x[filt],y[filt],values[filt]
        
        d=dict(zip(ids,list(np.arange(len(ids)))))        
        depth_map = dict(zip(nodes,[d[xi.split('_')[0]] for xi in nodes]))
        data =  nodes[np.vstack((x,y))].T
        for i in range(data.shape[0]):
            if d[data[i,0].split('_')[0]] > d[data[i,1].split('_')[0]]:
                data[i,:]=data[i,::-1]
        R = pd.DataFrame(data = data,columns=['source','target'])
        
        R['Value'] = values       
    else:
        d = M.values.copy()
        d[d<align_thr]=0
        x,y = d.nonzero()
        x,y = np.unique(np.sort(np.vstack((x,y)).T,axis=1),axis=0).T
        values = d[x,y]
        nodes = q(M.index)
        R = pd.DataFrame(data = nodes[np.vstack((x,y))].T,columns=['source','target'])
        R['Value'] = values
        depth_map=None
    
    try:
        from holoviews import dim
        #from bokeh.models import Label
        import holoviews as hv
        hv.extension('bokeh',logo=False)
        hv.output(size=100)        
    except:
        raise ImportError('Please install holoviews-samap with `!pip install holoviews-samap`.')

    def f(plot,element):
        plot.handles['plot'].sizing_mode='fixed'    
        plot.handles['plot'].x_range.start = -600    
        #plot.handles['plot'].add_layout(Label(x=plot.handles['plot'].x_range.end*0.78, y=plot.handles['plot'].y_range.end*0.96, text=id2))
        plot.handles['plot'].x_range.end = 1500    
        #plot.handles['plot'].add_layout(Label(x=0, y=plot.handles['plot'].y_range.end*0.96, text=id1))


    sankey1 = hv.Sankey(R, kdims=["source", "target"])#, vdims=["Value"])

    cmap = params.get('cmap','Colorblind')
    label_position = params.get('label_position','outer')
    edge_line_width = params.get('edge_line_width',0)
    show_values = params.get('show_values',False)
    node_padding = params.get('node_padding',4)
    node_alpha = params.get('node_alpha',1.0)
    node_width = params.get('node_width',40)
    node_sort = params.get('node_sort',True)
    frame_height = params.get('frame_height',1000)
    frame_width = params.get('frame_width',800)
    bgcolor = params.get('bgcolor','snow')
    apply_ranges = params.get('apply_ranges',True)


    sankey1.opts(cmap=cmap,label_position=label_position, edge_line_width=edge_line_width, show_values=show_values,
                 node_padding=node_padding,node_cmap=depth_map, node_alpha=node_alpha, node_width=node_width,
                 node_sort=node_sort,frame_height=frame_height,frame_width=frame_width,bgcolor=bgcolor,
                 apply_ranges=apply_ranges,hooks=[f],title=title)


    # Save the plot if save_to is specified
    if save_to:
        # Save as PNG, SVG, PDF, etc.
        hv.save(sankey1, save_to)
        print(f"Sankey plot saved to {save_to}")
    return sankey1


# In[8]:


#MappingTable=pd.read_csv("D:/111/rbAG-sgAG-MappingTable.csv",index_col = 0)


# In[9]:


#sankey_plot(MappingTable, align_thr=0.1, species_order=[sp1,sp2],save_to=f'{sp1}{gd1}-{sp2}{gd2}-samap-sankey.png',title=f'{sp1}{gd1}-{sp2}{gd2}-Cell Similarity')


# In[ ]:




