#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import plotly.graph_objects as go
import matplotlib.cm as cm
import matplotlib.colors as mcolors


# In[2]:


###global variable###


# In[3]:


df1=pd.read_csv("D:/pythonfile/all_csv/msMG-rbMG-newcelltype-epi-MappingTable.csv",index_col=0)
df2=pd.read_csv("D:/pythonfile/all_csv/msMG-sgMG-newcelltype-epi-MappingTable.csv",index_col=0)
dataset = 'newcelltype'


# In[4]:


df2


# In[5]:


def do_result(df1,sp1):
    df1.columns
    df_index=df1.index.tolist()
    df_column=df1.columns.tolist()
    data=df1.values
    result = [
        {"source":df_index[i],"target":df_column[j],"value":data[i,j]}
        for i in range(len(df_index))
        for j in range(len(df_column))
        if data[i,j]>0.1
    ]
    clean_result=[]
    for i in range(len(result)):
        if result[i]['source'].startswith(sp1):
            clean_result.append(result[i])
        else:
            continue
    return clean_result


# In[6]:


rbms=do_result(df1=df1,sp1='rb')
mssg=do_result(df1=df2,sp1='ms')


# In[7]:


mssg


# In[8]:


rbms


# In[9]:


nodes=list(sorted(set(d['source'] for d in rbms) | set(d['target'] for d in rbms) | set(d['target'] for d in mssg) | set(d['source'] for d in mssg)))


# In[10]:


nodes


# In[11]:


nodes_mapping={nodes[i]:i for i in range(len(nodes))}


# In[12]:


nodes_mapping


# In[13]:


rbms_links=[    
    {"source":nodes_mapping[d['source']],
    "target":nodes_mapping[d['target']],
    "value":d['value']}for d in rbms
]


# In[14]:


mssg_links=[    
    {"source":nodes_mapping[d['source']],
    "target":nodes_mapping[d['target']],
    "value":d['value']}for d in mssg
]


# In[15]:


rbms_links+mssg_links


# In[16]:


links=rbms_links+mssg_links


# In[17]:


links


# In[18]:


sankey_links = {
    'source':[d['source'] for d in links],
    'target':[d['target'] for d in links],
    'value':[d['value'] for d in links],   
}


# In[19]:


sankey_links


# In[20]:


fix_nodes=[a[3:] for a in nodes]


# In[21]:


fix_nodes


# In[22]:


ls1,ls2,ls3=[],[],[]
for point,node in enumerate(nodes_mapping):
    if node.startswith("rb"):
        ls1.append(point)
    if node.startswith("ms"):
        ls2.append(point)
    if node.startswith("sg"):
        ls3.append(point)
layers=dict(rb=ls1,
           ms=ls2,
           sg=ls3)


# In[34]:


layers


# In[27]:


node_colors=[None]*(max(max(indices) for indices in layers.values()) + 1)


# In[29]:


len(node_colors)


# In[32]:


def assign_simple_colors(layers, color_palette=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']):
    """为每个层（物种）分配简约的柔和颜色"""
    node_colors = [None] * (max(max(indices) for indices in layers.values()) + 1)
    for layer, indices in layers.items():
        # 按当前层分配柔和配色
        layer_colors = color_palette[:len(indices)]  # 使用前 len(indices) 个颜色
        for idx, color in zip(indices, layer_colors):
            node_colors[idx] = color
    return node_colors
node_colors = assign_simple_colors(layers)
node_colors


# In[33]:


link_colors = [
    f"rgba({int(c[1:3], 16)}, {int(c[3:5], 16)}, {int(c[5:7], 16)}, 0.5)"
    for c in [node_colors[source] for source in sankey_links["source"]]
]

# 修正和优化布局
fig = go.Figure(go.Sankey(
    node=dict(
        pad=30,  # 调整间距，使节点间的布局更美观
        thickness=25,  # 增加节点厚度，使视觉效果更协调
        line=dict(color="black", width=0.8),  # 调整节点边框
        label=fix_nodes,  # 节点标签
        color=node_colors,  # 使用已分配的节点颜色
    ),
    link=dict(
        source=sankey_links["source"],  # 链接起点
        target=sankey_links["target"],  # 链接终点
        value=sankey_links["value"],  # 权重
        color=link_colors,  # 链接颜色透明化
    )
))

# 确保每一列数据的 y 轴长度一致
fig.update_traces(arrangement="snap")  # 使用 snap 模式使布局整齐

# 更新布局
fig.update_layout(
    title_text=f"Cross-species MG SAMap Similarity Analysis {dataset}",  # 英文标题
    title_font=dict(size=24, color="darkblue"),  # 标题字体大小和颜色
    font_size=14,  # 全局字体大小
    paper_bgcolor="white",  # 图表背景颜色
    plot_bgcolor="white",  # 绘图区背景颜色
    margin=dict(l=60, r=60, t=60, b=60),  # 调整边距
    height=600,  # 调整图表高度
    width=900   # 调整图表宽度
)


# In[26]:


fig.write_html(f"sankey_{dataset}-3MG.html")


# In[27]:


##########################AG SANKEY PLOT####################################


# In[28]:


df3=pd.read_csv("D:/pythonfile/all_csv/rbAG-sgAG-newcelltype-epi-MappingTable.csv",index_col=0)


# In[29]:


rbsg=do_result(df1=df3,sp1='rb')


# In[30]:


rbsg


# In[31]:


nodes=list(sorted(set(d['source'] for d in rbsg) | set(d['target'] for d in rbsg)))
nodes


# In[32]:


nodes_mapping={nodes[i]:i for i in range(len(nodes))}
nodes_mapping


# In[33]:


rbsg_links=[    
    {"source":nodes_mapping[d['source']],
    "target":nodes_mapping[d['target']],
    "value":d['value']}for d in rbsg
]


# In[34]:


sankey_links = {
    'source':[d['source'] for d in rbsg_links],
    'target':[d['target'] for d in rbsg_links],
    'value':[d['value'] for d in rbsg_links],   
}
sankey_links


# In[35]:


fix_nodes=[a[3:] for a in nodes]
fix_nodes


# In[36]:


# 修正和优化布局
fig = go.Figure(go.Sankey(
    node=dict(
        pad=30,  # 调整间距，使节点间的布局更美观
        thickness=25,  # 增加节点厚度，使视觉效果更协调
        line=dict(color="black", width=0.8),  # 调整节点边框
        label=fix_nodes,  # 节点标签
    ),
    link=dict(
        source=sankey_links["source"],  # 链接起点
        target=sankey_links["target"],  # 链接终点
        value=sankey_links["value"],  # 权重
    )
))

# 确保每一列数据的 y 轴长度一致
fig.update_traces(arrangement="snap")  # 使用 snap 模式使布局整齐

# 更新布局
fig.update_layout(
    title_text="Cross-species AG SAMap Similarity Analysis",  # 英文标题
    title_font=dict(size=24, color="darkblue"),  # 标题字体大小和颜色
    font_size=14,  # 全局字体大小
    paper_bgcolor="white",  # 图表背景颜色
    plot_bgcolor="white",  # 绘图区背景颜色
    margin=dict(l=60, r=60, t=60, b=60),  # 调整边距
    height=600,  # 调整图表高度
    width=900   # 调整图表宽度
)
fig.write_html("sankey_diagram_interactive-AG.html")


# In[37]:


##########STAGE SANKEYPLOT##########################


# In[ ]:





# In[38]:


import pandas as pd
import plotly.graph_objects as go
import matplotlib.cm as cm
import matplotlib.colors as mcolors


# In[39]:


###global variable###


# In[42]:


df1=pd.read_csv("D:/pythonfile/all_csv/msMG-rbMG-stage-Epi-MappingTable.csv",index_col=0)
df2=pd.read_csv("D:/pythonfile/all_csv/msMG-sgMG-stage-Epi-MappingTable.csv",index_col=0)


# In[43]:


df2


# In[44]:


def do_result(df1,sp1):
    df1.columns
    df_index=df1.index.tolist()
    df_column=df1.columns.tolist()
    data=df1.values
    result = [
        {"source":df_index[i],"target":df_column[j],"value":data[i,j]}
        for i in range(len(df_index))
        for j in range(len(df_column))
        if data[i,j]>0.1
    ]
    clean_result=[]
    for i in range(len(result)):
        if result[i]['source'].startswith(sp1):
            clean_result.append(result[i])
        else:
            continue
    return clean_result


# In[45]:


rbms=do_result(df1=df1,sp1='rb')
mssg=do_result(df1=df2,sp1='ms')


# In[46]:


mssg


# In[47]:


rbms


# In[48]:


nodes=list(sorted(set(d['source'] for d in rbms) | set(d['target'] for d in rbms) | set(d['target'] for d in mssg) | set(d['source'] for d in mssg)))


# In[49]:


nodes


# In[50]:


nodes_mapping={nodes[i]:i for i in range(len(nodes))}


# In[51]:


nodes_mapping


# In[52]:


rbms_links=[    
    {"source":nodes_mapping[d['source']],
    "target":nodes_mapping[d['target']],
    "value":d['value']}for d in rbms
]


# In[53]:


mssg_links=[    
    {"source":nodes_mapping[d['source']],
    "target":nodes_mapping[d['target']],
    "value":d['value']}for d in mssg
]


# In[54]:


rbms_links+mssg_links


# In[55]:


links=rbms_links+mssg_links


# In[56]:


links


# In[57]:


sankey_links = {
    'source':[d['source'] for d in links],
    'target':[d['target'] for d in links],
    'value':[d['value'] for d in links],   
}


# In[58]:


sankey_links


# In[59]:


fix_nodes=[a[3:] for a in nodes]


# In[60]:


fix_nodes


# In[61]:


ls1,ls2,ls3=[],[],[]
for point,node in enumerate(nodes_mapping):
    if node.startswith("rb"):
        ls1.append(point)
    if node.startswith("ms"):
        ls2.append(point)
    if node.startswith("sg"):
        ls3.append(point)
layers=dict(rb=ls1,
           ms=ls2,
           sg=ls3)


# In[62]:


node_colors=[None]*(max(max(indices) for indices in layers.values()) + 1)


# In[63]:


def assign_simple_colors(layers, color_palette=['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']):
    """为每个层（物种）分配简约的柔和颜色"""
    node_colors = [None] * (max(max(indices) for indices in layers.values()) + 1)
    for layer, indices in layers.items():
        # 按当前层分配柔和配色
        layer_colors = color_palette[:len(indices)]  # 使用前 len(indices) 个颜色
        for idx, color in zip(indices, layer_colors):
            node_colors[idx] = color
    return node_colors
node_colors = assign_simple_colors(layers)
node_colors


# In[64]:


link_colors = [
    f"rgba({int(c[1:3], 16)}, {int(c[3:5], 16)}, {int(c[5:7], 16)}, 0.5)"
    for c in [node_colors[source] for source in sankey_links["source"]]
]

# 修正和优化布局
fig = go.Figure(go.Sankey(
    node=dict(
        pad=30,  # 调整间距，使节点间的布局更美观
        thickness=25,  # 增加节点厚度，使视觉效果更协调
        line=dict(color="black", width=0.8),  # 调整节点边框
        label=fix_nodes,  # 节点标签
        color=node_colors,  # 使用已分配的节点颜色
    ),
    link=dict(
        source=sankey_links["source"],  # 链接起点
        target=sankey_links["target"],  # 链接终点
        value=sankey_links["value"],  # 权重
        color=link_colors,  # 链接颜色透明化
    )
))

# 确保每一列数据的 y 轴长度一致
fig.update_traces(arrangement="snap")  # 使用 snap 模式使布局整齐

# 更新布局
fig.update_layout(
    title_text="Cross-species stage-MG-Epi SAMap Similarity Analysis",  # 英文标题
    title_font=dict(size=24, color="darkblue"),  # 标题字体大小和颜色
    font_size=14,  # 全局字体大小
    paper_bgcolor="white",  # 图表背景颜色
    plot_bgcolor="white",  # 绘图区背景颜色
    margin=dict(l=60, r=60, t=60, b=60),  # 调整边距
    height=600,  # 调整图表高度
    width=900   # 调整图表宽度
)


# In[65]:


fig.write_html("sankey_diagram_interactive-3MG-stage-Epi.html")


# In[ ]:




