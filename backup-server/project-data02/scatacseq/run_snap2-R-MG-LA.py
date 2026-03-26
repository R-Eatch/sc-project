#!/usr/bin/env python
# coding: utf-8

# In[1]:


import snapatac2 as snap
import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
import os
import re
import gzip
from typing import Iterable


# In[ ]:


import os
import tempfile

# 设置自定义 TMPDIR
custom_tmp = "/data02/sunxuebo/project/scatacseq/tmp"
os.makedirs(custom_tmp, exist_ok=True)

# 修改临时目录环境变量和 Python 的 tempfile 缓存
os.environ["TMPDIR"] = custom_tmp
tempfile.tempdir = custom_tmp

# 可选：打印确认
print(f"[INFO] Using custom tmp dir: {tempfile.gettempdir()}")


# In[7]:


sc._settings.ScanpyConfig(n_jobs=1)


# In[2]:


dataset= 'R-MG-LA'
species = 'rabbit'
marker_genes = marker_genes = [
    # Luminal epithelial markers
    'EPCAM',    # Luminal epithelial marker
    'KRT8',     # Luminal epithelial keratin
    'KRT18',    # Luminal epithelial keratin
    'LALBA',    # Lactating luminal (secretory) cell marker
    'ESR1',     # Estrogen receptor, hormone-responsive luminal cell
    'PIP',      # Secretory luminal cell
    'ELF5',     # Luminal progenitor and alveolar differentiation regulator
    'PIGR',     # Luminal secretory function
    'XBP1',     # Luminal secretory (milk biosynthesis)
    'AZGP1',    # Mature luminal secretory
    'CD36',     # Lipid metabolism in secretory cells

    # Basal/myoepithelial markers
    'KRT5',     # Basal keratin
    'KRT14',    # Basal keratin
    'ACTA2',    # Myoepithelial (smooth muscle actin)
    'TP63',     # Basal/myoepithelial transcription factor
    'CDH3',     # Basal lineage marker
    'KRT17',    # Basal keratin, progenitor
    'MMP14',    # Basal ECM interaction

    # Fibroblast / stromal
    'COL1A1',   # Fibroblast / stromal
    'DCN',      # Decorin, fibroblast ECM
    'PDGFRA',   # Fibroblast/stromal
    'VIM',      # Mesenchymal marker
    'FAP',      # Activated fibroblast
    'FN1',      # ECM component

    # Endothelial cells
    'VWF',      # Endothelial marker
    'PECAM1',   # CD31, pan-endothelial
    'CDH5',     # VE-cadherin

    # Immune cells
    'PTPRC',    # CD45, pan-leukocyte
    'CD68',     # Macrophage
    'CD3D',     # T cell
    'CD79A',    # B cell
    'NKG7',     # NK / cytotoxic lymphoid cells
    'LYZ',      # Myeloid (monocytes/macrophages)

    # Adipocytes
    'FABP4',    # Adipocyte marker
    'PLIN1',    # Lipid droplet protein
    'ADIPOQ',   # Adiponectin

    # Differentiation/others
    'KRT1',     # Differentiated suprabasal keratinocyte (epidermal contamination)
    'KRT7',     # Ductal epithelium
    'KRT19',    # Luminal epithelial
]
n_neighbors=30
res=0.5
celltype_dict={
    'Luminal':[2],
    'Basal':[9],
    'Fibro':[4,5],
    'Endothelial':[10],
    'c7-Krt5+':[7],
    'Immune':[0,1,3],
    'c6':[6],
    'c8':[8]
}
UPDATA_CELLTYPE=False
output_dir=f'./{dataset}/results'
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir


# In[ ]:


frag_file   = Path(f"./data/{dataset}/outs/fragments.tsv.gz")   # 10X fragment
out_h5ad    = Path(f"./{dataset}/{dataset}_raw.h5ad")               # backed 模式
chrom_fa    = Path(f"./data/genome/{species}.fa")      
gene_anno = Path(f"./data/genome/{species}.gtf")


# In[42]:


def clean_gtf(
    in_gtf: str | Path,
    out_gtf: str | Path,
    keep_keys: Iterable[str] = ("gene_id", "gene_name", "transcript_id"),
) -> Path:
    """
    清洗 GTF/GTF.GZ 文件
    - 如果 out_gtf 已存在，直接跳过。
    - 仅保留 keep_keys 指定的属性。
    - 若 gene_name 缺失，则用 gene_id 填补。

    Parameters
    ----------
    in_gtf : str | Path
        输入 GTF 路径，可为 .gtf 或 .gtf.gz
    out_gtf : str | Path
        输出文件路径，扩展名可自由指定
    keep_keys : Iterable[str], optional
        需要保留的 attribute 名称

    Returns
    -------
    Path
        输出文件的 Path 对象
    """
    in_gtf, out_gtf = Path(in_gtf), Path(out_gtf)

    # 若输出已存在则跳过
    if out_gtf.exists():
        print(f"输出文件已存在：{out_gtf}，跳过清洗。")
        return out_gtf

    # 根据扩展名选择打开方式
    def _open(path: Path, mode: str = "rt"):
        return gzip.open(path, mode) if path.suffix == ".gz" else open(path, mode)

    # 正则
    re_attr = re.compile(r'(\w+)\s+"([^"]*)"')
    re_fix_semicolon = re.compile(r";(?!\s)")

    keep_keys = set(keep_keys)

    with _open(in_gtf, "rt") as fin, _open(out_gtf, "wt") as fout:
        for line in fin:
            # 注释行或空行
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                fout.write(line)
                continue

            attrs = cols[8]

            # 修正分号后缺少空格的情况
            attrs = re_fix_semicolon.sub("; ", attrs)

            # 提取 key-value
            kv_pairs = re_attr.findall(attrs)

            # 过滤并记录 gene_id
            clean_pairs, gene_id_val = [], None
            for k, v in kv_pairs:
                if k == "gene_id":
                    gene_id_val = v
                if k in keep_keys:
                    clean_pairs.append((k, v))

            # 若缺少 gene_name，用 gene_id 填补
            if gene_id_val and not any(k == "gene_name" for k, _ in clean_pairs):
                clean_pairs.append(("gene_name", gene_id_val))

            # gene_id 缺失则写回原行
            if gene_id_val is None:
                fout.write(line)
                continue

            # 重新拼装 attributes
            attrs_new = "; ".join(f'{k} "{v}"' for k, v in clean_pairs) + ";"
            cols[8] = attrs_new
            fout.write("\t".join(cols) + "\n")

    print(f"GTF 清洗完成，输出文件：{out_gtf}")
    return out_gtf


# In[56]:


import pandas as pd
import anndata as ad
from pathlib import Path
from typing import Union, Iterable

def add_gene_names_from_gtf(
    adata: ad.AnnData,
    gtf_path: Union[str, Path],
    gtf_feature: str = "gene",
    gene_id_key: str = "gene_id",
    gene_name_key: str = "gene_name",
    keep_attrs: Iterable[str] = ("gene_id", "gene_name"),
) -> ad.AnnData:
    """
    读取 GTF，提取 gene_id 与 gene_name 映射，写入 adata.var['gene_name']。
    var_names 必须是 gene_id。
    """
    gtf_path = Path(gtf_path)
    cols = [
        "seqname", "source", "feature",
        "start", "end", "score",
        "strand", "frame", "attribute",
    ]

    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=cols,
        low_memory=False,
    )
    gtf = gtf[gtf["feature"] == gtf_feature]

    # 解析第 9 列属性
    def parse_attr(attr: str) -> dict:
        out = {}
        for kv in attr.strip().split(";"):
            kv = kv.strip()
            if not kv:
                continue
            key, val = kv.split(" ", 1)
            out[key] = val.strip('"')
        return out

    attrs = (
        gtf["attribute"]
        .apply(parse_attr)
        .apply(pd.Series)
        .loc[:, list(keep_attrs)]
        .rename(columns={gene_id_key: "gene_id", gene_name_key: "gene_name"})
    )

    mapping = (
        attrs.drop_duplicates(subset="gene_id")
             .set_index("gene_id")
    )

    # 用索引填补缺失 gene_name
    missing = mapping["gene_name"].isna()
    if missing.any():
        mapping.loc[missing, "gene_name"] = mapping.index[missing]

    # 保证 var_names 是 str
    adata.var_names = adata.var_names.astype(str)

    # 合并进 .var
    adata.var = adata.var.join(mapping, how="left")

    # 再次检查 .var 中是否仍缺 gene_name
    mask = adata.var["gene_name"].isna()
    if mask.any():
        adata.var.loc[mask, "gene_name"] = adata.var.index[mask]

    return adata


# In[74]:


def invert_annotation_dict(cluster_to_cells: dict) -> dict:
    return {
        str(idx): celltype
        for celltype, idx_list in cluster_to_cells.items()
        for idx in idx_list
    }


# In[5]:


if not os.path.exists(f'./{dataset}'):
    os.mkdir(f'./{dataset}')
else:
    print(f'dir {dataset} existed')


# In[14]:


from Bio import SeqIO
chrom_sizes = {rec.id: len(rec.seq) 
               for rec in SeqIO.parse(chrom_fa, "fasta")}


# In[8]:


frag_file


# In[25]:


adata = snap.pp.import_fragments(
    frag_file,
    chrom_sizes=chrom_sizes,
    sorted_by_barcode=False  # 10X 输出默认按位置排
)
adata


# In[10]:


adata


# In[11]:


snap.pl.frag_size_distr(adata, interactive=False)


# In[86]:


fig = snap.pl.frag_size_distr(adata, show=False)
fig.update_yaxes(type="log")
fig.write_image(f"{output_dir}/{dataset}_frag_plot.png",scale=2)
fig.show()


# In[13]:


snap.metrics.tsse(adata=adata, gene_anno=gene_anno)


# In[82]:


tsse = snap.pl.tsse(adata, interactive=False,show=False, out_file=None)


# In[88]:


tsse.write_image(f"{output_dir}/{dataset}_tsse_plot.png", scale=2)
tsse.show()


# In[15]:


snap.pp.filter_cells(adata, min_counts=5000, min_tsse=5, max_counts=100000)
adata


# In[16]:


snap.pp.add_tile_matrix(adata)


# In[17]:


snap.pp.select_features(adata, n_features=250000)


# In[21]:


mask = adata.var['selected'].to_numpy(bool)  # 或者 .values.astype(bool)
snap.preprocessing.scrublet(
    adata,
    features=mask,
    n_comps=50,
    sim_doublet_ratio=2,
    expected_doublet_rate=0.05,
    inplace=True
)


# In[22]:


snap.pp.filter_doublets(adata)


# In[23]:


snap.tl.spectral(adata)


# In[24]:


snap.tl.umap(adata,random_state=42)


# In[26]:


snap.pp.knn(adata,n_neighbors=n_neighbors)
snap.tl.leiden(adata,resolution=res)


# In[27]:


snap.pl.umap(adata, color='leiden', interactive=False, height=500)


# In[30]:


adata


# In[29]:


adata.write(f'./{dataset}/{dataset}_cleaned.h5ad')


# In[3]:


#adata=snap.read('R-LA-MG-ATAC_cleaned.h5ad')


# In[31]:


gene_anno


# In[44]:


fixed_gtf = clean_gtf(gene_anno, f"./{species}.clean.gtf.gz")


# In[45]:


gene_matrix = snap.pp.make_gene_matrix(adata,gene_anno=fixed_gtf,gene_id_key='gene_id',gene_name_key='gene_name')
gene_matrix


# In[46]:


sc.pp.filter_genes(gene_matrix, min_cells= 5)
sc.pp.normalize_total(gene_matrix)
sc.pp.log1p(gene_matrix)


# In[47]:


sc.external.pp.magic(gene_matrix, solver="approximate")


# In[48]:


gene_matrix.obsm["X_umap"] = adata.obsm["X_umap"]


# In[57]:


gene_matrix = add_gene_names_from_gtf(
    adata=gene_matrix,
    gtf_path=fixed_gtf
)


# In[60]:


#gene_matrix.var


# In[71]:


gene_matrix.obs['cell_type']=gene_matrix.obs['leiden']


# In[ ]:


if UPDATA_CELLTYPE:
    anno_dict = invert_annotation_dict(celltype_dict)
    gene_matrix.obs['cell_type']=gene_matrix.obs['leiden'].map(anno_dict)
    adata.obs['cell_type'] = gene_matrix.obs.loc[adata.obs_names, 'cell_type'].to_numpy()
else:
    gene_matrix.obs['cell_type']=gene_matrix.obs['leiden']


# In[94]:


filtered_genes = [g for g in marker_genes if g in gene_matrix.var["gene_name"].values]
sc.pl.umap(gene_matrix, use_raw=False, color=["leiden","cell_type"] + filtered_genes,gene_symbols="gene_name",legend_loc='on data',save=f'{dataset}_fp.png')


# In[102]:


gene_matrix.var['gene_id']=gene_matrix.var_names


# In[103]:


gene_matrix.var_names=gene_matrix.var['gene_name']


# In[ ]:


#sc.tl.rank_genes_groups(gene_matrix, groupby="leiden", method="wilcoxon", use_raw=False, pts=True)
#deg_df = sc.get.rank_genes_groups_df(gene_matrix, group=None)
#deg_df.to_csv(f"{output_dir}/{dataset}_deg_by_leiden.csv", index=False)#15GB 内存爆炸


# In[105]:


# sc.pl.rank_genes_groups_dotplot(
#     gene_matrix,
#     n_genes=5,              # 每个聚类显示前 5 个 marker
#     show=False,             # 不在 notebook 显示
#     return_fig=True
# ).savefig(f"{output_dir}/{dataset}_dotplot_cluster.png", dpi=300)


# In[ ]:


gene_matrix.write(f"./{dataset}/{dataset}_gene_mat.h5ad", compression='gzip')


# In[3]:


adata = snap.read(f'./{dataset}/{dataset}_cleaned.h5ad', backed=None)


# In[4]:


gene_matrix=snap.read(f"./{dataset}/{dataset}_gene_mat.h5ad",backed=None)


# In[8]:


adata.obs['cell_type'] = gene_matrix.obs['cell_type']


# In[ ]:


##############################################################################


# In[9]:


snap.tl.macs3(adata, groupby='cell_type',n_jobs=1)


# In[16]:


peaks = snap.tl.merge_peaks(adata.uns['macs3'],chrom_sizes=chrom_sizes)
#peaks.head()


# In[17]:


peak_mat = snap.pp.make_peak_matrix(adata, use_rep=peaks['Peaks'])
#peak_mat


# In[18]:


marker_peaks = snap.tl.marker_regions(peak_mat, groupby='cell_type', pvalue=0.01)


# In[20]:


reg = snap.pl.regions(peak_mat, groupby='cell_type', peaks=marker_peaks, interactive=False,show=False)
reg.write_image(f'{output_dir}/{dataset}_reg.png')


# In[21]:


motifs = snap.tl.motif_enrichment(
    motifs=snap.datasets.cis_bp(unique=True),
    regions=marker_peaks,
    genome_fasta=chrom_fa,
)


# In[ ]:


motif_fig = snap.pl.motif_enrichment(
    motifs,
    max_fdr=0.0001,
    interactive=False,
    show=False
)


# In[ ]:


motif_fig.write_image(f'{output_dir}/{dataset}_motif.png',scale = 2)


# In[ ]:


# group1 = "Basal"
# group2 = "LumHR"
# naive_B = adata.obs['cell_type'] == group1
# memory_B = adata.obs['cell_type'] == group2
# peaks_selected = np.logical_or(
#     peaks[group1].to_numpy(),
#     peaks[group2].to_numpy(),
# )

# peaks_selected

# peak_mat

# %%time
# diff_peaks = snap.tl.diff_test(
#     peak_mat,
#     cell_group1=naive_B,
#     cell_group2=memory_B,
#     features=peaks_selected,
# )

# diff_peaks = diff_peaks.filter(pl.col('adjusted p-value') < 0.01)
# diff_peaks.head()

# snap.pl.regions(
#     peak_mat,
#     groupby = 'cell_type',
#     peaks = {
#         group1: diff_peaks.filter(pl.col("log2(fold_change)") > 0)['feature name'].to_numpy(),
#         group2: diff_peaks.filter(pl.col("log2(fold_change)") < 0)['feature name'].to_numpy(),
#     },
#     interactive = False)

# snap.ex.export_coverage(
#     adata,                       # 你的 h5ad
#     groupby="cell_type",         # 也可以是 'leiden' 等
#     bin_size=10,                 # 10-bp 分辨率足够平滑
#     normalization="RPKM",        # 默认；也可 CPM/None
#     out_dir="bw_tracks",         # 输出文件夹
#     suffix=f"{dataset}.bw",                # 每个群体会得到 CellType.bw
#     n_jobs=8
# )


# In[22]:


import pandas as pd
from pathlib import Path

def export_macs3_peaks(macs3_dict: dict,
                       out_dir: str = "peak_tracks",
                       fmt: str = "narrowPeak"):
    """
    把 adata.uns['macs3'] 导出为多个 peak 文件
    ------------------------------------------------
    Parameters
    ----------
    macs3_dict : dict[str, pandas.DataFrame]
        键是 cell_type / cluster 名，值是 DataFrame（MACS3 输出列）
    out_dir : str
        输出文件夹，默认 'peak_tracks'
    fmt : {'narrowPeak', 'bed3'}
        narrowPeak: 10 列完整格式
        bed3      : 3 列 (chrom, start, end)

    Returns
    -------
    List[Path] : 导出的文件路径列表
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    required_cols = ["chrom", "start", "end",
                     "name", "score", "strand",
                     "signal_value", "p_value", "q_value", "peak"]

    written = []
    for label, df in macs3_dict.items():
        # --- 补缺列（少数版本只留前三列） ---
        if fmt == "narrowPeak":
            for col in required_cols:
                if col not in df.columns:
                    # 用占位符填充
                    fill = "." if col in {"name", "strand"} else 0
                    df[col] = fill
            out_file = out_path / f"{label}.narrowPeak"
            df[required_cols].to_csv(out_file, sep="\t",
                                     header=False, index=False)
        elif fmt == "bed3":
            out_file = out_path / f"{label}.bed"
            df.iloc[:, :3].to_csv(out_file, sep="\t",
                                  header=False, index=False)
        else:
            raise ValueError("fmt must be 'narrowPeak' or 'bed3'")
        written.append(out_file)

    print(f"[INFO] Exported {len(written)} peak files → {out_path}")
    return written


# In[23]:


peak_files = export_macs3_peaks(
    adata.uns["macs3"],
    out_dir=f"{output_dir}/peak_tracks",   # 输出目录
    fmt="narrowPeak"        # 或 'bed3'
)


# In[ ]:




