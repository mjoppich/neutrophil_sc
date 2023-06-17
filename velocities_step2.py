import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
from natsort import natsorted
import itertools
import sys
import argparse

scv.set_figure_params()

LW = 0.3

import math
import matplotlib.pyplot as plt
from matplotlib.colors import ColorConverter, to_hex


def gamma_correct(u):
    # Standard CRT Gamma
    GAMMA = 2.4
    if u > 0.00304:
        u = (1.055*u ** (1/GAMMA)) - 0.055
    else:
        u = 12.92*u
    return u


def hcl2rgb(h,c,l):
    # ADAPTED FOR PYTHON BY MARKUS JOPPICH
    # 
    # HCL2RGB Convert a HCL (i.e., CIELUV) color space value to one
    #   in sRGB space.
    #   RGB = HCL2RGB(H, C, L) will convert the color (H, C, L) in
    #   HCL color space to RGB = [R, G, B] in sRGB color space.
    #   Values that lie outside sRGB space will be silently corrected.
    # Code written by Nicholas J. Hughes, 2014, released under the following
    # licence.
    #
    # The MIT License (MIT)
    #
    # Copyright (c) 2014 Nicholas J. Hughes
    # 
    # Permission is hereby granted, free of charge, to any person obtaining a copy
    # of this software and associated documentation files (the "Software"), to deal
    # in the Software without restriction, including without limitation the rights
    # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    # copies of the Software, and to permit persons to whom the Software is
    # furnished to do so, subject to the following conditions:
    # 
    # The above copyright notice and this permission notice shall be included in
    # all copies or substantial portions of the Software.
    # 
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    # THE SOFTWARE.
    # D65 White Point
    WHITE_Y = 100.000
    WHITE_u = 0.1978398
    WHITE_v = 0.4683363
    if l < 0 or l > WHITE_Y or c < 0:
        print("Invalid CIE-HCL color.")
        assert(False)
    L = float(l)
    U = c * math.cos(math.radians(h))
    V = c * math.sin(math.radians(h))
    if L <= 0 and U == 0 and V == 0:
        X = 0
        Y = 0
        Z = 0
    else:
        Y = WHITE_Y
        if L > 7.999592:
            Y = Y*(((L + 16)/116) ** 3.0)
        else:
            Y = Y*L/903.3
        
        u = U/(13*L) + WHITE_u
        v = V/(13*L) + WHITE_v
        X = (9.0*Y*u)/(4*v)
        Z = -X/3 - 5*Y + 3*Y/v
    # Now convert to sRGB
    r = gamma_correct((3.240479*X - 1.537150*Y - 0.498535*Z)/WHITE_Y)
    g = gamma_correct((-0.969256*X + 1.875992*Y + 0.041556*Z)/WHITE_Y)
    b = gamma_correct((0.055648*X - 0.204043*Y + 1.057311*Z)/WHITE_Y)
    # Round to integers and correct
    r = max([min([round(255 * r), 255]), 0])
    g = max([min([round(255 * g), 255]), 0])
    b = max([min([round(255 * b), 255]), 0])   
    rgb = [x/255.0 for x in [r, g, b]]
    #rgb = [r,g,b]
    print(rgb)
    return rgb


def hue_pal(n=1, h = [15, 375], c = 100, l = 65):
    print(h)
    assert(len(h) == 2)
    if ((h[1]-h[0] % 360) < 1):
        h[1] = h[1]-360/n
    print(h)
    hues = []
    curH = h[0]
    while curH < h[1]:
        hues.append((curH, c, l))
        curH += (h[1]-h[0])/n      
    hexColors = []
    for x in hues:
        rgbColor = hcl2rgb( x[0], x[1], x[2] )
        hexColor = to_hex(rgbColor)
        hexColors.append(hexColor)
    print(hexColors)
    return hexColors


def getClusterColors(numClusters):
    useColors = hue_pal(numClusters)
    cluster2color = {}
    for x in range(0, numClusters):
        cluster2color[x] = useColors[x]
    return cluster2color



parser = argparse.ArgumentParser(description='Automagic Velocity Analysis.')
parser.add_argument("-i", "--inprefix", type=str, required=True)
parser.add_argument("-o", "--outprefix", type=str, required=True)
parser.add_argument("-c", "--clusters", type=str, nargs="+", default=None, required=False)
args = parser.parse_args()

#
## Setup Object for velocity analysis
#
adata = scv.read_loom("{}.loom".format(args.inprefix))

u = scv.read("{}_unspliced.csv".format(args.inprefix))
s = scv.read("{}_RNA.csv".format(args.inprefix))
X_umap = scv.load('{}_umap_embeddings.csv'.format(args.inprefix), index_col=0)

cells = list(adata.obs.index)
genes = list(adata.var.index)

usubset = u[u.obs.index.isin(cells)]
ssubset = s[s.obs.index.isin(cells)]

adata.obsm['X_umap'] = X_umap.loc[adata.obs_names].values
adata.layers['unspliced'] = usubset.X.T[u.var.index.isin(genes)].T
adata.layers['spliced'] = adata.X #ssubset.X.T[s.var.index.isin(genes)].T

if not min(set(adata.obs.identsstr.astype(str).astype(int))) == 0:
    adata.obs["idents"] = [str(x-1) for x in adata.obs.identsstr.astype(str).astype(int)]
else:
    adata.obs["idents"] = [str(x) for x in adata.obs.identsstr]


numClusters = len(set(adata.obs.identsstr))
print("Identified", numClusters, "clusters:", set(adata.obs.identsstr))
ident2col = getClusterColors(numClusters)
print(ident2col)
sortedClusters = natsorted(set([str(x) for x in adata.obs.identsstr]))
print(sortedClusters)
adata.uns['idents_colors']={sortedClusters[ix]: ident2col[x] for ix, x in enumerate(ident2col)}

adata.write("{}_adata.h5ad".format(args.outprefix))
# adata_full = sc.read("patint_thr_adata_final.h5ad")
# adata = adata_full[adata_full.obs.identsstr.isin(["0","1","2","3","9"]), :]
#
## Prepare object for velocity analysis
#

if not args.clusters is None:
    print(args.clusters)
    adata = adata[adata.obs.identsstr.isin(args.clusters), :].copy()


sc.pp.highly_variable_genes(adata)

adata_orig = adata
adata = adata_orig[:, adata_orig.var_names[adata_orig.var['highly_variable']]].copy()


#
## Velocity Analysis
#

class ArgEmulator:
    def __init__(self):
        self.outprefix = "patint_reduced01239"


args = ArgEmulator()

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)


scv.pl.velocity_embedding_stream(adata, basis="umap", color="idents", dpi=300, save="{}_scvelo_idents.png".format(args.outprefix))
scv.pl.velocity_embedding_stream(adata, basis="umap", color="idents", dpi=300, save="{}_scvelo_idents.svg".format(args.outprefix))


#
## Velocity Confidence
#

scv.pl.proportions(adata, groupby='idents', dpi=300, save="{}_scvelo_proportions.png".format(args.outprefix))
scv.pl.proportions(adata, groupby='idents', dpi=300, save="{}_scvelo_proportions.svg".format(args.outprefix))

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
df = adata.obs.groupby('idents')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], dpi=300, save="{}_scvelo_velocity_confidence.png".format(args.outprefix))
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], dpi=300, save="{}_scvelo_velocity_confidence.svg".format(args.outprefix))


#
## Rank Velocity Genes
#
scv.tl.rank_velocity_genes(adata, groupby="idents" )
rvgDF = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
rvgGenes = set(list(itertools.chain(*rvgDF.values)))
print("Ranked velocity genes", rvgGenes)
if len(rvgGenes) > 0:
    scv.pl.velocity(adata, rvgGenes,  ncols=2, color="idents", add_outline=True, save="{}_scvelo_velocity_rvggenes.png".format(args.outprefix))
    scv.pl.velocity(adata, rvgGenes,  ncols=2, color="idents", add_outline=True, save="{}_scvelo_velocity_rvggenes.svg".format(args.outprefix))


#
## Gene-wise velocity
#
scv.pl.velocity(adata, ['CXCR4', 'TLR4', 'PTPRC', 'MME', 'CXCR2', 'ITGAM', 'ICAM1', 'FCGR3A', 'FCGR3B'], color="idents", ncols=1, add_outline=True, save="{}_scvelo_velocity_genes.png".format(args.outprefix))
scv.pl.velocity(adata, ['CXCR4', 'TLR4', 'PTPRC', 'MME', 'CXCR2', 'ITGAM', 'ICAM1', 'FCGR3A', 'FCGR3B'], color="idents", ncols=1, add_outline=True, save="{}_scvelo_velocity_genes.svg".format(args.outprefix))


x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=1000,n_steps=50,random_state=1)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)
plt.savefig("figures/{}_scvelo_cell_transitions.png".format(args.outprefix))
plt.savefig("figures/{}_scvelo_cell_transitions.svg".format(args.outprefix))

#
## Pseudotime
#

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save="{}_scvelo_pseudotimes.png".format(args.outprefix))
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save="{}_scvelo_pseudotimes.svg".format(args.outprefix))

adata.obs["groups"] = adata.obs.idents

#
## PAGA
#

# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='idents')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, dpi=300, save="{}_scvelo_paga.png".format(args.outprefix))
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, dpi=300, save="{}_scvelo_paga.svg".format(args.outprefix))


#
## Recovering dynamics + PAGA
#
goi = ['LSEL', 'CXCR4', 'TLR4', 'PTPRC', 'MME', 'CXCR2', 'ITGAM', 'ICAM1', 'FCGR3A', 'FCGR3B']
scv.tl.recover_dynamics(adata)#, var_names=goi)
scv.pl.scatter(adata, basis=goi, vkey='dynamics', color="idents", linewidth=5, frameon=False, size=120, legend_loc='none', fontsize=16, save="{}_scvelo_recoverdynamics.png".format(args.outprefix))
scv.pl.scatter(adata, basis=goi, vkey='dynamics', color="idents", linewidth=5, frameon=False, size=120, legend_loc='none', fontsize=16, save="{}_scvelo_recoverdynamics.svg".format(args.outprefix))

scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="{}_scvelo_rd_latenttime.png".format(args.outprefix))
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="{}_scvelo_rd_latenttime.svg".format(args.outprefix))

scv.tl.velocity(adata, diff_kinetics=True)
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding(adata, dpi=120, arrow_size=4, arrow_length=2, color="idents", save="{}_scvelo_rd_velocity.png".format(args.outprefix))
scv.pl.velocity_embedding(adata, dpi=120, arrow_size=4, arrow_length=2, color="idents", save="{}_scvelo_rd_velocity.svg".format(args.outprefix))

scv.pl.velocity_embedding_stream(adata, basis="umap", color="idents", dpi=300, save="{}_scvelo_rd_velocitydynamics.png".format(args.outprefix))
scv.pl.velocity_embedding_stream(adata, basis="umap", color="idents", dpi=300, save="{}_scvelo_rd_velocitydynamics.svg".format(args.outprefix))


# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='idents')

df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save="{}_scvelo_rd_paga.png".format(args.outprefix))
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save="{}_scvelo_rd_paga.svg".format(args.outprefix))


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='idents', n_convolve=100, save="{}_scvelo_rd_heatmap.png".format(args.outprefix), figsize = (8, 16))
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='idents', n_convolve=100, save="{}_scvelo_rd_heatmap.svg".format(args.outprefix), figsize = (8, 16))

if len(rvgGenes) > 0:
    scv.pl.velocity(adata, rvgGenes,  ncols=2, color="idents", add_outline=True, save="{}_scvelo_rd_rvggenes2.png".format(args.outprefix))
    scv.pl.velocity(adata, rvgGenes,  ncols=2, color="idents", add_outline=True, save="{}_scvelo_rd_rvggenes2.svg".format(args.outprefix))

adata.write("{}_adata_final.h5ad".format(args.outprefix))



#adata = scv.read("patint_thr_adata_final.h5ad")
#scv.pl.velocity(adata, ["CD177", "CXCR4"],  ncols=2, color="idents", add_outline=True, save="targeted_velocities2.png")
#scv.pl.velocity(adata, ["CD177", "CXCR4"],  ncols=2, color="idents", add_outline=True, save="targeted_velocities.svg")