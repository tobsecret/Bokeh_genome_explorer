'''
This is a simple in-browser app that:
1. loads up a fasta file with coding sequences, 
2. calculates the z-score-normalized codon distribution in all coding sequences
3. applies PCA to the dataframe with the codon distributions
4. plots a scatter plot of the PCA data
5. serves the codon distribution of a single selected data point in form of a bar diagram
'''

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.SeqUtils import CodonUsage, CodonUsageIndices
import textwrap
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from bokeh.layouts import row
from bokeh.models import Spacer, CategoricalColorMapper
from bokeh.models.sources import ColumnDataSource
from bokeh.plotting import figure, curdoc



# create three normal population samples with different parameters

#General settings:
TOOLS="pan,tap,zoom_in,zoom_out,reset"
codons_dict = CodonUsage.CodonsDict.copy()
codon_list = list(codons_dict.keys())
synonymous_codons = CodonUsage.SynonymousCodons.copy()
codon_aa_dict = dict( (v,k) for k in synonymous_codons for v in synonymous_codons[k] )
color_palette21 = ["#771155", "#114477", "#117777", 
                   "#117744", "#777711", "#774411",   
                   "#771122", "#AA4455","#AA7744", 
                   "#AAAA44",  "#44AA77","#44AAAA", 
                   "#4477AA", "#AA4488","#DD7788", 
                   "#CC99BB", "#77AADD", "#77CCCC", 
                   "#88CCAA", "#DDDD77","#DDAA77"] #categorical color palette with 21 colors - one for each codon.

#Importing and pre-processing the data
smeg_cds = list(SeqIO.parse('Mycobacterium_smegmatis_coding_sequences_NC_008596.fasta.txt', 'fasta'))
df = pd.DataFrame([pd.Series(textwrap.wrap(str(record.seq), width=3)).value_counts() 
                for record in smeg_cds if (len(record.seq)%3==0)]
                ).fillna(0) #Splitting all the records into their individual codons and counting how often each codon comes up in each record. NaNs are replaced with 0.
df = df.loc[pd.Series([len(record.seq)%3 for record in smeg_cds])==0].astype(int) #dropping all coding sequences that have a number of bases that is not divisible by 3.
#df.loc['amino acids'] = pd.Series(df.columns.tolist()).apply(lambda x: codon_aa_dict[x['codon']])
array_stand = StandardScaler().fit_transform(df) #Doing a Z-score normalization of our data.
df_pca = pd.DataFrame(PCA(n_components=10).fit_transform(array_stand))
df_pca.columns = [str(col) for col in df_pca.columns.tolist()]
joined = df_pca.join(df, how='outer')


# create the scatter plot
pca_scatter = figure(tools=TOOLS, plot_width=600, plot_height=600, min_border=10, min_border_left=50,
           toolbar_location="above", x_axis_location=None, y_axis_location=None,
           title="PCA of codon distributions")
pca_scatter.background_fill_color = "#fafafa"

scatter_points = pca_scatter.scatter(joined['0'], joined['1'], size=5, color="#3A5785", alpha=0.6)

'''
# create the horizontal histogram
hhist, hedges = np.histogram(x, bins=20)
hzeros = np.zeros(len(hedges)-1)
hmax = max(hhist)*1.1

LINE_ARGS = dict(color="#3A5785", line_color=None)

ph = figure(toolbar_location=None, plot_width=p.plot_width, plot_height=200, x_range=p.x_range,
            y_range=(-hmax, hmax), min_border=10, min_border_left=50, y_axis_location="right")
ph.xgrid.grid_line_color = None
ph.yaxis.major_label_orientation = np.pi/4
ph.background_fill_color = "#fafafa"

ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hhist, color="white", line_color="#3A5785")
hh1 = ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hzeros, alpha=0.5, **LINE_ARGS)
hh2 = ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hzeros, alpha=0.1, **LINE_ARGS)
'''

# create the vertical histogram
datapoint = joined.loc[:,codon_list].iloc[0,:].astype(int).reset_index().rename(columns={'index':'codon',0:'counts'})
datapoint['amino acid'] = datapoint.apply(lambda x: codon_aa_dict[x['codon']], axis=1)
#datapoint['counts'] = 0
print(datapoint.shape)
codon_dist_cd = ColumnDataSource(datapoint.sort_values(by='amino acid'))
codon_dist = figure(y_range=datapoint.sort_values(by='amino acid')['codon'].tolist(),
          toolbar_location=None, tools='', plot_height=pca_scatter.plot_height, plot_width=400)

colormapper = CategoricalColorMapper(factors=pd.Series(codon_dist_cd.data['amino acid'].tolist()).unique().tolist(), 
                                     palette=color_palette21)
codon_hbar = codon_dist.hbar(y='codon', height=0.5, left=0, right='counts', 
          color={'field':'amino acid', 'transform':colormapper}, 
          source=codon_dist_cd,
         legend='amino acid'
         )
codon_dist.legend.location = (0,70)
codon_dist.legend.border_line_color = None
codon_dist.right.append(codon_dist.legend[0])

layout = row(pca_scatter, codon_dist)

curdoc().add_root(layout)
curdoc().title = "Selection Histogram"
print('joined.shape',joined.shape)
#Updating the data
def update(attr, old, new):
    inds = list(new['1d']['indices']) #the index of the new row of data in the data source.
    print(inds)
    if len(inds) == 0:
        new_counts = np.zeros(64)
    elif len(inds)==1:
        print('inds[0] in joined.index',inds[0] in joined.index)
        temp = joined.loc[:,codon_list].iloc[inds,:].sum().astype(int).reset_index().rename(columns={'index':'codon',0:'counts'})
        temp['amino acid'] =  temp.apply(lambda x: codon_aa_dict[x['codon']], axis=1)
        temp.sort_values(by='amino acid')
        print(temp.head())
        new_counts = temp['counts'].values
    elif len(inds) > 1:
        temp = joined.loc[:,codon_list].iloc[inds,:].sum().astype(int).reset_index().rename(columns={'index':'codon',0:'counts'})
        temp['amino acid'] =  temp.apply(lambda x: codon_aa_dict[x['codon']], axis=1)
        temp.sort_values(by='amino acid')
        new_counts = temp['counts'].values
    codon_hbar.data_source.data['counts']   =  new_counts
    

scatter_points.data_source.on_change('selected', update)