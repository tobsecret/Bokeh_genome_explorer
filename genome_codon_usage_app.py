'''
This is a simple in-browser app that:
1. loads up a fasta file with coding sequences, 
2. calculates the z-score-normalized codon distribution in all coding sequences
3. applies PCA to the dataframe with the codon distributions
4. plots a scatter plot of the PCA data
5. serves the codon distribution of a single selected data point in form of a bar diagram
'''
#Imports from standard library
from math import pi

#Imports from other libraries
import numpy as np
import pandas as pd


from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO #this line should be removed once I have started using the SimpleFastaParser.
from Bio.SeqUtils import CodonUsage, CodonUsageIndices
from Bio.SeqUtils.MeltingTemp import Tm_GC


from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from bokeh.layouts import row
from bokeh.models import LinearColorMapper, CategoricalColorMapper, HoverTool
from bokeh.models.sources import ColumnDataSource
from bokeh.plotting import figure, curdoc
from bokeh.palettes import Magma10


def simple_string_wrap(string, width=50):
    '''
    Takes a single string and returns a generator of consecutive substrings
    with a length of width.
    '''
    return (string[i:i+width] for i in range(0, len(string), width))

#General settings:
TOOLS="pan,tap,box_select,zoom_in,zoom_out,reset"
codons_dict = CodonUsage.CodonsDict.copy()
codon_list = list(codons_dict.keys())
synonymous_codons = CodonUsage.SynonymousCodons.copy()
codon_aa_dict = dict( (v,k) for k in synonymous_codons for v in synonymous_codons[k] )
codon_aa_df = pd.DataFrame(data={'codon'    : list(codon_aa_dict.keys()), 
                  'amino acid': list(codon_aa_dict.values())}).sort_values(by='amino acid')
color_palette21 = ["#771155", "#114477", "#117777", 
                   "#117744", "#777711", "#774411",   
                   "#771122", "#AA4455","#AA7744", 
                   "#AAAA44",  "#44AA77","#44AAAA", 
                   "#4477AA", "#AA4488","#DD7788", 
                   "#CC99BB", "#77AADD", "#77CCCC", #categorical color palette with 21 colors - one for each codon. 
                   "#88CCAA", "#DDDD77","#DDAA77"]  #Source: https://www.r-bloggers.com/the-paul-tol-21-color-salute/ (I scrambled the palette)

#Importing and pre-processing the data
smeg_cds = list(SeqIO.parse('Mycobacterium_smegmatis_coding_sequences_NC_008596.fasta.txt', 'fasta'))
df = pd.DataFrame([pd.Series(simple_string_wrap(str(record.seq), width=3)).value_counts() 
                for record in smeg_cds if (len(record.seq)%3==0)]
                ).fillna(0) #Splitting all the records into their individual codons and counting how often each codon comes up in each record. NaNs are replaced with 0.
df = df.loc[pd.Series([len(record.seq)%3 for record in smeg_cds])==0].astype(int) #dropping all coding sequences that have a number of bases that is not divisible by 3.
df = df.reindex_axis(codon_aa_df['codon'].tolist(),axis=1)


array_stand = StandardScaler().fit_transform(df) #Doing a Z-score normalization of our data.
df_pca = pd.DataFrame(PCA(n_components=10).fit_transform(array_stand))
df_pca.columns = [str(col) for col in df_pca.columns.tolist()]

joined = df_pca.join(df, how='outer')
joined['gene'] = pd.Series([record.description.split('gene=')[1].split(']')[0] for record in smeg_cds])
joined['gene_length'] = pd.Series([len(record.seq) for record in smeg_cds])

smeg_genome = list(SeqIO.parse('Mycobacterium_smegmatis_whole_genome_NC_008596.fasta', 'fasta'))


# Create the scatter plot
pca_scatter = figure(tools=TOOLS, plot_width=600, plot_height=600, min_border=10, min_border_left=50,
           toolbar_location="above", x_axis_location=None, y_axis_location=None,
           title="PCA of codon distributions")
pca_scatter.background_fill_color = "#fafafa"

hover = HoverTool(tooltips=[
                  ('gene name','@gene'),
                  ('length', '@gene_length')])
scatter_points = pca_scatter.scatter(x='0', y='1', size=5, color="#3A5785", alpha=0.6, source=joined)
pca_scatter.add_tools(hover)

# Create the bar chart
datapoint = joined.loc[:,codon_list].iloc[0,:].astype(int).reset_index().rename(columns={'index':'codon',0:'counts'})
datapoint['amino acid'] = datapoint.apply(lambda x: codon_aa_dict[x['codon']], axis=1)
datapoint.sort_values(by='amino acid', inplace=True)
datapoint['counts'] = 0

codon_dist_cd = ColumnDataSource(datapoint)
codon_dist = figure(y_range=datapoint.sort_values(by='amino acid')['codon'].tolist(),
          toolbar_location=None, tools='', plot_height=pca_scatter.plot_height, plot_width=400)

colormapper = CategoricalColorMapper(factors=datapoint['amino acid'].unique().tolist(), 
                                     palette=color_palette21)
codon_hbar = codon_dist.hbar(
                            y='codon', height=0.5, 
                            left=0, right='counts', 
                            color={'field':'amino acid', 'transform':colormapper}, 
                            source=codon_dist_cd,
                            legend='amino acid'
                            )
codon_dist.legend.location = (0,70)
codon_dist.legend.border_line_color = None
codon_dist.right.append(codon_dist.legend[0])

#Create the ring plot
seq_partition_length = 50000
genome50mers = pd.Series(simple_string_wrap(str(smeg_genome[0].seq), width=seq_partition_length))
tms_50mers =  genome50mers.apply(Tm_GC).reset_index().rename(columns={'index':'50mer_number',0:'melting temperature'})

angle_increment = 2*pi/len(tms_50mers['50mer_number'])

tms_50mers['end_angle'] = 0.5*pi-tms_50mers['50mer_number']*angle_increment
tms_50mers['start_angle'] = tms_50mers['end_angle']-angle_increment
tms_50mers['start_base'] = tms_50mers['50mer_number']*seq_partition_length+1
tms_50mers['end_base'] = tms_50mers['start_base']+seq_partition_length-1
tms_50mers.loc[tms_50mers.index[-1],'end_base'] = len(smeg_genome[0].seq) #The last end_base value should be identical to
                                                                          #the length of the genome.

tms_50mers_cd = ColumnDataSource(tms_50mers)
tm_min = tms_50mers['melting temperature'].min()
tm_max = tms_50mers['melting temperature'].max()

melting_temp_color_mapper = LinearColorMapper(low=tm_min, high=tm_max,palette=Magma10)

ring_plot = figure(
                  tools='zoom_in,zoom_out,reset,pan', 
                  plot_height=pca_scatter.plot_height, plot_width=400,
                  title='Melting temperature of genome', 
                  x_axis_location=None, y_axis_location=None
                  )
ring_plot.xgrid.grid_line_color = None
ring_plot.ygrid.grid_line_color = None

ring_plot.annular_wedge(x=1, y=1, inner_radius=0.25, outer_radius=0.5,
                        start_angle='start_angle', 
                        end_angle='end_angle', 
                        direction='anticlock',
                        color={'field':'melting temperature', 'transform':melting_temp_color_mapper}, 
                        alpha=0.6, source=tms_50mers_cd)

ring_hover = HoverTool(tooltips=[('start base','@start_base'), 
                                 ('end base', '@end_base')
                                ])
ring_plot.add_tools(ring_hover)

#Create the layout of the plots in the browser
layout = row(ring_plot ,pca_scatter, codon_dist)

#adding layout to the server environment
curdoc().add_root(layout)
curdoc().title = "Genome Explorer"

#Updating the data
def update(attr, old, new):
    inds = list(new['1d']['indices']) #the index of the new row of data in the data source.
    if len(inds) == 0:
        new_counts = np.zeros(64)
    else:
        new_counts = joined.loc[:,codon_list].iloc[inds,:].sum().reset_index().rename(columns={'index':'codon',0:'counts'})['counts'].tolist()

    codon_hbar.data_source.data['counts'] = new_counts
    

scatter_points.data_source.on_change('selected', update)