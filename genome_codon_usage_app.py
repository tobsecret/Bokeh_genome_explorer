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


from Bio import SeqIO
from Bio.SeqUtils import CodonUsage, CodonUsageIndices
from Bio.SeqUtils.MeltingTemp import Tm_GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis


#from sklearn.decomposition import PCA
#from sklearn.preprocessing import StandardScaler

from bokeh.layouts import row
from bokeh.models import LinearColorMapper, CategoricalColorMapper, HoverTool, ColorBar
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
ordered_codon_list = codon_aa_df['codon'].tolist()
color_palette21 = ["#771155", "#114477", "#117777", 
                   "#117744", "#777711", "#774411",   
                   "#771122", "#AA4455","#AA7744", 
                   "#AAAA44",  "#44AA77","#44AAAA", 
                   "#4477AA", "#AA4488","#DD7788", 
                   "#CC99BB", "#77AADD", "#77CCCC", #categorical color palette with 21 colors - one for each codon. 
                   "#88CCAA", "#DDDD77","#DDAA77"]  #Source: https://www.r-bloggers.com/the-paul-tol-21-color-salute/ (I scrambled the palette)

#Importing and pre-processing the data

mtb_HN024 =SeqIO.read('Mycobacterium_tuberculosis_HN_024_genbank_AP018033.1.gb','genbank')
cds = [feature for feature in mtb_HN024.features if feature.type == 'CDS']
cds_extract = [feature.extract(mtb_HN024.seq) for feature in cds]
mtb_cds =  pd.DataFrame(
    {
        'name':[feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else '' for feature in cds],
        'sequence':[str(feature) for feature in cds_extract],
        'translation':[str(feature.translate(to_stop=True, table=11)) for feature in cds_extract],
        'ref':mtb_HN024.annotations['references'][0].title,
        'start':[feature.location.nofuzzy_start for feature in cds],
        'end':[feature.location.nofuzzy_end for feature in cds],
        'strand':[feature.strand for feature in cds],
        'product':[feature.qualifiers['product'][0] for feature in cds],
        'locus_tag':[feature.qualifiers['locus_tag'][0] for feature in cds],
    }, 
    columns=['name', 'product', 'locus_tag','sequence','translation', 'start', 'end', 'strand','ref']
    )

mtb_cds['gene_length'] = mtb_cds['sequence'].apply(len)
mtb_cds['aromaticity'] = mtb_cds['translation'].apply(lambda x: ProteinAnalysis(str(x)).aromaticity())
mtb_cds['isoelectric_point'] = mtb_cds['translation'].apply(lambda x: ProteinAnalysis(str(x)).isoelectric_point())
mtb_cds['instability'] = mtb_cds['translation'].apply(lambda x: ProteinAnalysis(str(x)).instability_index())
mtb_cds['middle_point'] = (mtb_cds['start']+(mtb_cds['end']- mtb_cds['start'])/2).astype(int)
mtb_cds = mtb_cds.loc[mtb_cds['gene_length']%3==0]
#print(mtb_cds.head())
codon_mtb_df = pd.DataFrame(mtb_cds['sequence'].apply(
  lambda sequence: pd.Series(
    simple_string_wrap(sequence, width=3)).value_counts() 
  )
).reindex_axis(ordered_codon_list, axis=1)
#print(codon_mtb_df.head())
codon_mtb_df = codon_mtb_df.fillna(0).astype(int)

mtb_cds = mtb_cds.join(codon_mtb_df, how='outer')



# Create the scatter plot
gene_scatter = figure(tools=TOOLS, plot_width=550, plot_height=600, min_border=10, min_border_left=50,
           toolbar_location="above", #x_axis_location=None, y_axis_location=None,
           title="Properties of coding sequences", x_axis_label='instability', y_axis_label='isoelectric_point')
gene_scatter.background_fill_color = "#fafafa"

hover = HoverTool(tooltips=[
                  ('gene name','@name'),
                  ('length', '@gene_length')])
scatter_points = gene_scatter.scatter(x='instability', y='isoelectric_point', size=5, color="#3A5785", alpha=0.6, source=mtb_cds)
gene_scatter.add_tools(hover)

# Create the bar chart
datapoint = mtb_cds.loc[:,ordered_codon_list].iloc[0,:].astype(int).reset_index().rename(columns={'index':'codon',0:'counts'})
datapoint['amino acid'] = datapoint.apply(lambda x: codon_aa_dict[x['codon']], axis=1)
datapoint.sort_values(by='amino acid', inplace=True)
datapoint['counts'] = 0

codon_dist_cd = ColumnDataSource(datapoint)
codon_dist = figure(y_range=datapoint.sort_values(by='amino acid')['codon'].tolist(),
          toolbar_location=None, tools='', plot_height=gene_scatter.plot_height, plot_width=400)

colormapper = CategoricalColorMapper(factors=datapoint['amino acid'].unique().tolist(), 
                                     palette=color_palette21)
codon_hbar = codon_dist.hbar(
                            y='codon', height=0.5, 
                            left=0, right='counts', 
                            color={'field':'amino acid', 'transform':colormapper}, 
                            source=codon_dist_cd,
                            legend='amino acid'
                            )
codon_dist.legend.location = (0,-18)
codon_dist.legend.border_line_color = None
codon_dist.right.append(codon_dist.legend[0])

#Create the donut plot
seq_partition_length = round(len(mtb_HN024)/1500)
mtb_genome_nmers = pd.Series(simple_string_wrap(str(mtb_HN024.seq), width=seq_partition_length))
nmers =  mtb_genome_nmers.apply(Tm_GC).reset_index().rename(columns={'index':'nmer_number',0:'melting_temperature'})
angle_increment = 2*pi/len(nmers['nmer_number'])

nmers['end_angle'] = 0.5*pi-nmers['nmer_number']*angle_increment
nmers['start_angle'] = nmers['end_angle']-angle_increment
nmers['start_base'] = nmers['nmer_number']*seq_partition_length+1
nmers['end_base'] = nmers['start_base']+seq_partition_length-1
nmers.loc[nmers.index[-1],'end_base'] = len(mtb_HN024) #The last end_base value should be identical to
                                                                          #the length of the genome.

nmers_cd = ColumnDataSource(nmers)
tm_min = nmers['melting_temperature'].min()
tm_max = nmers['melting_temperature'].max()

melting_temp_color_mapper = LinearColorMapper(low=tm_min, high=tm_max,palette=Magma10)
#Making some colormappers for the mtb_cds dataframe columns
linearcolor_columns = ['aromaticity', 'isoelectric_point', 'instability']
linearcolormapperdict = dict([[column, 
                              LinearColorMapper(low=mtb_cds[column].min(), high=mtb_cds[column].max(), palette="Oranges4")]
                             for column in linearcolor_columns])

#Configuring the figure axes
donut_plot = figure(tools='zoom_in,zoom_out,box_zoom,reset,pan' ,plot_width=500, plot_height=gene_scatter.plot_height,
                    title='Annotated Genome', 
                    x_axis_location=None, y_axis_location=None, 
                    background_fill_color='black',
                    toolbar_location='above'
                   #background_fill_alpha=0.75
                  )
donut_plot.xgrid.grid_line_color = None
donut_plot.ygrid.grid_line_color = None

#Making the donut plot with the melting temperature and the annotations.
melt_temp_renderer = donut_plot.annular_wedge(x=0, y=0, inner_radius=0.25, outer_radius=0.5,
                                             start_angle='start_angle',
                                             end_angle='end_angle', 
                                             direction='anticlock',
                                             color={'field':'melting_temperature',
                                                    'transform':melting_temp_color_mapper},
                                             alpha=0.6, 
                                             source= nmers_cd,
                                             legend='melting temperature'
                                             )

#configuring the HoverTool
melt_temp_hover = HoverTool(tooltips=[('start base','@start_base'), 
                                 ('end base', '@end_base'),
                                 ('melting temp', '@melting_temperature')
                                 
                                ],
                           renderers=[melt_temp_renderer])
donut_plot.add_tools(melt_temp_hover)

#Configuring the ColorBar
melt_temp_cbar = ColorBar(color_mapper=melting_temp_color_mapper, label_standoff=12, border_line_color=None, location=(0,0))
donut_plot.add_layout(melt_temp_cbar, 'right')

angle_increment = 2*pi/len(mtb_HN024)
mtb_cds['end_angle'] = 0.5*pi - mtb_cds['start']*angle_increment
mtb_cds['start_angle'] = 0.5*pi - mtb_cds['end']*angle_increment
mtb_cds['inner_radius'] = 0.55 + 0.1 + mtb_cds['strand']*0.11
mtb_cds['outer_radius'] = mtb_cds['inner_radius'] +0.2
mtb_anno_cd = ColumnDataSource(mtb_cds)

anno_renderer = donut_plot.annular_wedge(x=0, 
                                        y=0,
                                        inner_radius='inner_radius',
                                        outer_radius='outer_radius',
                                        start_angle='start_angle',
                                        end_angle='end_angle', 
                                        direction='anticlock',
                                        source= mtb_anno_cd,
                                        alpha=(0.3),
                                        color={'field':'aromaticity', 
                                               'transform':linearcolormapperdict['aromaticity']},
                                        legend='annotations'
                                       )
anno_hover = HoverTool(tooltips=[('name','@name'),
                                 ('locus tag', '@locus_tag'),
                                 ('product', '@product'),
                                 ('start', '@start'),
                                 ('end', '@end')
                                 
                                ], renderers=[anno_renderer])
donut_plot.add_tools(anno_hover)
donut_plot.legend.click_policy = 'hide'


#Create the layout of the plots in the browser
layout = row(donut_plot ,gene_scatter, codon_dist)

#adding layout to the server environment
curdoc().add_root(layout)
curdoc().title = "Genome Explorer"

#Updating the data
def update(attr, old, new):
    inds = list(new['1d']['indices']) #the index of the new row of data in the data source.
    if len(inds) == 0:
        new_counts = np.zeros(64)
    else:
        new_counts = mtb_cds.loc[:,ordered_codon_list].iloc[inds,:].sum().reset_index().rename(columns={'index':'codon',0:'counts'})['counts'].tolist()

    codon_hbar.data_source.data['counts'] = new_counts
    

scatter_points.data_source.on_change('selected', update)