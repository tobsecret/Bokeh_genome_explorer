# Bokeh_genome_explorer
This is my first Bokeh server app.

Here is what it does:
1. loads up a fasta file with coding sequences (DNA) 
2. calculates the z-score-normalized codon distribution in all coding sequences
3. applies PCA to the dataframe with the codon distributions
4. plots a scatter plot of the PCA data
5. serves the codon distribution of a single selected data point in form of a bar diagram

If you make a copy of the [environment](https://github.com/tobsecret/Bokeh_genome_explorer/blob/master/package_requirements_for_humans.txt) 
in which I worked on this app, you can run `bokeh serve genome_codon_usage_app.py --show` in the directory to which you cloned the app.

That should give you something like the screenshot below in which I used the zoom tool to zoom in a little bit and marked a single point.
![alt text](https://github.com/tobsecret/Bokeh_genome_explorer/blob/master/screenshot2.png "Browser Screenshot")
