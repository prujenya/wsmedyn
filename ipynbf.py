__author__ = 'jenya'

import graphlab

sales = graphlab.SFrame('home_data.gl/')

sales

#graphlab.canvas.set_target('ipynb')#
sales.show(view='Scatter Plot',x='sqft_living',y='price')
