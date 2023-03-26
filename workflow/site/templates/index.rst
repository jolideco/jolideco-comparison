Jolideco Comparison
===================

This webpage shows the results from an extensive comparison of the Jolideco method to
other deconvolution algorithms. 


.. toctree::
   :maxdepth: 2

   {% for filename in filenames_toctree %}
   {{ filename|indent(3, False) }}
   {% endfor %}{{ '\n' }}