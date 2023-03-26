Jolideco Comparison
===================

This webpage shows the results from an extensive comparison of the Jolideco method to
other deconvolution algorithms. 


{% for subtitle, filenames in filenames_toctree.items() %}
.. toctree::
   :maxdepth: 1
   :caption: {{ subtitle }}
   {% for filename in filenames %}
   {{ filename|indent(3, False) }}
   {% endfor %}{{ '\n' }}
{% endfor %}
