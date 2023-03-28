Jolideco Comparison
===================

This webpage shows the results from an extensive comparison of the Jolideco method to
other deconvolution algorithms. For the comparison we used:

- Pylira: https://github.com/astrostat/pylira
- Jolideco: https://github.com/jolideco/jolideco

You can navigate through the results with the menu on the left.

{% for subtitle, filenames in filenames_toctree.items() %}
.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: {{ subtitle }}
   {% for filename in filenames %}
   {{ filename|indent(3, False) }}
   {% endfor %}{{ '\n' }}
{% endfor %}
