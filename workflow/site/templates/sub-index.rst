{{ title }}
{{ "=" * title|length }}

.. toctree::
   :maxdepth: 1

   {% for filename in filenames_toctree %}
   {{ filename|indent(3, False) }}
   {% endfor %}{{ '\n' }}