{{ title }}
{{ "=" * title|length }}

.. toctree::
   :maxdepth: 2

   {% for filename in filenames_toctree %}
   {{ filename|indent(3, False) }}
   {% endfor %}{{ '\n' }}