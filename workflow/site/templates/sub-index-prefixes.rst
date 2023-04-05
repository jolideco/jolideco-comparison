{{ title }}
{{ "=" * title|length }}


.. list-table:: 
   :header-rows: 1

   * - Instrument
     - Ground Truth
     {% for method in methods %}
     - {{ method }}
     {% endfor %}

   {% for prefix in prefixes %}

   * - {{ prefix }}
     - .. image:: ../images/flux-true-thumbnail.png
          :align: center
          :width: 200px

     {% for method in methods %}
     - .. image:: {{ "{prefix}/{method}/images/flux-thumbnail.png".format(method=method, prefix=prefix) }}
          :align: center
          :width: 200px
     {% endfor %}

   {% endfor %}


.. toctree::
   :maxdepth: 1

   {% for prefix in prefixes %}
   {% set filename ="{prefix}/index.rst".format(prefix=prefix) %}
   {{ filename|indent(3, False) }}
   {% endfor %}
