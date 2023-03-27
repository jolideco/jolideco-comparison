{{ title }}
{{ "=" * title|length }}

Dataset
-------

TODO: Add dataset image

Results Overview
----------------
.. list-table:: 
   :header-rows: 1

   * - Ground Truth
     {% for method in methods %}
     - {{ method }}
     {% endfor %}

   * - .. image:: ../../flux-thumbnail.png
          :align: center

     {% for method in methods %}
     - .. image:: {{ "{method}/images/flux-thumbnail.png".format(method=method) }}
          :align: center
          :target: {{ "{method}/index.html".format(method=method) }}
     {% endfor %}

   * - 
     {% for method in methods %}
     - Value = 1234
     {% endfor %}


.. toctree::
   :hidden:
   :maxdepth: 1

   {% for method in methods %}
   {% set filename ="{method}/index.rst".format(method=method) %}
   {{ filename|indent(3, False) }}
   {% endfor %}