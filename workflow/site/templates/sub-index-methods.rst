{{ title }}
{{ "=" * title|length }}


Results Overview
----------------
.. list-table:: 
   :header-rows: 1

   * - Ground Truth
     {% for method in methods %}
     - {{ method }}
     {% endfor %}

   * - .. image:: ../../images/flux-true-thumbnail.png
          :align: center
          :width: 200px

     {% for method in methods %}
     - .. image:: {{ "{method}/images/flux-thumbnail.png".format(method=method) }}
          :align: center
          :width: 200px
     {% endfor %}

   * -
     {% for method in methods %}
     - :doc:`Details <{{ "{method}/index.rst".format(method=method) }}>`  
     {% endfor %}
  
   * - 
     {% for method in methods %}
     - Value = 1234
     {% endfor %}


Dataset
-------

.. list-table:: 
   :header-rows: 1

   * - Counts
     - PSF
     - Exposure
     - Background

   * - .. figure:: images/counts.png
          :width: 200px
          :align: center
     
     - .. figure:: images/psf.png
          :width: 200px
          :align: center
     
     - .. figure:: images/exposure.png
          :width: 200px
          :align: center
     
     - .. figure:: images/background.png
          :width: 200px
          :align: center

.. toctree::
   :hidden:
   :maxdepth: 1

   {% for method in methods %}
   {% set filename ="{method}/index.rst".format(method=method) %}
   {{ filename|indent(3, False) }}
   {% endfor %}
