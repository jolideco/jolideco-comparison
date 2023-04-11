{{ title }}
{{ "=" * title|length }}


Results Overview
----------------
.. list-table:: 
   :header-rows: 1

   * - 
     - Ground Truth
     - Counts
     {% for method in methods %}
     - {{ method }}
     {% endfor %}

   * - Flux
     - .. image:: ../../images/flux-true-thumbnail.png
          :target: ../../images/flux-true-thumbnail.png
          :align: center
          :width: 200px

     - .. image:: images/counts-thumbnail.png
          :target: images/counts-thumbnail.png
          :align: center
          :width: 200px

     {% for method in methods %}
     - .. image:: {{ "{method_}/images/flux-thumbnail.png".format(method_=method) }}
          :target: {{ "{method_}/index.html".format(method_=method) }}
          :align: center
          :width: 200px
     {% endfor %}

   * - SSI
     - 1
     - 
     {% for method in methods %}
     - {{ metrics[method]["SSI"] }}
     {% endfor %}

   * - MSE
     - 0
     - 
     {% for method in methods %}
     - {{ metrics[method]["MSE"] }}
     {% endfor %}

   * - NRMSE
     - 0
     - 
     {% for method in methods %}
     - {{ metrics[method]["NRMSE"] }}
     {% endfor %}

   * - NMI
     - 0
     - 
     {% for method in methods %}
     - {{ metrics[method]["NMI"] }}
     {% endfor %}


MSE = Mean Squared Error,
SSI = Structural Similarity Index,
NMI = Normalized Mutual Information,
NRMSE = Normalized Root Mean Squared Error


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
   :maxdepth: 1
   :hidden:

   {% for method in methods %}
   {% set filename ="{method}/index.rst".format(method=method) %}
   {{ filename|indent(3, False) }}
   {% endfor %}
