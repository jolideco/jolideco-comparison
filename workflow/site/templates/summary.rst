{{ title }}
{{ "=" * title|length }}


Configuration
-------------
This is the configuration the deconvolution was run with:

::

    {{ configuration|indent(4, False) }}


{% if model_configuration %}
Model Configuration
-------------------
This is the configuration the Jolideco model:

::

    {{ model_configuration|indent(4, False) }}

{% endif %}


Flux
----


.. figure:: images/flux.png
    :width: 1000

    The plot shows the deconvolved flux image, the reference flux and residuals.


Predicted Counts
----------------

.. figure:: images/npred.png
    :width: 1000

    The plot shows the predicted counts correspoding to the deconvolved flux image,
    the reference flux and residuals.


Metrics
-------
MSE = Mean Squared Error, SSI = Structural Similarity Index

.. list-table:: 
   :header-rows: 1

   * - Metric
     - Value

   {% for name, value in metrics.items() %}
   * - {{ name }}
     - {{ value }}
   {% endfor %}




Traces
------

.. figure:: images/trace.png
    :width: 1000


Files
-----

Results files for download:

:download:`{{ filename_result }}`