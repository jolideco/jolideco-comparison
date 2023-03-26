{{ title }}
{{ "=" * title|length }}


Configuration
-------------
This is the configuration the deconvolution was run with:

::

    {{ configuration|indent(4, False) }}


Flux
----


.. figure:: images/flux.png
    :width: 600

    The plot shows the deconvolved flux image, the reference flux and residuals.


Predicted Counts
----------------

.. figure:: images/npred.png
    :width: 600

    The plot shows the predicted counts correspoding to the deconvolved flux image,
    the reference flux and residuals.


Trace
-----

.. figure:: images/trace.png
    :width: 600


Files
-----

Results files for download:

:download:`{{ filename_result }}`