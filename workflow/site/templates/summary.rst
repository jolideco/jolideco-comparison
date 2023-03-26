{{ title }}
{{ "=" * title|length }}


Configuration
-------------
This is the configuration the deconvolution was run with:

::

    {{ configuration|indent(4, False) }}


Flux
----


.. image:: images/flux.png
    :width: 600

    The plot shows the deconvolved flux image, the reference flux and residuals.


Predicted Counts
----------------

.. image:: images/npred.png
    :width: 600



Trace
-----

.. image:: images/trace.png
    :width: 600


Files
-----

Results files for download:

:download:`{{ filename_result }}`