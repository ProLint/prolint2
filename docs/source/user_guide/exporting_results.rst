Exporting the results
=====================
You can export contacts information using the **export** method as below:

.. code-block:: python

    target_system.contacts.export('results.csv') 
    # this will export two csv files, one with the contacts information 
    # stored in target_system.contacts.contacts and a second one 
    # ('results_metrics.csv') with the metrics information.