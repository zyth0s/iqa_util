# iqa_util
Some scripts to analyze IQA results. IQAJSON, we shall call it, is an approach to store
all IQA information using JSON format.

Start IPython and execute `%run iqanalysis.py`. A window will open
requesting you to select the IQAJSON file. Then, you will return back to the
interpreter. A new variable `m` will be defined, that contains all extracted
information. `m.TAB` will show you all accessible properties or methods.

# Dependencies

* json package
* matplotlib
* numpy
* PyQt5 (not strictly necessary but better to have it)
