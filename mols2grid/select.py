selection = {}

def reset_selection():
    global selection
    selection.clear()

def set_selection(_id, smiles):
    global selection
    selection[_id] = smiles

def del_selection(_id):
    global selection
    del selection[_id]

try:
    from google import colab
except (ModuleNotFoundError, ImportError):
    pass
else:
    colab.output.register_callback('m2g.reset_selection', reset_selection)
    colab.output.register_callback('m2g.set_selection', set_selection)
    colab.output.register_callback('m2g.del_selection', del_selection)
