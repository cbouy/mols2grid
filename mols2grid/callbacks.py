from .utils import env


def make_popup_callback(title, html, js="", style=""):
    """Creates a JavaScript callback that displays a popup window

    Parameters
    ----------
    title : str
        Title of the popup. Use ``title='${data["Name"]}'`` to use the value
        of the column "Name" as a title
    html : str
        Content of the popup window
    js : str
        JavaScript code executed before making the content of the popup window.
        This allows you to create variables and reuse them later in the `html`
        content of the popup, using the ``${my_variable}`` syntax
    style : str
        CSS style assigned to the popup window

    Returns
    -------
    js_callback : str
        JavaScript code that allows to display a popup window
    """
    return (env.get_template('js/popup.js')
               .render(js=js,
                       html=html,
                       title=title,
                       style=style))


def info(title="SMILES", img_size=(400, 300), style="max-width: 80%;"):
    """Displays a bigger image of the molecule, alongside some useful descriptors:
    molecular weight, number of Hydrogen bond donors and acceptors, TPSA and Crippen
    ClogP.

    Parameters
    ----------
    title : str or None
        Data field used to set the title of the modal window. If ``None``, no title is
        displayed.
    img_size : tuple
        Width and height of the molecule depiction.
    style : str
        CSS style applied to the modal window.
    """
    title_field = "${data['" + title + "']}" if title else None
    return make_popup_callback(
        title=title_field,
        js=f"""
            var mol = RDKit.get_mol(data["SMILES"]);
            var svg = mol.get_svg({img_size[0]}, {img_size[1]});
            var desc = JSON.parse(mol.get_descriptors());
            mol.delete();
        """,
        html="""
            <div class="row">
            <div class="col">${svg}</div>
            <div class="col">
                <b>Molecular weight</b>: ${desc.exactmw}<br/>
                <b>HBond Acceptors</b>: ${desc.NumHBA}<br/>
                <b>HBond Donors</b>: ${desc.NumHBD}<br/>
                <b>TPSA</b>: ${desc.tpsa}<br/>
                <b>ClogP</b>: ${desc.CrippenClogP}<br/>
            </div>
            </div>""",
        style=style,
    )