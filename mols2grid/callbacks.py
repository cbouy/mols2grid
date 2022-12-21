from typing import Optional, NamedTuple

from .utils import env


class _JSCallback(NamedTuple):
    code: str
    library_src: Optional[str] = None


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


def _get_title_field(title):
    return "${data['" + title + "']}" if title else None


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
    return make_popup_callback(
        title=_get_title_field(title),
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


def show_3d(
        title="SMILES",
        query=["pubchem", "cactus"],
        height="350px",
        style="max-width: 80%"
    ):
    """Queries the API listed in ``resolver`` using the SMILES of the structure, to
    fetch the 3D structure and display it with 3Dmol.js

    Parameters
    ----------
    title : str or None
        Data field used to set the title of the modal window. If ``None``, no title is
        displayed.
    query : list or dict
        List of APIs used to fetch the 3D structure from (by order of priority). To use
        a custom API, use a dict with the following format::

            {
                "url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{}/SDF?record_type=3d",
                "format": "sdf",
                "field": "SMILES",
                "encode": true,
            }

        In this example, the value in the ``SMILES`` field will be URI-encoded and
        replace the curly-braces in the ``url`` strings to create the URL to fetch
        the 3D structure from in the SDF format.

    height : str
        Height of the 3Dmol.js viewer in the modal.
    style : str
        CSS style applied to the modal window.
    """
    js_script = env.get_template('js/callbacks/show_3d.js').render(query=query)
    code = make_popup_callback(
        title=_get_title_field(title),
        js=js_script,
        html=f'<div id="molviewer" style="width: 100%; height: {height};"></div>',
        style=style,
    )
    library = """<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.8.0/3Dmol-nojquery-min.js" integrity="sha512-9iiTgstim185ZZPL/nZ+t+MLMmIbZEMfoZ1swSBUhxt4AukOPY34yyO2217X1dN5ziVMKi+YLmp/JBj+KyEaUQ==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>"""
    return _JSCallback(code=code, library_src=library)


def external_link(
    url="https://leruli.com/search/{}/home",
    field="SMILES",
    url_encode=False,
    b64_encode=True,
):
    if url_encode and b64_encode:
        raise ValueError("Setting both URL and B64 encoding is not supported")
    return env.get_template('js/callbacks/external_link.js').render(
        url=url,
        field=field,
        url_encode=url_encode,
        b64_encode=b64_encode,
    )
