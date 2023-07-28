from typing import NamedTuple, Optional

from .utils import env


class _JSCallback(NamedTuple):
    """Class that holds JavaScript code for running a callback function. If an external
    library is required for the callback to function correctly, it can be passed in
    the optional ``library_src`` as a ``<script>`` tag.
    """

    code: str
    library_src: Optional[str] = None


def make_popup_callback(title=None, subtitle=None, svg=None, html="", js="", style=""):
    """Creates a JavaScript callback that displays a popup window

    Parameters
    ----------
    title : str
        Title of the popup. Use ``title='${data["Name"]}'`` to use the value
        of the column "Name" as a title
    subtitle : str
        Secondary title which works the same as title.
    svg : str
        SVG depiction of the molecule
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
    return env.get_template("js/popup.js").render(
        title=title, subtitle=subtitle, html=html, svg=svg, js=js, style=style
    )


def _get_title_field(title):
    return "${data['" + title + "']}" if title else None


def info(title="SMILES", subtitle=None, img_size=(400, 300), style="") -> _JSCallback:
    """Displays a bigger image of the molecule, alongside some useful descriptors:
    molecular weight, number of Hydrogen bond donors and acceptors, TPSA, Crippen ClogP
    and InChIKey.

    Parameters
    ----------
    title : str
        Title of the popup. Use ``title='${data["Name"]}'`` to use the value
        of the column "Name" as a title
    subtitle : str
        Secondary title which works the same as title.
    img_size : tuple
        Width and height of the molecule depiction.
    style : str
        CSS style applied to the modal window.
    """
    code = make_popup_callback(
        title=_get_title_field(title),
        subtitle=_get_title_field(subtitle),
        svg="${svg}",
        js=f"""
            let mol = RDKit.get_mol(data["SMILES"]);
            let svg = mol.get_svg({img_size[0]}, {img_size[1]});
            let desc = JSON.parse(mol.get_descriptors());
            let inchikey = RDKit.get_inchikey_for_inchi(mol.get_inchi());
            mol.delete();
        """,
        html="""
            <b>Molecular weight</b>: ${desc.exactmw}<br/>
            <b>HBond Acceptors</b>: ${desc.NumHBA}<br/>
            <b>HBond Donors</b>: ${desc.NumHBD}<br/>
            <b>TPSA</b>: ${desc.tpsa}<br/>
            <b>ClogP</b>: ${desc.CrippenClogP}<br/>
            <hr>
            <b>InChIKey</b>: ${inchikey}
            """,
        style=style,
    )
    return _JSCallback(code=code)


def show_3d(
    title="SMILES",
    subtitle=None,
    query=["pubchem", "cactus"],
    height="100%",
    style="width:100%;height:100%",
) -> _JSCallback:
    """Queries the API(s) listed in ``query`` using the SMILES of the structure, to
    fetch the 3D structure and display it with ``3Dmol.js``

    Parameters
    ----------
    title : str
        Title of the popup. Use ``title='${data["Name"]}'`` to use the value
        of the column "Name" as a title
    subtitle : str
        Secondary title which works the same as title.
    query : list or dict
        List of APIs used to fetch the 3D structure from (by order of priority). To use
        a custom API, use a dict with the following format::

            {
                "url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{}/SDF?record_type=3d",
                "format": "sdf",
                "field": "SMILES",
                "encode": True,
            }

        In this example, the value in the ``SMILES`` field will be URI-encoded and
        replace the curly-braces in the ``url`` strings to create the URL to fetch
        the 3D structure from in the SDF format.
    height : str
        Height of the 3Dmol.js viewer in the modal.
    style : str
        CSS style applied to the modal window.
    """
    js_script = env.get_template("js/callbacks/show_3d.js").render(query=query)
    code = make_popup_callback(
        title=_get_title_field(title),
        subtitle=_get_title_field(subtitle),
        html=f'<div id="molviewer" style="width: 100%; height: {height};"></div>',
        js=js_script,
        style=style,
    )
    library = """<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.8.0/3Dmol-nojquery-min.js" integrity="sha512-9iiTgstim185ZZPL/nZ+t+MLMmIbZEMfoZ1swSBUhxt4AukOPY34yyO2217X1dN5ziVMKi+YLmp/JBj+KyEaUQ==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>"""
    return _JSCallback(code=code, library_src=library)


def external_link(
    url="https://leruli.com/search/{}/home",
    field="SMILES",
    url_encode=False,
    b64_encode=True,
) -> _JSCallback:
    """Opens an external link using ``url`` as a template string and the value in the
    corresponding ``field``. The value can be URL-encoded or base64-encoded if needed.

    Parameters
    ----------
    url : str
        Template string used to generate the URL that will be opened.
    field : str
        Field name used to generate the URL that will be opened. The value can be
        encoded (see below).
    url_encode : bool
        Encode the value fetched from the specified field to replace characters that are
        not allowed in a URL e.g., spaces become ``%20``.
    b64_encode : bool
        Base64-encode the value fetched from the field.

    Raises
    ------
    ValueError : Both ``url_encode`` and ``b64_encode`` have been specified.
    """
    if url_encode and b64_encode:
        raise ValueError("Setting both URL and B64 encoding is not supported")
    code = env.get_template("js/callbacks/external_link.js").render(
        url=url,
        field=field,
        url_encode=url_encode,
        b64_encode=b64_encode,
    )
    return _JSCallback(code=code)
