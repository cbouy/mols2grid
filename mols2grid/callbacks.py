from mols2grid.utils import env


def make_popup_callback(
    title=None, subtitle=None, svg=None, html=None, js=None, style=None, identifier=None
):
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
    identifier : str
        Grid element identifier

    Returns
    -------
    js_callback : str
        JavaScript code that allows to display a popup window
    """
    return env.get_template("js/callbacks/popup.js.j2").render(
        title=title,
        subtitle=subtitle,
        html=html or "",
        svg=svg,
        js=js or "",
        style=style or "",
        identifier=f"#{identifier}" or "",
    )


def _get_title_field(title):
    return "${data['" + title + "']}" if title else None


def info(
    title="SMILES", subtitle=None, img_size=(-1, -1), style=None, scaling_factor=25
) -> str:
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
    scaling_factor : int
        Scaling factor used for the drawings.
    """
    return make_popup_callback(
        title=_get_title_field(title),
        subtitle=_get_title_field(subtitle),
        svg="${svg}",
        js=f"""
            let mol = RDKit.get_mol(data["SMILES"]);
            let options = {{
                width: {img_size[0]},
                height: {img_size[0]},
                scalingFactor: {scaling_factor}
            }}
            let svg = mol.get_svg_with_highlights(JSON.stringify(options));
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


def show_3d(
    title="SMILES",
    subtitle=None,
    query=None,
    height="100%",
    style="width:100%; height:100%",
) -> str:
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
    if query is None:
        query = ["pubchem", "cactus"]
    js_script = env.get_template("js/callbacks/show_3d.js.j2").render(query=query)
    return make_popup_callback(
        title=_get_title_field(title),
        subtitle=_get_title_field(subtitle),
        html=f'<div id="molviewer" style="width: 100%; height: {height};"></div>',
        js=js_script,
        style=style,
    )


def external_link(
    url="https://leruli.com/search/{}/home",
    field="SMILES",
    url_encode=False,
    b64_encode=True,
) -> str:
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
    return env.get_template("js/callbacks/external_link.js.j2").render(
        url=url,
        field=field,
        url_encode=url_encode,
        b64_encode=b64_encode,
    )
