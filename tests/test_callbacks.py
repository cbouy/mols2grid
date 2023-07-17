import pytest

from mols2grid import callbacks


def test_make_popup_callback():
    popup = callbacks.make_popup_callback(
        title="${title}",
        subtitle="${title}",
        svg="<svg><rect width='10' height='10'/></svg>",
        html='<span id="foo">${title}</span>',
        js='var title = "FOOBAR";',
        style="max-width: 42%;",
    )
    assert '<h2>${title}</h2>' in popup
    assert '<p>${title}</p>' in popup
    assert '<span id="foo">${title}</span>' in popup
    assert '// Prerequisite JavaScript code.\n// prettier-ignore\nvar title = "FOOBAR";' in popup
    assert '<div id="m2g-modal" tabindex="-1" style="max-width: 42%;">' in popup


@pytest.mark.parametrize(
    "title, expected",
    [
        ("SMILES", "${data['SMILES']}"),
        (None, None),
    ],
)
def test_title_field(title, expected):
    assert callbacks._get_title_field(title) == expected
