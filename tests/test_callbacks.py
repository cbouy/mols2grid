from mols2grid import callbacks

def test_make_popup_callback():
    popup = callbacks.make_popup_callback(
        title="${title}",
        html='<span id="foo">${title}</span>',
        js='var title = "FOOBAR";',
        style="max-width: 42%;",
    )
    assert '<h5 class="modal-title">${title}</h5>' in popup
    assert '<span id="foo">${title}</span>' in popup
    assert '// prerequisite JavaScript code\nvar title = "FOOBAR";' in popup
    assert '<div class="modal-dialog" style="max-width: 42%;">' in popup