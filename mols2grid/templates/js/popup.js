// Prerequisite JavaScript code.
// prettier-ignore
{{ js }}

// HTML template for the popup.
var html = `
<div id="m2g-modal" tabindex="-1" style="{{ style }}">

    <div class="m2g-modal-header">
        {% if title %}<h2>{{ title }}</h2>{% endif %}
        {% if subtitle %}<p>{{ subtitle }}</p>{% endif %}
        <button class="close">&times;</button>
    </div>

    <div class="m2g-modal-body">
        {% if svg %}<div class="svg-wrap">{{ svg }}</div>{% endif %}
        {{ html }}
    </div>

</div>
`
// Create container element where the popup will be inserted.
if ($('#m2g-modal-container').length === 0) {
    $('<span id="m2g-modal-container"></span>').insertAfter('#mols2grid')
}

// Insert the code inside the container element.
$('#m2g-modal-container').html(html)

// Show modal.
setTimeout(function () {
    $('#m2g-modal-container').addClass('show')
}, 0)

// Hide modal on close / ESC key.
$('#m2g-modal-container').click(function (e) {
    if (e.target.id == 'm2g-modal-container' || e.target.className == 'close') {
        closeModal()
    }
})
$(document).keydown(function (e) {
    if (e.key == 'Escape') {
        closeModal()
        e.preventDefault()
    }
})
function closeModal() {
    $('#m2g-modal-container').removeClass('show')
    setTimeout(function () {
        $('#m2g-modal-container').remove()
    }, 150 + 10)
}
