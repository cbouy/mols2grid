// prerequisite JavaScript code
{{ js }}
// HTML template for the popup
var html = `
<div class="modal fade" id="m2g-modal" tabindex="-1">
  <div class="modal-dialog" style="{{ style }}">
    <div class="modal-content">
      <div class="modal-header">
        <h5 class="modal-title">{{ title }}</h5>
        <button type="button" class="close" data-dismiss="modal">
          <span>&times;</span>
        </button>
      </div>
      <div class="modal-body">
      <!-- provided HTML code -->
        {{ html }}
      <!-- end of provided code -->
      </div>
    </div>
  </div>
</div>
`
// create container element where the popup will be inserted
if ($("#modal-container").length === 0) {
  $('<span id="modal-container"></span>').insertAfter("#mols2grid");
}
// insert the code inside the container element
$('#modal-container').html(html);
// trigger display of popup window
$('#m2g-modal').modal('show'); 