const field = "{{ field }}";
let value = data[field];
let url = "{{ url }}";

{% if url_encode %}
value = encodeURIComponent(value);
{% elif b64_encode %}
value = window.btoa(unescape(encodeURIComponent(value)));
{% endif %}

url = url.replace("{}", value)
window.open(url, '_blank').focus();