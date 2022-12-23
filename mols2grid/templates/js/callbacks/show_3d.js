// fetch file and display with 3Dmol.js
let resolvers = {
    "pubchem": {
        "url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{}/SDF?record_type=3d",
        "format": "sdf",
        "field": "SMILES",
        "encode": true,
    },
    "cactus": {
        "url": "https://cactus.nci.nih.gov/chemical/structure/{}/file?format=sdf",
        "format": "sdf",
        "field": "SMILES",
        "encode": true,
    }
}
function get_url(resolver, data) {
    let value = data[resolver.field];
    if (resolver.encode) {
        value = encodeURIComponent(value)
    }
    return resolver.url.replace("{}", value);
}
// fetch and display
function show_3d(data, apis_or_custom_resolver, viewer) {
    let resolver, api;
    if (Array.isArray(apis_or_custom_resolver)) {
        // go through list of APIs available
        api = apis_or_custom_resolver.shift();
        if (api === undefined) {
            console.error("No API left to query...")
            return
        }
        resolver = resolvers[api];
    } else {
        // user passed a custom resolver
        api = "custom"
        resolver = apis_or_custom_resolver;
    }
    let sdf_url = get_url(resolver, data);
    $.ajax(sdf_url, { 
        success: function(data) {
            viewer.addModel(data, resolver.format);
            viewer.setStyle({}, {stick: {}});
            viewer.zoomTo();
            viewer.setHoverable(
                {}, 1,
                function(atom, viewer, event, container) {
                    if (!atom.label) {
                        atom.label = viewer.addLabel(
                            atom.serial + ':' + atom.atom, {
                                position: atom,
                                backgroundColor: 'mintcream',
                                fontColor:'black'
                            }
                        );
                    }
                },
                function(atom, viewer) { 
                    if (atom.label) {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                    }
                }
            );
            viewer.render();
        },
        error: function(hdr, status, err) {
            console.error(
                "Failed to load SDF with " + api
            );
            show_3d(data, apis_or_custom_resolver, viewer);
        },
    });
}
$(document).ready(function() {
    // 3Dmol.js options
    let element = $('#molviewer');
    let config = { backgroundColor: 'white' };
    let viewer = $3Dmol.createViewer(element, config);
    // prepare query to fetch 3D SDF from SMILES
    let apis_or_custom_resolver = {{ query | tojson }};
    show_3d(data, apis_or_custom_resolver, viewer);
});