function toggleVisibility(obj) {
    if ( obj.css('visibility') == 'hidden' ) 
	obj.css('visibility', 'visible');
    else
	obj.css('visibility', 'hidden');
}

function displayReaction(rxninfo) {
    // display reaction
    toggleVisibility( $('.Info1') );
    $('#front').submit();
}

function rxnIO() {
     $('#rxn').change( function() {
	 if ($('#rxn').val() != '') {
	     $('#smarts').val('');
	     displayReaction(this.value);
	     }
     });
}

function smartsIO() {
    $('#smarts').change( function() {
	if ($('#smarts').val() != '') {
	    $('#rxn').val('');
	    displayReaction(this.value);
	   }
    });
}


// Main jQuery call
$(document)
    .ready(
	function() {
	    rxnIO();
	    smartsIO();
	    $('#options').submit(
		function( event ) {
		    toggleVisibility( $('.Info2') );
		    });
	   $('#front').submit(
	    function( event ) {
		$.ajax({	
		    data : new FormData( this ),
		    processData: false,
		    contentType: false,
		    type : 'POST',
		    url : '/display',
		    success : function(serverdata) {
			data = JSON.parse(serverdata);
			if (data['status'] == true) {
			    $('#upload').prop('disabled', false);
			}
			if (data['success'] == true) {
			    $('.canvas img').attr('src', data['data']);
//			    $('.canvas').html(data['data']);
//			    svg = $('svg')[0];
//			    var bbox = svg.getBBox();
//			    var viewBox = [bbox.x, bbox.y, bbox.width, bbox.height].join(" ");
//			    svg.setAttribute("viewBox", viewBox);
//			    svg.removeAttribute('height');
			    $('.canvas').dialog({width: 600, title: 'Reaction query'});
			    toggleVisibility( $('.Info1') );
			}
		    },
		    error : function() {
			alert('Unknown reaction format');
			console.log('Sorry, no luck');
		    }
		});
		event.preventDefault();
	    }
	   );
	    
	});
