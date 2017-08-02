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

function customAlert( message, title ) {
    if (!title ) 
	title = 'Alert';
    if ( !message ) 
	message = 'System error';
    $( '<div></div>' ).html( message ).dialog( {
	title: title,
	resizable: false,
	modal: true,
	buttons: {
	    'Ok': function() {
		$( this ).dialog( 'close' )
	    }
	   }
	});
}

function resetForm( message ) {
    customAlert( message );
    $('.Info1').css('visibility', 'hidden' ) ;
    $('.Info2').css('visibility', 'hidden' ) ;
    $('#upload').prop('disabled', true);
    $('#rxn').val('');
    $('#smarts').val('');
    if ( $('.canvas').hasClass('ui-dialog-content') ) {
	$('.canvas').dialog('close');
    }
    $( '.status' ).css('color', 'red').html('No reaction.');

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
			toggleVisibility( $('.Info1') );
			data = JSON.parse(serverdata);
			if (data['success'] == true) {
			    $('#upload').prop('disabled', false);
			    if (data['svg'] == false) {
				$('.canvas img').attr('src', data['data']).attr('width', 600);
			    } else {
				$('.canvas').html(data['data']);
				svg = $('svg')[0];
				var bbox = svg.getBBox();
				var viewBox = [bbox.x, bbox.y, bbox.width, bbox.height].join(" ");
				svg.setAttribute("viewBox", viewBox);
				svg.removeAttribute('height');
			    }
			    $('.canvas').dialog({
				title: 'Reaction query',
				width: 640,
				height: data['size'][1],
				resize: function(event, ui) {
				    $('.canvas img').attr('width', $('.canvas').width() );
				    }
			    });
			    $( '.status' ).css('color', 'green').html('Validated.');
			} else {
			    resetForm( 'Unknown reaction format!<br/>Supported formats: SMILES, SMARTS, MDL RXN.' );
			}
		    },
		    error : function() {
			resetForm( 'Bad request' );
		    }
		});
		event.preventDefault();
	    }
	   );
	    
	});
