function customConfirm( message, title ) {
    if (!title ) 
	title = 'Alert';
    if ( !message ) 
	message = 'System error';
    $( '<div></div>' ).html( message ).dialog( {
	title: title,
	resizable: false,
	modal: true,
	buttons: {
	    'Confirm': function() {
		deleteRows();
		$( this ).dialog( 'close' )
	    },
	    'Cancel': function() {
		$( this ).dialog( 'close' )
	    }
	   }
	});
}


function formatTable(csvdata) {
    $('.dataframe').addClass('Selenzy');
    $('.Selenzy th').each( function(index, value) {
	$( this ).addClass('active');
	});

    addSorted();
    addLinks();
    addSelection();
    $( '.Remove' ).addClass('disabledbutton');
    $( '.Remove').click( function(event) {
	if (selectedRows().length > 0) {
	    customConfirm('Do you want to remove permanently selected rows?');
	  }
    });
}

function formatHeader(filt) {
    for (i=0; i < filt.length; i++) {
	xval = filt[i];
	if (xval > 0) {
	    $( $('thead th')[xval]).attr("sort", "up").attr("title", "Sort by column");
	   }
    }
}

function addFilter(nFilters) {
    var columns = [];
    $('.dataframe thead th').each(
			 function( index, value ) {
			     columns.push(value.innerHTML);
			 }

	);
    var selectors = []
    for (i=0; i < nFilters; i++) {
	selectors[i] = $('<select>' ).addClass('prior item').attr('id', 'sel'+parseInt(i));
	for (x = 0; x < columns.length; x++) {
	    selectors[i].append( $('<option>', {value:x, text:columns[x]}) );
	}
	$('.Filter').append(selectors[i]);
    }
    return selectors;
}

function sortTable(valueSelected) {
    // post and retrieve data from server to sort table
    $.ajax({	
	data : {
	    filter : JSON.stringify(valueSelected),
	    csv : JSON.stringify(csvlink),
	    session: JSON.stringify(sessionid),
	    event: JSON.stringify(event.timeStamp)
	},
	type : 'POST',
	url : '/sorter',
	success : function(serverdata) {
	    data = JSON.parse(serverdata);
	    $('.Selenzy').replaceWith(data.data.csv);
	    formatTable();
	    formatHeader(data.data.filter);
	},
	error : function() {
	    console.log('Sorry, no luck');
	}
    });
}

function deleteRows() {
    // Delete selected rows
    $.ajax({	
	data : {
	    filter : JSON.stringify(selectedRows()),
	    csv : JSON.stringify(csvlink),
	    session: JSON.stringify(sessionid),
	    event: JSON.stringify(event.timeStamp)
	},
	type : 'POST',
	url : '/remover',
	success : function(serverdata) {
	    data = JSON.parse(serverdata);
	    $('.Selenzy').replaceWith(data.data.csv);
	    // Avoid propagation of the click event
	    $( '.Remove').off("click");
	    formatTable();

	},
	error : function() {
	    console.log('Sorry, no luck');
	}
    });
}

function selectedRows() {
    // Returns an array with the indices of selected rows
    var index = [];
    $.each( $( 'tbody th.selected'), function( i, value) {
	index[i] = $( value ).text();
    });
    return(index);

}

function addSorted() {
       $('thead th').each(
	   function( index, value) {
	       $( this ).attr( "index", index.toString());
	       $( this ).attr( "sort", "down");
	       $( this ).attr( "title", "Sort by column");
	       $( this ).click(
		   function() {
		       if (this.getAttribute("sort") == "up") {
			   $( this ).attr("sort", "down");
			   sortTable([-parseInt(this.getAttribute("index"))]);
		      } else {
			   $( this ).attr("sort", "up");
			   sortTable([parseInt(this.getAttribute("index"))]);
		      }
		   }
		   );
	   }
	   );
}

function addNavigation() {
    var csvTag = $( '<a>' ).attr('href', csvlink)
	.addClass('item').append( "[Download CSV]");
    var fastaTag = $( '<a>' ).attr('href', fastalink).attr('target', '_blank')
	.addClass('item').append("[Download FASTA]");
    var msaTag = $( '<a>' ).attr('href', msafastalink).attr('target', '_blank')
	.addClass('item msa').append("[Download MSA]");
    var msaViewTag = $( '<a>' ).attr('href', msaviewlink).attr('target', '_blank')
	.addClass('item msa').append("[View MSA]");
    if (flagFasta == "False") {
	fastaTag.addClass('disabledbutton');
    }
    if (flagMSA == "False") {
	msaTag.addClass('disabledbutton');
	msaViewTag.addClass('disabledbutton');
    }

    $( '.Navigate' ).append(csvTag).append(fastaTag).append(msaTag).append(msaViewTag);

}


function addLinks() {
    var rows = document.getElementsByTagName('tr');
    for (var i = 1; i < rows.length; i++)	{

	var seqID = rows[i].getElementsByTagName('td')[0];
	var rxnID = rows[i].getElementsByTagName('td')[4];

	var seqlink = "//www.uniprot.org/uniprot/" + rows[i].getElementsByTagName('td')[0].innerHTML;

	var rxnlink = "http://www.metanetx.org/cgi-bin/mnxweb/equa_info?equa=" + rows[i].getElementsByTagName('td')[4].innerHTML;

	var aTag = document.createElement('a');
	var bTag = document.createElement('a');
	

	aTag.setAttribute('href', seqlink);
	aTag.setAttribute('target', '_blank');
	bTag.setAttribute('href', rxnlink);
	bTag.setAttribute('target', '_blank');

	aTag.innerHTML = rows[i].getElementsByTagName('td')[0].innerHTML;
	bTag.innerHTML = rows[i].getElementsByTagName('td')[4].innerHTML;

	seqID.replaceChild(aTag, seqID.childNodes[0]);
	rxnID.replaceChild(bTag, rxnID.childNodes[0]);

    }
}


function addSelection() {
    $('.Selenzy tr th').click(
	function(event) {
	    if ($( this ).hasClass('selected')) {
		$( this ).removeClass('selected');
		$( this ).parent().removeClass('selected');
		if ( selectedRows().length == 0 ) {
		    $( '.Remove' ).addClass('disabledbutton');
		}
	    } else {
		$( this ).addClass('selected');
		$( this ).parent().addClass('selected');			
		$( '.Remove' ).removeClass('disabledbutton');
	    }
	}
    );
    
    $('.Selenzy tbody tr th').attr('title', 'Select row');

}


// Main jQuery call
$(document)
    .ready(
	function() {
	    addNavigation();

	    formatTable();

	    // selectors = addFilter(3);
	    // $('.prior').selectmenu({
	    // 	change: function(event) {
	    // 	    var optionsFilter = [];
	    // 	    for (i=0; i < selectors.length; i++) {
	    // 		var optionSelected = $("option:selected", selectors[i][0]);
	    // 		var valueSelected = selectors[i][0].value;
	    // 		optionsFilter[i] = valueSelected;
	    // 	    }
	    // 	    sortTable(optionsFilter);
	    // 	    }
	    // 	});


	}
    );
