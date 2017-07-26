

function addFilter(nFilters) {
    var columns = [];
    $('.dataframe thead th').each(
			 function( index, value ) {
			     columns.push(value.innerHTML);
			 }

	);
    var selectors = []
    for (i=0; i < nFilters; i++) {
	selectors[i] = $('<select>' ).addClass('prior').attr('id', 'sel'+parseInt(i));
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
	    $('.dataframe').addClass('Selenzy');
	    addArrows();
	    addLinks();
	},
	error : function() {
	    console.log('Sorry, no luck');
	}
    });
}

function addArrows() {
    $('thead th').each(
	function( index, value ) {
	    if (index > 0 ){
		value.innerHTML = '<div class="cell">' +
		    '<div class="head">' + value.innerHTML+'</div>'+
		    '<div class="uarr" index="'+index.toString()+'">' 
		    + '&uarr;' + '</div> &nbsp; &nbsp;' +
		    '<div class="darr" index="'+index.toString()+'">' 
		    + '&darr;' + '</div>'
		    + '<div class="footcell"> </div> </div>';
	    }
	}	
    );

    $('.uarr').click( function () {
	sortTable([parseInt(this.getAttribute("index"))]);
    }
		    );
    $('.darr').click( function () {
	sortTable([-parseInt(this.getAttribute("index"))]);
    }
		    );


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





// Main jQuery call
$(document)
    .ready(
	function() {

	    $('.dataframe').addClass('Selenzy');

	    addLinks();

	    selectors = addFilter(3);

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

	    addArrows();


	}
    );
