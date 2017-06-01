

function addFilter() {
    columns = [];
    $('.dataframe thead th').each(
			 function( index, value ) {
			     columns.push(value.innerHTML);
			 }

	);
    var sel = $('.Filter').append('select');
    sel.addClass('sel1');
    for (x = 0; x < columns.length; x++) {
	$('.sel1').append('option').attr(value,x).text(columns[x]);
    }
    return columns;
}

// Main jQuery call
$(document)
    .ready(
	function() {

	    addFilter();

	    $('#dataframe').addClass('Selenzy');

	    

	}
    );
