// Main jQuery call
$(document)
    .ready(
	function() {

	    $.ajax({	
		data : {
		    sessionid: JSON.stringify(sessionid),
		},
		type : 'POST',
		url : '/msa',
		success : function(serverdata) {
		    var data = JSON.parse(serverdata);
		    var seqs = msa.io.fasta.parse(data.msa);
		    var opts = {
			el: document.getElementById("msa"),
			seqs: seqs,
			vis: {
			    conserv: false,
			    overviewbox: true,
			    seqlogo: true
			},
			zoomer : {
			    rowHeight: 20,
			},
			// smaller menu for JSBin
			menu: "small",
			bootstrapMenu: true
		    };
		    var m = msa(opts);
		    m.render();
		},
		error : function() {
	    console.log('Sorry, no luck');
		}
	    });
		
});
// To do: display tree 
