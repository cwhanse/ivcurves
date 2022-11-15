// Disable ordering, pagination, and search by default.
// This matches the functionality of a normal Sphinx table.
$.extend($.fn.dataTable.defaults, {
    ordering: false,
    paging:  false,
    searching: false
});

$(document).ready( function () {
    // Enable ordering for leaderboard and compare-submissions tables
    $('#leaderboard > table').DataTable({ ordering: true });
    $('#compare-submissions > table').DataTable({ ordering: true });
    // DataTable does not create a DataTable for elements that are already
    // a DataTable, so trying to create a DataTable for all tables with this
    // following line is ok. If this is not the case anymore, swap this line
    // with the previous two lines.
    $('table').DataTable();
} );

