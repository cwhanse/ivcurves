// Disable ordering, pagination, and search by default.
// This matches the functionality of a normal Sphinx table.
$.extend($.fn.dataTable.defaults, {
    ordering: false,
    paging:  false,
    searching: false,
    bInfo: false
});

$(document).ready(() => {
    // Enable ordering scoreboard table
    cs_t = $('#scoreboard > table').DataTable({ ordering: true });
    // Sort the scoreboard table by Overall Score
    cs_t.column(2).order('asc').draw();
});

