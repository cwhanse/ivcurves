// Disable ordering, pagination, and search by default.
// This matches the functionality of a normal Sphinx table.
$.extend($.fn.dataTable.defaults, {
    ordering: false,
    paging:  false,
    searching: false,
    bInfo: false
});

$(document).ready(() => {
    // Enable ordering compare-submissions table
    cs_t = $('#compare-submissions > table').DataTable({ ordering: true });
    // Sort the compare-submissions table by the first column
    cs_t.column(1).order('asc').draw();
});

