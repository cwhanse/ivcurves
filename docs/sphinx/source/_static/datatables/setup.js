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
    $('#compare-submissions > table').DataTable({ ordering: true });
});

