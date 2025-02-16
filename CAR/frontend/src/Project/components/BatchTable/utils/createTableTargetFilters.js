/**
 * Creates a table filter. Checks only for the target (parent) row and returns all of its children
 */
const filterTableTargetRows = (rowFilter, rows, ids, filterValue) => {
  return (
    rows
      .filter(row => {
        if (row.depth === 0) {
          return rowFilter(row, ids, filterValue);
        }

        return true;
      })
      // A bug in the library results in loosing subRows when applying filtering. To avoid it a copy of each row
      // has to be returned from the filter method.
      .map(row => ({ ...row }))
  );
};

export const createTableTargetAutocompleteFilter = rowFilter => (rows, ids, filterValue) => {
  if (!filterValue || !filterValue.length) {
    return [...rows];
  }

  return filterTableTargetRows(rowFilter, rows, ids, filterValue);
};

export const createTableTargetYesNoFilter = rowFilter => (rows, ids, filterValue) => {
  if (filterValue === '') {
    return [...rows];
  }

  return filterTableTargetRows(rowFilter, rows, ids, filterValue);
};

export const createTableTargetRangeFilter = rowFilter => (rows, ids, filterValue) => {
  if (!filterValue || !filterValue.length) {
    return [...rows];
  }

  return filterTableTargetRows(rowFilter, rows, ids, filterValue);
};
