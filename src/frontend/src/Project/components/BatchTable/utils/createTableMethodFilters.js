/**
 * Creates a table filter. If a target (parent) row is being iterated, it checks its children if at least one of them match
 * the filter and if so, returns the parent row. It repeats the check for each method row.
 */
const filterTableMethodRows = (rowFilter, rows, ids, filterValue) => {
  return (
    rows
      .filter(row => {
        if (row.depth === 0) {
          return row.subRows.some(subRow => rowFilter(subRow, ids, filterValue));
        }

        return rowFilter(row, ids, filterValue);
      })
      // A bug in the library results in loosing subRows when applying filtering. To avoid it a copy of each row
      // has to be returned from the filter method.
      .map(row => ({ ...row }))
  );
};

export const createTableMethodAutocompleteFilter = rowFilter => (rows, ids, filterValue) => {
  if (!filterValue || !filterValue.length) {
    return [...rows];
  }

  return filterTableMethodRows(rowFilter, rows, ids, filterValue);
};

export const createTableMethodYesNoFilter = rowFilter => (rows, ids, filterValue) => {
  if (filterValue === '') {
    return [...rows];
  }

  return filterTableMethodRows(rowFilter, rows, ids, filterValue);
};

export const createTableMethodRangeFilter = rowFilter => (rows, ids, filterValue) => {
  if (!filterValue || !filterValue.length) {
    return [...rows];
  }

  return filterTableMethodRows(rowFilter, rows, ids, filterValue);
};

export const createTableMethodSmilesFilter = rowFilter => (rows, ids, filterValue) => {
  if (!filterValue || !filterValue.length) {
    return [...rows];
  }

  return filterTableMethodRows(rowFilter, rows, ids, filterValue);
};
