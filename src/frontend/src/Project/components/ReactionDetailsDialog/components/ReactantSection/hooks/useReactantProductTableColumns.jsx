import { Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import React, { useMemo } from 'react';

const StyledTypography = styled(Typography)(({ theme, centered }) => ({
  width: '100%',
  textAlign: centered ? 'center' : 'inherit'
}));

const ColumnHeader = styled(Typography)(({ theme }) => ({
  fontWeight: theme.typography.fontWeightBold
}));

export const useReactantProductTableColumns = () => {
  return useMemo(() => {
    return [
      {
        accessor: 'vendor',
        sortLabel: 'vendor',
        Header: () => (
          <ColumnHeader>Vendor</ColumnHeader>
        ),
        Cell: ({ value }) => (
          <StyledTypography component="span" noWrap>
            {value}
          </StyledTypography>
        )
      },
      {
        accessor: 'catalogid',
        sortLabel: 'catalog ID',
        Header: () => (
          <ColumnHeader>Catalog ID</ColumnHeader>
        ),
        Cell: ({ value }) => (
          <StyledTypography component="span" noWrap>
            {value}
          </StyledTypography>
        )
      },
      {
        accessor: 'upperprice',
        sortLabel: 'upper price',
        Header: () => (
          <ColumnHeader>Upper price</ColumnHeader>
        ),
        Cell: ({ value }) => (
          <StyledTypography centered component="span" noWrap>
            {value}
          </StyledTypography>
        ),
        sortType: 'number'
      },
      {
        accessor: 'priceinfo',
        sortLabel: 'price info',
        Header: () => (
          <ColumnHeader>Price info</ColumnHeader>
        ),
        Cell: ({ value }) => (
          <StyledTypography centered component="span" noWrap>
            {value}
          </StyledTypography>
        ),
        sortType: 'number'
      },
      {
        accessor: 'leadtime',
        sortLabel: 'lead time',
        Header: () => (
          <ColumnHeader>Lead time</ColumnHeader>
        ),
        Cell: ({ value }) => (
          <StyledTypography centered component="span" noWrap>
            {value}
          </StyledTypography>
        ),
        sortType: 'number'
      }
    ];
  }, []);
};
