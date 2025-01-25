import { makeStyles, Typography } from '@material-ui/core';
import classNames from 'classnames';
import React, { useMemo } from 'react';

const useStyles = makeStyles(theme => ({
  text: {
    width: '100%'
  },
  centered: {
    textAlign: 'center'
  }
}));

export const useReactantProductTableColumns = () => {
  const classes = useStyles();

  return useMemo(() => {
    return [
      {
        accessor: 'vendor',
        sortLabel: 'vendor',
        Header: () => {
          return (
            <Typography>
              <strong>Vendor</strong>
            </Typography>
          );
        },
        Cell: ({ value }) => {
          return (
            <Typography className={classes.text} component="span" noWrap>
              {value}
            </Typography>
          );
        }
      },
      {
        accessor: 'catalogid',
        sortLabel: 'catalog ID',
        Header: () => {
          return (
            <Typography>
              <strong>Catalog ID</strong>
            </Typography>
          );
        },
        Cell: ({ value }) => {
          return (
            <Typography className={classes.text} component="span" noWrap>
              {value}
            </Typography>
          );
        }
      },
      {
        accessor: 'upperprice',
        sortLabel: 'upper price',
        Header: () => {
          return (
            <Typography>
              <strong>Upper price</strong>
            </Typography>
          );
        },
        Cell: ({ value }) => {
          return (
            <Typography className={classNames(classes.text, classes.centered)} component="span" noWrap>
              {value}
            </Typography>
          );
        },
        sortType: 'number'
      },
      {
        accessor: 'priceinfo',
        sortLabel: 'price info',
        Header: () => {
          return (
            <Typography>
              <strong>Price info</strong>
            </Typography>
          );
        },
        Cell: ({ value }) => {
          return (
            <Typography className={classNames(classes.text, classes.centered)} component="span" noWrap>
              {value}
            </Typography>
          );
        },
        sortType: 'number'
      },
      {
        accessor: 'leadtime',
        sortLabel: 'lead time',
        Header: () => {
          return (
            <Typography>
              <strong>Lead time</strong>
            </Typography>
          );
        },
        Cell: ({ value }) => {
          return (
            <Typography className={classNames(classes.text, classes.centered)} component="span" noWrap>
              {value}
            </Typography>
          );
        },
        sortType: 'number'
      }
    ];
  }, [classes.centered, classes.text]);
};
