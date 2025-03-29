import { colors, makeStyles, Tooltip } from '@material-ui/core';
import classNames from 'classnames';
import React from 'react';
import { PREFERRED_LEAD_TIME, PREFERRED_PRICE, PREFERRED_VENDORS } from '../../../../constants/preferredContstants';
import { formatPreferredVendorsString } from '../../../../utils/formatPreferredVendorsString';

const useStyles = makeStyles(theme => ({
  circle: {
    borderRadius: '50%',
    width: theme.spacing(),
    height: theme.spacing(),
    display: 'block'
  },
  green: {
    backgroundColor: colors.lime[500]
  },
  red: {
    backgroundColor: colors.red[500]
  }
}));

const typeMapping = {
  vendor: 'preferredVendor',
  leadTime: 'preferredLeadTime',
  price: 'preferredPrice'
};

const colorMapping = {
  true: 'green',
  false: 'red'
};

const tooltipMapping = {
  vendor: {
    true: `Both reaction's reactants are available from ${formatPreferredVendorsString(PREFERRED_VENDORS)}`,
    false: `Not all reactions's reactants are available from ${formatPreferredVendorsString(PREFERRED_VENDORS)}`
  },
  leadTime: {
    true: `Lead time for reaction's reactants is within ${PREFERRED_LEAD_TIME} weeks`,
    false: `Lead time for reaction's reactants exceeds ${PREFERRED_LEAD_TIME} weeks`
  },
  price: {
    true: `Price for both reaction's reactants in total is within ${PREFERRED_PRICE}`,
    false: `Price for both reaction's reactants in total exceeds ${PREFERRED_PRICE}`
  }
};

export const PreferredFlagIndicator = ({ reaction, type }) => {
  const classes = useStyles();

  const value = reaction[typeMapping[type]];
  const colorClass = classes[colorMapping[value]];

  return (
    <Tooltip title={tooltipMapping[type][value]}>
      <span className={classNames(classes.circle, colorClass)} />
    </Tooltip>
  );
};
