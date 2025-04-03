import React from 'react';
import { colors, Tooltip } from '@mui/material';
import { styled } from '@mui/material/styles';
import { PREFERRED_LEAD_TIME, PREFERRED_PRICE, PREFERRED_VENDORS } from '../../../../constants/preferredContstants';
import { formatPreferredVendorsString } from '../../../../utils/formatPreferredVendorsString';

const Circle = styled('span')(({ theme, color }) => ({
  borderRadius: '50%',
  width: theme.spacing(1),
  height: theme.spacing(1),
  display: 'block',
  backgroundColor: color === 'green' ? colors.lime[500] : colors.red[500]
}));

const typeMapping = {
  vendor: 'preferredVendor',
  leadTime: 'preferredLeadTime',
  price: 'preferredPrice'
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
  const value = reaction[typeMapping[type]];
  const color = value ? 'green' : 'red';

  return (
    <Tooltip title={tooltipMapping[type][value]}>
      <Circle color={color} />
    </Tooltip>
  );
};

PreferredFlagIndicator.displayName = 'PreferredFlagIndicator';
