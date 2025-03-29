import { useMemo } from 'react';
import { useGetTargets } from '../../common/hooks/useGetTargets';
import { PREFERRED_LEAD_TIME, PREFERRED_PRICE, PREFERRED_VENDORS } from '../constants/preferredContstants';
import { useBatchContext } from './useBatchContext';

const getPreferredVendorsFlag = reactants => {
  const reactantVendors1 = reactants?.[0]?.catalogentries?.map(({ vendor }) => vendor) || [];
  const reactantVendors2 = reactants?.[1]?.catalogentries?.map(({ vendor }) => vendor) || [];

  const presentInBoth =
    reactantVendors1.some(vendor => PREFERRED_VENDORS.includes(vendor)) &&
    reactantVendors2.some(vendor => PREFERRED_VENDORS.includes(vendor));

  return presentInBoth;
};

const getPreferredLeadTimeFlag = reactants => {
  const reactantLeadTimes1 = reactants?.[0]?.catalogentries?.map(({ leadtime }) => leadtime) || [];
  const reactantLeadTimes2 = reactants?.[1]?.catalogentries?.map(({ leadtime }) => leadtime) || [];

  if (!reactantLeadTimes1.length || !reactantLeadTimes2.length) {
    return false;
  }

  const minimalLeadTime1 = Math.min(...reactantLeadTimes1);
  const minimalLeadTime2 = Math.min(...reactantLeadTimes2);
  const maxFromMinimals = Math.max(minimalLeadTime1, minimalLeadTime2);

  const withinRange = maxFromMinimals <= PREFERRED_LEAD_TIME;

  return withinRange;
};

const getPreferredPriceFlag = reactants => {
  const reactantPrices1 = reactants?.[0]?.catalogentries?.map(({ upperprice }) => upperprice) || [];
  const reactantPrices2 = reactants?.[1]?.catalogentries?.map(({ upperprice }) => upperprice) || [];

  if (!reactantPrices1.length || !reactantPrices2.length) {
    return false;
  }

  const minimalPrice1 = Math.min(...reactantPrices1);
  const minimalPrice2 = Math.min(...reactantPrices2);
  const sum = minimalPrice1 + minimalPrice2;

  const withinRange = sum <= PREFERRED_PRICE;

  return withinRange;
};

export const useGetTableData = () => {
  const batch = useBatchContext();

  const { data: targets } = useGetTargets({ batch_id: batch.id, fetchall: 'yes' });

  const tableData = useMemo(() => {
    if (!targets) {
      return [];
    }

    // There's a bug in react-table library which prevents subRows from being selected when parent is being selected if
    // the subRows aren't located in the subRows field.
    return targets.map(target => {
      const { methods, ...rest } = target;
      return {
        ...rest,
        subRows: methods?.map((method, index) => ({
          ...method,
          position: index + 1, // This should match with attribute method_no in exported data from ExportProjectDialog
          reactions: method.reactions?.map(reaction => ({
            ...reaction,
            // TODO this might be better to do on the BE?
            preferredVendor: getPreferredVendorsFlag(reaction.reactants),
            preferredLeadTime: getPreferredLeadTimeFlag(reaction.reactants),
            preferredPrice: getPreferredPriceFlag(reaction.reactants)
          }))
        }))
      };
    });
  }, [targets]);

  return tableData;
};
