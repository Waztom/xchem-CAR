import { useCallback, useState } from 'react';
import { Parser } from '@json2csv/plainjs';
import { saveAs } from 'file-saver';
import { getTargetsQueryKey } from '../../../../../../common/api/targetsQueryKeys';
import { axiosGet } from '../../../../../../common/utils/axiosFunctions';

/**
 * Returns a callback that fetches targets for a single batch,
 * transforms them into CSV rows, and triggers a browser download.
 */
export const useDownloadBatchCsv = (batch) => {
  const [loading, setLoading] = useState(false);

  const download = useCallback(async () => {
    if (!batch?.id) return;
    setLoading(true);

    try {
      const queryKey = getTargetsQueryKey({ batch_id: batch.id, fetchall: 'yes' });
      const targets = await axiosGet(queryKey);

      // Flatten targets → methods → reactions into one row per method
      const data = targets
        .map(target =>
          target.methods?.map((method, j) => ({
            'method-no': j + 1,
            'target-names': target.name,
            'no-steps': method.reactions.length,
            'target-SMILES': target.smiles,
            'batch-tag': batch.batchtag,
            ...Object.fromEntries(
              method.reactions
                .map((reaction, k) => [
                  ...reaction.reactants.map((reactant, l) => [
                    `reactant-${l + 1}-${k + 1}`,
                    reactant.smiles,
                  ]),
                  [`reaction-product-smiles-${k + 1}`, reaction.products[0].smiles],
                  [`reaction-name-${k + 1}`, reaction.reactionclass],
                  [`reaction-recipe-${k + 1}`, reaction.recipe],
                ])
                .flat()
            ),
          }))
        )
        .flat(2)
        .filter(Boolean);

      if (!data.length) {
        throw new Error('No data to export for this batch');
      }

      const fields = data.reduce((best, row) => {
        const keys = Object.keys(row);
        return keys.length > best.length ? keys : best;
      }, []);

      const parser = new Parser({ fields });
      const csvData = parser.parse(data);

      const blob = new Blob([csvData], { type: 'text/plain' });
      saveAs(blob, `${batch.batchtag}-${new Date().toISOString()}.csv`);
    } catch (err) {
      console.error('Batch CSV export error:', err);
    } finally {
      setLoading(false);
    }
  }, [batch]);

  return { download, loading };
};
