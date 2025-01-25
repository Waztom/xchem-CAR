import React from 'react';
import { useQueryClient } from 'react-query';
import { SubmitDialog } from '../../../common/components/SubmitDialog';
import { Typography } from '@material-ui/core';
import { DialogSection } from '../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../common/components/DialogSectionHeading';
import { getTargetsQueryKey } from '../../../common/api/targetsQueryKeys';
import { Form, Formik } from 'formik';
import * as yup from 'yup';
import { FormExportBatchSelector } from './components/FormExportBatchSelector';
import { parseAsync } from 'json2csv';
import { axiosGet } from '../../../common/utils/axiosFunctions';
import { useGlobalSnackbar } from '../../../common/hooks/useGlobalSnackbar';
import { useProjectSnackbar } from '../../../common/hooks/useProjectSnackbar';
import { CloseSnackbarButton } from '../../../common/components/CloseSnackbarButton';
import { DownloadButton } from './components/DownloadButton/DownloadButton';
import { getBatchesQueryKey } from '../../../common/api/batchesQueryKeys';
import { useCurrentProjectStore } from '../../../common/stores/currentProjectStore';

export const ExportProjectDialog = ({ open, onClose }) => {
  const queryClient = useQueryClient();

  const project = useCurrentProjectStore.useCurrentProject();

  const { enqueueSnackbarInfo } = useProjectSnackbar();
  const { enqueueSnackbarSuccess, enqueueSnackbarError, closeSnackbar } = useGlobalSnackbar();

  return (
    <Formik
      initialValues={{
        selectedBatchesMap: {}
      }}
      validationSchema={yup.object().shape({
        selectedBatchesMap: yup
          .object()
          .test('one-or-more-selected', 'Select at least one batch', value => Object.values(value).some(val => val))
      })}
      onSubmit={async ({ selectedBatchesMap }) => {
        onClose();

        const messageId = enqueueSnackbarInfo('Your download is being prepared...');

        try {
          // Pair the selected batch IDs with batches
          const batches = await queryClient.fetchQuery(getBatchesQueryKey({ project_id: project.id }));
          const selectedBatches = batches.filter(batch => !!selectedBatchesMap[batch.id]);

          const batchTargets = await Promise.all(
            selectedBatches.map(batch => {
              const queryKey = getTargetsQueryKey({ batch_id: batch.id, fetchall: 'yes' });

              return queryClient.fetchQuery({ queryKey, queryFn: async () => (await axiosGet(queryKey)).results });
            })
          );

          // This transforms the data into an array where an item represents a method
          const data = batchTargets
            .map((targets, i) =>
              targets.map(target =>
                target.methods?.map((method, j) => ({
                  batch_name: selectedBatches[i].batchtag,
                  target_SMILES: target.smiles,
                  method_no: j + 1, // This should match with useGetTableData position attribute for sub rows
                  ...Object.fromEntries(
                    method.reactions
                      .map((reaction, k) => [
                        ...reaction.reactants.map((reactant, l) => [`react${k + 1}.${l + 1}_SMILES`, reactant.smiles]),
                        [`prod${k + 1}_SMILES`, reaction.products[0].smiles]
                      ])
                      .flat()
                  )
                }))
              )
            )
            .flat(2);

          /**
           * Derive the columns from data. A row with the most amount of columns should have all the columns other rows have.
           * This isn't very stable because in case there is more data which varies across rows the mechanism might not
           * resolve all of the unique rows. The reason why it's done this way is to ensure the order of the columns.
           */
          const fields = data.reduce((currentBest, row) => {
            const keys = Object.keys(row);
            if (keys.length > currentBest.length) {
              return keys;
            }
            return currentBest;
          }, []);

          const csvData = await parseAsync(data, { fields });

          closeSnackbar(messageId);
          enqueueSnackbarSuccess('Your download is ready', {
            action: key => (
              <>
                <DownloadButton messageId={key} csvData={csvData} />
                <CloseSnackbarButton messageId={key} />
              </>
            )
          });
        } catch (err) {
          closeSnackbar(messageId);
          enqueueSnackbarError(err.message);
        }
      }}
    >
      {({ isSubmitting, resetForm }) => (
        <SubmitDialog
          id="export-project-dialog"
          open={open}
          title="Export project"
          content={
            <Form id="export-project-form">
              <DialogSection>
                <DialogSectionHeading>Batches</DialogSectionHeading>
                <Typography>Please select batches for export:</Typography>

                <FormExportBatchSelector name="selectedBatchesMap" label="Batch selector" />
              </DialogSection>
            </Form>
          }
          onClose={onClose}
          submitDisabled={isSubmitting}
          SubmitButtonProps={{
            type: 'submit',
            form: 'export-project-form',
            children: 'Download'
          }}
          TransitionProps={{
            onExited: () => {
              resetForm();
            }
          }}
        />
      )}
    </Formik>
  );
};
