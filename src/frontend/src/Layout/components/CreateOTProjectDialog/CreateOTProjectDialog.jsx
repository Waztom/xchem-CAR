import React from 'react';
import { SubmitDialog } from '../../../common/components/SubmitDialog';
import { Typography } from '@mui/material';
import { useCreateOTProject } from './hooks/useCreateOTProject';
import { DialogSection } from '../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../common/components/DialogSectionHeading';
import { Form, Formik } from 'formik';
import * as yup from 'yup';
import { FormBatchSelector } from './components/FormBatchSelector';

export const CreateOTProjectDialog = ({ open, onClose }) => {
  const { mutate: createOTProject } = useCreateOTProject();

  const handleSubmit = async ({ selectedBatchesMap, startingMaterialFiles = {} }) => {
    const selectedBatchesIds = Object.entries(selectedBatchesMap)
      .filter(([_, value]) => value)
      .map(([key]) => {
        const batchId = Number(key);
        return isNaN(batchId) ? null : batchId;
      })
      .filter(Boolean);

    if (selectedBatchesIds.length > 0) {
      // Create FormData to send files
      const formData = new FormData();
      formData.append('protocol_name', `OT Protocol ${new Date().toLocaleString()}`);
      formData.append('batchids', JSON.stringify(selectedBatchesIds));
      
      // Append starting material files (if any)
      let hasFiles = false;
      
      Object.entries(startingMaterialFiles).forEach(([batchId, file]) => {
        if (file && file instanceof File) {
          hasFiles = true;
          formData.append(`starting_materials_batch_${batchId}`, file);
          console.log(`Added file for batch ${batchId}: ${file.name}`);
        }
      });
      
      formData.append('has_custom_starting_materials', hasFiles ? 'true' : 'false');

      createOTProject(formData);
      onClose();
    }
  };

  return (
    <Formik
      initialValues={{
        selectedBatchesMap: {},
        startingMaterialFiles: {}
      }}
      validationSchema={yup.object().shape({
        selectedBatchesMap: yup
          .object()
          .test('one-or-more-selected', 'Select at least one batch', value => 
            Object.values(value).some(val => val))
        // No validation for startingMaterialFiles since it's optional
      })}
      onSubmit={handleSubmit}
    >
      {({ isSubmitting, resetForm }) => (
        <SubmitDialog
          id="create-ot-protocol-dialog"
          open={open}
          title="Create OT protocol"
          content={
            <Form id="create-ot-protocol-form">
              <DialogSection>
                <DialogSectionHeading>Batches</DialogSectionHeading>
                <Typography>Please select batches for OT protocol:</Typography>
                <FormBatchSelector 
                  name="selectedBatchesMap" 
                  label="Batch selector"
                />
              </DialogSection>
            </Form>
          }
          onClose={onClose}
          submitDisabled={isSubmitting}
          SubmitButtonProps={{
            type: 'submit',
            form: 'create-ot-protocol-form'
          }}
          TransitionProps={{
            onExited: resetForm
          }}
        />
      )}
    </Formik>
  );
};
