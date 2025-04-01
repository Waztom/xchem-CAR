import React from 'react';
import * as yup from 'yup';
import { Formik } from 'formik';
import { Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import { DialogSection } from '../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../common/components/DialogSectionHeading';
import { FormSelectField } from '../../../common/components/FormSelectField';
import { FormTextField } from '../../../common/components/FormTextField';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';
import { CreateSubBatchSelectedTargetsList } from './components/CreateSubBatchSelectedTargetsList';
import { SubmitDialog } from '../../../common/components/SubmitDialog';
import { useCreateSubBatch } from './hooks/useCreateSubBatch';

const StyledForm = styled('form')({
  display: 'contents'
});

const FormFields = styled('div')(({ theme }) => ({
  display: 'grid',
  gap: theme.spacing(2)
}));

const possibleTags = [
  { value: 'OT-chem', text: 'OT-chem' },
  { value: 'Low-cost', text: 'Low-cost' },
  { value: 'For-review', text: 'For-review' },
  { value: 'custom', text: 'Custom tag' }
];

export const CreateSubBatchDialog = ({ open, onClose, selectedMethodsIds }) => {
  const { mutate: createSubBatch } = useCreateSubBatch({
    onSuccess: () => {
      onClose();
    },
    onError: (error) => {
      // Handle error case if needed
      console.error('Failed to create subbatch:', error);
    }
  });

  const handleSubmit = async (values, { setSubmitting }) => {
    try {
      const finalTag = values.batchtag === 'custom' ? values.customTag : values.batchtag;
      await createSubBatch({ 
        batchtag: finalTag, 
        methodids: selectedMethodsIds 
      });
    } catch (error) {
      console.error('Form submission failed:', error);
    } finally {
      setSubmitting(false);
    }
  };

  return (
    <SuspenseWithBoundary>
      <Formik
        initialValues={{
          batchtag: '',
          customTag: ''
        }}
        validationSchema={yup.object().shape({
          batchtag: yup.string().required('Required'),
          customTag: yup.string().when('batchtag', {
            is: val => val === 'custom',
            then: () => yup.string()
              .required('Required')
              .max(20, 'Maximum 20 characters'),
            otherwise: () => yup.string()
          })
        })}
        onSubmit={handleSubmit}
      >
        {({ isSubmitting, resetForm, values }) => (
          <SubmitDialog
            id="create-subbatch-dialog"
            open={open}
            title="Create subbatch"
            content={
              <StyledForm id="create-subbatch-form">
                <DialogSection>
                  <DialogSectionHeading>Batch information</DialogSectionHeading>
                  <Typography>Please provide following information:</Typography>
                  <FormFields>
                    <FormSelectField 
                      id="create-subbatch-batchtag" 
                      name="batchtag" 
                      items={possibleTags} 
                    />
                    {values.batchtag === 'custom' && (
                      <FormTextField
                        id="create-subbatch-custom-tag"
                        name="customTag"
                        label="Custom tag"
                        inputProps={{ maxLength: 20 }}
                        helperText="Maximum 20 characters"
                      />
                    )}
                  </FormFields>
                </DialogSection>

                <DialogSection>
                  <DialogSectionHeading>Targets</DialogSectionHeading>
                  <Typography>
                    These targets (with methods) will be added to the batch:
                  </Typography>
                  <SuspenseWithBoundary>
                    <CreateSubBatchSelectedTargetsList />
                  </SuspenseWithBoundary>
                </DialogSection>
              </StyledForm>
            }
            onClose={onClose}
            submitDisabled={isSubmitting}
            SubmitButtonProps={{
              type: 'submit',
              form: 'create-subbatch-form'
            }}
            TransitionProps={{
              onExited: resetForm
            }}
          />
        )}
      </Formik>
    </SuspenseWithBoundary>
  );
};

CreateSubBatchDialog.displayName = 'CreateSubBatchDialog';
