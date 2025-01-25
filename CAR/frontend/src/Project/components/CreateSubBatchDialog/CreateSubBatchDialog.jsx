import { makeStyles, Typography } from '@material-ui/core';
import React from 'react';
import { SubmitDialog } from '../../../common/components/SubmitDialog';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';
import { useCreateSubBatch } from './hooks/useCreateSubBatch';
import { Formik, Form } from 'formik';
import * as yup from 'yup';
import { CreateSubBatchSelectedTargetsList } from './components/CreateSubBatchSelectedTargetsList';
import { DialogSection } from '../../..//common/components/DialogSection';
import { DialogSectionHeading } from '../../../common/components/DialogSectionHeading';
import { FormSelectField } from '../../../common/components/FormSelectField';

const useStyles = makeStyles(theme => ({
  form: {
    display: 'grid',
    gap: theme.spacing(2)
  }
}));

const possibleTags = [
  { value: 'OT-chem', text: 'OT-chem' },
  { value: 'Low-cost', text: 'Low-cost' },
  { value: 'For-review', text: 'For-review' }
];

export const CreateSubBatchDialog = ({ open, onClose, selectedMethodsIds }) => {
  const classes = useStyles();

  const { mutate: createSubBatch } = useCreateSubBatch();

  return (
    <Formik
      initialValues={{
        batchtag: ''
      }}
      validationSchema={yup.object().shape({
        batchtag: yup.string().required('Required')
      })}
      onSubmit={({ batchtag }) => {
        createSubBatch({ batchtag: batchtag, methodids: selectedMethodsIds });
        onClose();
      }}
    >
      {({ isSubmitting, resetForm }) => (
        <SubmitDialog
          id="create-subbatch-dialog"
          open={open}
          title="Create subbatch"
          content={
            <Form className={classes.form} id="create-subbatch-form">
              <DialogSection>
                <DialogSectionHeading>Batch information</DialogSectionHeading>
                <Typography>Please provide following information:</Typography>
                <FormSelectField id="create-subbatch-batchtag" name="batchtag" label="Name" items={possibleTags} />
              </DialogSection>

              <DialogSection>
                <DialogSectionHeading>Targets</DialogSectionHeading>
                <Typography>These targets (with methods) will be added to the batch:</Typography>
                <SuspenseWithBoundary>
                  <CreateSubBatchSelectedTargetsList />
                </SuspenseWithBoundary>
              </DialogSection>
            </Form>
          }
          onClose={onClose}
          submitDisabled={isSubmitting}
          SubmitButtonProps={{
            type: 'submit',
            form: 'create-subbatch-form'
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
