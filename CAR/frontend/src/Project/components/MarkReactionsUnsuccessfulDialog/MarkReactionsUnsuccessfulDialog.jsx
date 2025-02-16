import React from 'react';
import { Form, Formik } from 'formik';
import * as yup from 'yup';
import { SubmitDialog } from '../../../common/components/SubmitDialog';
import { DialogSection } from '../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../common/components/DialogSectionHeading';
import { Typography } from '@material-ui/core';
import { FormReactionIdsSelector } from './components/FormReactionIdsSelector';
import { useGetReactionsData } from './hooks/useGetReactionsData';
import { useMarkReactionsUnsuccessful } from './hooks/useMarkReactionsUnsuccessful';

export const MarkReactionsUnsuccessfulDialog = ({ open, onClose }) => {
  const { reactionsIds, reactionsMap } = useGetReactionsData();

  const { mutate: markReactionsUnsuccessful } = useMarkReactionsUnsuccessful();

  return (
    <Formik
      initialValues={{
        reactionIds: []
      }}
      validationSchema={yup.object().shape({
        reactionIds: yup
          .array(yup.string())
          .min(1, 'At least one reaction ID has to be selected')
          .test('existing-ids', 'Selection contains invalid or nonexistent reaction IDs', value =>
            value.every(id => !!reactionsMap[id])
          )
      })}
      onSubmit={({ reactionIds }) => {
        markReactionsUnsuccessful({ reaction_ids: reactionIds.map(reactionId => Number(reactionId)) });
        onClose();
      }}
    >
      {({ isSubmitting, resetForm }) => (
        <SubmitDialog
          id="mark-reactions-unsuccessful-dialog"
          open={open}
          title="Mark reactions unsuccessful"
          content={
            <Form id="mark-reactions-unsuccessful-form">
              <DialogSection>
                <DialogSectionHeading>Reactions</DialogSectionHeading>
                <Typography>Please select reactions IDs or upload them from a CSV file:</Typography>

                <FormReactionIdsSelector
                  name="reactionIds"
                  label="Reaction IDs"
                  autocompleteId="mark-reactions-unsuccessful-autocomplete"
                  filePickerId="mark-reactions-unsuccessful-file"
                  reactionsIds={reactionsIds}
                  reactionsMap={reactionsMap}
                />
              </DialogSection>
            </Form>
          }
          onClose={onClose}
          submitDisabled={isSubmitting}
          SubmitButtonProps={{
            type: 'submit',
            form: 'mark-reactions-unsuccessful-form'
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
